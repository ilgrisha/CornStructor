# File: backend/app/core/assembly/ga_overlap_selector.py
# Version: v1.6.0

"""
Genetic Algorithm to select orthogonal, specific overlaps between neighboring oligo bodies.

What's new in v1.6.0 — Early-stop & run-summary of stagnation
- **Optional early-stop** when best fitness does not improve for a long window
  (default: `3 × stagnation_patience`). Tunable via `GAParameters.enable_early_stop`
  and `GAParameters.early_stop_multiplier`.
- **Run summary at INFO** after `evolve()` finishes: reports generations run,
  best fitness, early-stop flag, count of boosts/restarts, total boosted gens,
  and last-generation diversity.
- **Stagnation events surfaced**: whenever a stagnation boost or partial restart
  occurs, it is recorded and listed in the run summary.
- Keeps v1.5.x anti-stagnation (mutation & immigrant boost, optional partial restart)
  and v1.4.x diagnostics (rich per-gen stats) and logging of resolved runtime knobs.
"""

from __future__ import annotations

import logging
import multiprocessing
import random
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, fields
from statistics import mean, median, pstdev
from typing import Dict, List, Optional, Tuple

from .batch_fitness import (
    BatchOptions,
    pairwise_norm_edit_stats,
    rc_kmer_hits_total,
)
from .fitness_utils import (
    count_3prime_mismatch_events,
    count_motif_occurrences,
    gc_fraction,  # returns GC **percent** [0..100]
    contains_any_motif,
    max_same_base_run,
    reverse_complement,
    compute_tm,  # PRIMER3 by default (falls back automatically)
)

logger = logging.getLogger(__name__)

# === Fitness weights (can be overridden) ===
ORTHO_AVG_WEIGHT = 100.0
ORTHO_WORST_WEIGHT = 200.0
MIS3_PENALTY = 2.0
MISK_PENALTY = 1.0
OVERLAP_REWARD = 10.0
ALLOWED_MOTIF_SCALE = 1.0

# === GA defaults ===
POPULATION_SIZE = 200
NUM_GENERATIONS = 3
MUTATION_RATE_INITIAL = 0.35
CROSSOVER_RATE = 0.6
ELITISM_COUNT = 1
RANDOM_INJECTION_RATE = 0.05
TOURNAMENT_SIZE = 3

# Local search (kept off by default generations < period)
ENABLE_HILL_CLIMB = True
HILL_CLIMB_PERIOD = 5
HILL_CLIMB_JUNC = 2
HILL_CLIMB_CANDS = 10


# ---------------------------------------------------------------------
# Dataclasses / Errors
# ---------------------------------------------------------------------

@dataclass
class GAParameters:
    """
    GA knobs, fitness scales, perf + anti-stagnation / early-stop options.
    """
    # Population / operators
    population_size: int = POPULATION_SIZE
    num_generations: int = NUM_GENERATIONS
    mutation_rate_initial: float = MUTATION_RATE_INITIAL
    crossover_rate: float = CROSSOVER_RATE
    elitism_count: int = ELITISM_COUNT
    random_injection_rate: float = RANDOM_INJECTION_RATE
    tournament_size: int = TOURNAMENT_SIZE

    # Local search
    enable_hill_climb: bool = ENABLE_HILL_CLIMB
    hill_climb_period: int = HILL_CLIMB_PERIOD
    hill_climb_junc: int = HILL_CLIMB_JUNC
    hill_climb_cands: int = HILL_CLIMB_CANDS

    # Fitness extras
    allowed_motif_scale: float = ALLOWED_MOTIF_SCALE

    # Performance
    workers: int = 0                # 0 → auto (≈0.9× CPUs)
    batch_chunk_size: int = 128

    # Anti-stagnation
    stagnation_patience: int = 5            # gens with no improvement before action
    stagnation_improve_epsilon: float = 1e-6
    stagnation_mutation_mult: float = 2.0   # × current mutation when firing
    mutation_rate_cap: float = 0.9
    stagnation_inject_add: float = 0.15     # add to random injection during boost
    stagnation_boost_gens: int = 2          # gens to keep boost active
    enable_restart: bool = True
    stagnation_restart_fraction: float = 0.5  # replace bottom fraction when firing
    stagnation_max_restarts: int = 3
    diversity_min_unique_pct: float = 20.0  # preemptive boost if unique% < threshold

    # Early-stop (new)
    enable_early_stop: bool = True
    early_stop_multiplier: float = 3.0      # stop after multiplier × stagnation_patience gens


@dataclass
class RunSummary:
    """
    Programmatic summary of the last GA run (also logged at INFO).
    """
    generations_run: int
    early_stopped: bool
    early_stop_gen: Optional[int]
    best_fitness: float
    boosts: int
    restarts: int
    boost_gens_total: int
    unique_pct_last: float
    boost_gens: List[int]
    restart_gens: List[int]


def _coerce_ga_params(obj: Optional[object]) -> GAParameters:
    """
    Accept GAParameters or any GA-like object (e.g., GAParamsView) and return a GAParameters
    with known fields copied and missing fields defaulted.
    """
    if obj is None or isinstance(obj, GAParameters):
        return obj or GAParameters()
    coerced = GAParameters()
    for f in fields(GAParameters):
        if hasattr(obj, f.name):
            try:
                setattr(coerced, f.name, getattr(obj, f.name))
            except Exception:
                pass
    return coerced


class NoOverlapCandidatesError(RuntimeError):
    """Raised when a junction has zero viable overlap candidates after filtering."""

    def __init__(self, junction_index: int, left_end: int, right_start: int, reasons: Dict[str, int]):
        self.junction_index = junction_index
        self.left_end = left_end
        self.right_start = right_start
        self.reasons = dict(sorted(reasons.items(), key=lambda kv: kv[1], reverse=True))
        super().__init__(self.formatted_reasons())

    def formatted_reasons(self) -> str:
        if not self.reasons:
            return "no candidates (no reason counters recorded)"
        items = [f"{k}:{v}" for k, v in self.reasons.items()]
        return ", ".join(items)


class OverlapChromosome:
    """
    Candidate solution: list of overlaps per junction. Each overlap is (sequence, abs_start, abs_end).
    """
    def __init__(self, overlaps: List[Tuple[str, int, int]]):
        self.overlaps = overlaps
        self.fitness: Optional[float] = None


# ---------------------------------------------------------------------
# Core GA
# ---------------------------------------------------------------------

class GAOverlapSelector:
    """
    Optimize overlaps at each junction via a genetic algorithm.

    Behavior summary:
    - Generates overlap candidates per junction with biophysical filters.
    - Evolves with tournament selection, single-point crossover, mutation and optional hill-climb.
    - Uses thread-pooled evaluation leveraging C-accelerated hot paths.
    - Detects stagnation (no Δbest > epsilon for `stagnation_patience` gens) and triggers:
        * temporary mutation & immigrant boost (for `stagnation_boost_gens`)
        * optional partial restart (keep top elites; refill remaining population randomly)
    - **Early-stop**: optional termination when no improvement persists for
      `stagnation_patience * early_stop_multiplier` generations.
    """

    # Public: after `evolve()`, this holds a RunSummary
    last_run_summary: Optional[RunSummary] = None

    def __init__(
        self,
        full_sequence: str,
        oligo_positions: List[Tuple[int, int]],
        oligo_seqs: List[str],
        overlap_min_size: int,
        overlap_max_size: int,
        tm_min: Optional[float] = None,
        tm_max: Optional[float] = None,
        gc_min: Optional[float] = None,  # percent [0..100]
        gc_max: Optional[float] = None,  # percent [0..100]
        run_min: Optional[int] = None,
        run_max: Optional[int] = None,
        disallowed_motifs: Optional[List[str]] = None,
        allowed_motifs: Optional[Dict[str, float]] = None,
        max_gap_allowed: int = 0,
        enforce_span: bool = False,
        tm_method: str = "PRIMER3",
        tm_params: Optional[Dict[str, float]] = None,
        ga_params: Optional[object] = None,
    ):
        self.full_sequence = full_sequence
        self.oligo_positions = oligo_positions
        self.oligo_seqs = oligo_seqs
        self.num_overlaps = len(oligo_positions) - 1
        self.overlap_min_size = overlap_min_size
        self.overlap_max_size = overlap_max_size

        # Filters / preferences
        self.tm_min = tm_min
        self.tm_max = tm_max
        self.gc_min = gc_min
        self.gc_max = gc_max
        self.run_min = run_min
        self.run_max = run_max
        self.disallowed_motifs = [m.upper() for m in (disallowed_motifs or [])]
        self.allowed_motifs = {m.upper(): float(w) for m, w in (allowed_motifs or {}).items()}

        # Thermo
        self.tm_method = tm_method
        self.tm_params = tm_params or {}

        # GA knobs (coerced to our dataclass with defaults)
        self.ga: GAParameters = _coerce_ga_params(ga_params)
        self.mutation_rate = self.ga.mutation_rate_initial

        # Workers
        auto_workers = max(1, int(multiprocessing.cpu_count() * 0.9))
        self.workers = self.ga.workers or auto_workers
        self.batch_opts = BatchOptions(workers=self.workers, chunk_size=self.ga.batch_chunk_size)

        # Anti-stagnation runtime state
        self._no_improve = 0              # recent (since last improvement/restart)
        self._no_improve_total = 0        # continuous total for early-stop
        self._restart_count = 0
        self._boost_until_gen = 0         # inclusive generation index when boost ends
        self._boost_events: List[int] = []
        self._restart_events: List[int] = []

        # Log resolved runtime knobs
        logger.info(
            "GAOverlapSelector init: workers=%d, batch_chunk_size=%d, pop=%d, gens=%d, tm=%s",
            self.workers, self.ga.batch_chunk_size, self.ga.population_size,
            self.ga.num_generations, self.tm_method,
        )

        # Precompute valid candidates per junction
        self.valid_overlap_sets: List[List[Tuple[str, int, int]]] = self._build_all_candidates()

    # ---------------- Utility / Setup ----------------

    def set_cpu_fraction(self, frac: float):
        """Limit worker threads to a fraction of available CPUs (0<frac≤1)."""
        frac = max(0.05, min(1.0, float(frac)))
        self.workers = max(1, int(multiprocessing.cpu_count() * frac))
        self.batch_opts = BatchOptions(workers=self.workers, chunk_size=self.ga.batch_chunk_size)
        logger.info("GAOverlapSelector set_cpu_fraction: workers=%d (frac=%.2f)", self.workers, frac)

    def _build_all_candidates(self) -> List[List[Tuple[str, int, int]]]:
        """Generate and filter candidates around each junction."""
        N = len(self.full_sequence)
        all_sets: List[List[Tuple[str, int, int]]] = []

        for j in range(self.num_overlaps):
            left_s, left_e = self.oligo_positions[j]
            right_s, right_e = self.oligo_positions[j + 1]

            ws = max(0, left_e - self.overlap_max_size)
            we = min(N, right_s + self.overlap_max_size)
            window = self.full_sequence[ws:we]

            reasons: Dict[str, int] = {
                "length_window": 0, "tm_low": 0, "tm_high": 0,
                "gc_low": 0, "gc_high": 0, "run_low": 0, "run_high": 0,
                "motif_disallowed": 0,
            }

            cands: List[Tuple[str, int, int]] = []
            for L in range(self.overlap_min_size, self.overlap_max_size + 1):
                if L <= 0 or L > len(window):
                    continue
                for k in range(0, len(window) - L + 1):
                    abs_s = ws + k
                    abs_e = abs_s + L
                    seq = self.full_sequence[abs_s:abs_e]

                    # Tm
                    if (self.tm_min is not None) or (self.tm_max is not None):
                        t = compute_tm(seq, method=self.tm_method, **self.tm_params)
                        if (self.tm_min is not None) and (t < self.tm_min):
                            reasons["tm_low"] += 1; continue
                        if (self.tm_max is not None) and (t > self.tm_max):
                            reasons["tm_high"] += 1; continue

                    # GC %
                    if (self.gc_min is not None) or (self.gc_max is not None):
                        gcp = gc_fraction(seq)
                        if (self.gc_min is not None) and (gcp < self.gc_min):
                            reasons["gc_low"] += 1; continue
                        if (self.gc_max is not None) and (gcp > self.gc_max):
                            reasons["gc_high"] += 1; continue

                    # Homopolymers
                    if (self.run_min is not None) or (self.run_max is not None):
                        r = max_same_base_run(seq)
                        if (self.run_min is not None) and (r < self.run_min):
                            reasons["run_low"] += 1; continue
                        if (self.run_max is not None) and (r > self.run_max):
                            reasons["run_high"] += 1; continue

                    # Disallowed motifs
                    if self.disallowed_motifs and contains_any_motif(seq, self.disallowed_motifs):
                        reasons["motif_disallowed"] += 1; continue

                    cands.append((seq, abs_s, abs_e))

            if not cands:
                raise NoOverlapCandidatesError(j, left_e, right_s, reasons)

            all_sets.append(cands)

        return all_sets

    # ---------------- GA primitives ----------------

    def overlaps_disjoint(self, overlaps: List[Tuple[str, int, int]]) -> bool:
        for i in range(len(overlaps) - 1):
            _, _, e1 = overlaps[i]
            _, s2, _ = overlaps[i + 1]
            if e1 > s2:
                return False
        return True

    def initial_population(self) -> List[OverlapChromosome]:
        pop: List[OverlapChromosome] = []
        greedy = [max(cands, key=lambda x: len(x[0])) for cands in self.valid_overlap_sets]
        pop.append(OverlapChromosome(greedy))
        while len(pop) < self.ga.population_size:
            ovrs = [random.choice(cands) for cands in self.valid_overlap_sets]
            if self.overlaps_disjoint(ovrs):
                pop.append(OverlapChromosome(ovrs))
        return pop

    def _build_full_oligos(self, chromo: OverlapChromosome) -> List[str]:
        bodies = self.oligo_seqs
        positions = self.oligo_positions
        ovs = chromo.overlaps
        full_oligos: List[str] = []
        for idx, ((s0, e0), _) in enumerate(zip(positions, bodies)):
            prev = ovs[idx - 1] if idx > 0 else None
            nxt = ovs[idx] if idx < len(bodies) - 1 else None
            start = prev[1] if prev else s0
            end = nxt[2] if nxt else e0
            frag = self.full_sequence[start:end]
            seq = frag if idx % 2 == 0 else reverse_complement(frag)
            full_oligos.append(seq)
        return full_oligos

    def _evaluate_fitness(self, chromo: OverlapChromosome) -> float:
        full_oligos = self._build_full_oligos(chromo)
        overlap_seqs = [ov[0] for ov in chromo.overlaps]

        avg_nd, min_nd = pairwise_norm_edit_stats(overlap_seqs, options=self.batch_opts)
        score = ORTHO_AVG_WEIGHT * avg_nd + ORTHO_WORST_WEIGHT * min_nd

        total_events3 = 0
        total_hitsk = 0
        for ov in overlap_seqs:
            rc_ov = reverse_complement(ov)
            total_events3 += count_3prime_mismatch_events(ov, full_oligos)
            total_events3 += count_3prime_mismatch_events(rc_ov, full_oligos)
            total_hitsk   += rc_kmer_hits_total(ov, full_oligos, k=10, options=self.batch_opts)
            total_hitsk   += rc_kmer_hits_total(rc_ov, full_oligos, k=10, options=self.batch_opts)

        score -= MIS3_PENALTY * total_events3
        score -= MISK_PENALTY * total_hitsk
        score += OVERLAP_REWARD * len(overlap_seqs)

        if self.allowed_motifs:
            motif_reward = 0.0
            for ov in overlap_seqs:
                for m, w in self.allowed_motifs.items():
                    if not m:
                        continue
                    motif_reward += w * count_motif_occurrences(ov, m)
            score += self.ga.allowed_motif_scale * motif_reward

        chromo.fitness = score
        return score

    def tournament(self, population: List[OverlapChromosome]) -> OverlapChromosome:
        k = max(2, min(self.ga.tournament_size, len(population)))
        cohort = random.sample(population, k)
        return max(cohort, key=lambda c: c.fitness or -1e9)

    def select_parents(self, population: List[OverlapChromosome]) -> Tuple[OverlapChromosome, OverlapChromosome]:
        return self.tournament(population), self.tournament(population)

    def crossover(self, p1: OverlapChromosome, p2: OverlapChromosome) -> OverlapChromosome:
        if self.num_overlaps < 2 or random.random() > self.ga.crossover_rate:
            return OverlapChromosome(p1.overlaps[:])
        pt = random.randint(1, self.num_overlaps - 1)
        child_ov = p1.overlaps[:pt] + p2.overlaps[pt:]
        if not self.overlaps_disjoint(child_ov):
            for i in range(len(child_ov) - 1):
                _, _, e1 = child_ov[i]
                _, s2, _ = child_ov[i + 1]
                if e1 > s2:
                    child_ov[i]   = random.choice(self.valid_overlap_sets[i])
                    child_ov[i+1] = random.choice(self.valid_overlap_sets[i + 1])
        if not self.overlaps_disjoint(child_ov):
            return OverlapChromosome(p1.overlaps[:])
        return OverlapChromosome(child_ov)

    def mutate(self, chromo: OverlapChromosome):
        for i in range(self.num_overlaps):
            if random.random() < self.mutation_rate:
                orig = chromo.overlaps[i]
                for _ in range(10):
                    cand = random.choice(self.valid_overlap_sets[i])
                    chromo.overlaps[i] = cand
                    if self.overlaps_disjoint(chromo.overlaps):
                        break
                else:
                    chromo.overlaps[i] = orig

    # ---------------- Main loop ----------------

    def evolve(self) -> Tuple[OverlapChromosome, List[float]]:
        """
        Run the GA and return the best chromosome and a per-generation best-fitness log.
        Also stores a `RunSummary` in `self.last_run_summary` and logs it at INFO.
        """
        population = self.initial_population()
        progress_log: List[float] = []
        prev_best = float("-inf")
        gens = max(1, int(self.ga.num_generations))
        self.mutation_rate = self.ga.mutation_rate_initial

        # Counters for run summary
        boosts = 0
        boost_gens_total = 0
        early_stopped = False
        early_stop_gen: Optional[int] = None

        # Early-stop patience in generations
        early_patience = int(max(1, round(self.ga.stagnation_patience * self.ga.early_stop_multiplier)))

        for gen in range(1, gens + 1):
            # Evaluate fitness (threads)
            with ThreadPoolExecutor(max_workers=self.workers) as exe:
                futures = {exe.submit(self._evaluate_fitness, c): c for c in population}
                for fut in as_completed(futures):
                    futures[fut].fitness = fut.result()

            # Aggregate stats
            fvals = [c.fitness if c.fitness is not None else float("-inf") for c in population]
            pop_size = len(population)
            best = max(population, key=lambda c: c.fitness or -1e9)
            best_f = best.fitness or float("-inf")
            worst_f = min(fvals)
            mean_f = mean(fvals) if pop_size else float("nan")
            median_f = median(fvals) if pop_size else float("nan")
            std_f = pstdev(fvals) if pop_size > 1 else 0.0
            unique_sets = len({tuple(ov[0] for ov in c.overlaps) for c in population})
            uniq_pct = 100.0 * unique_sets / max(1, pop_size)
            impr = (best_f - prev_best) if prev_best != float("-inf") else 0.0
            progress_log.append(best_f)

            # Log generation summary
            logger.info(
                "Gen%02d: size=%d uniq=%d (%.1f%%) best=%.3f mean=%.3f median=%.3f worst=%.3f std=%.3f mut=%.2f Δbest=%.3f",
                gen, pop_size, unique_sets, uniq_pct, best_f, mean_f, median_f, worst_f, std_f,
                self.mutation_rate, impr,
            )

            # --- Stagnation detection & handling ---
            if best_f > prev_best + self.ga.stagnation_improve_epsilon:
                self._no_improve = 0
                self._no_improve_total = 0
                prev_best = best_f
                # gently cool mutation
                self.mutation_rate = max(self.ga.mutation_rate_initial, self.mutation_rate * 0.95)
            else:
                self._no_improve += 1
                self._no_improve_total += 1
                diversity_collapsed = (uniq_pct < self.ga.diversity_min_unique_pct)
                if self._no_improve >= self.ga.stagnation_patience or diversity_collapsed:
                    # temporary boosts
                    self.mutation_rate = min(self.ga.mutation_rate_cap,
                                             self.mutation_rate * self.ga.stagnation_mutation_mult)
                    self._boost_until_gen = max(self._boost_until_gen, gen + self.ga.stagnation_boost_gens)
                    boosts += 1
                    self._boost_events.append(gen)
                    logger.info(
                        "  stagnation detected (no_improve=%d, uniq=%.1f%%) → boost: mut=%.2f; boost_until_gen=%d",
                        self._no_improve, uniq_pct, self.mutation_rate, self._boost_until_gen,
                    )
                    # optional partial restart
                    if self.ga.enable_restart and self._restart_count < self.ga.stagnation_max_restarts:
                        keep = max(1, self.ga.elitism_count)
                        target = max(keep, int(self.ga.population_size * (1.0 - self.ga.stagnation_restart_fraction)))
                        elites = sorted(population, key=lambda c: c.fitness or -1e9, reverse=True)[:keep]
                        new_pop = elites[:]
                        while len(new_pop) < target:
                            ovrs = [random.choice(c) for c in self.valid_overlap_sets]
                            if self.overlaps_disjoint(ovrs):
                                new_pop.append(OverlapChromosome(ovrs))
                        population = new_pop
                        self._restart_count += 1
                        self._restart_events.append(gen)
                        self._no_improve = 0            # reset recent patience after restart
                        # keep total counter for early-stop

                        logger.info(
                            "  partial restart: kept=%d, rebuilt_to=%d, restarts=%d",
                            keep, len(population), self._restart_count
                        )

            # Count boosted gens
            if gen <= self._boost_until_gen:
                boost_gens_total += 1

            # --- Early-stop check ---
            if self.ga.enable_early_stop and self._no_improve_total >= early_patience:
                early_stopped = True
                early_stop_gen = gen
                logger.info(
                    "Early-stop triggered at Gen%02d (no improvement for %d consecutive gens ≥ patience %d).",
                    gen, self._no_improve_total, early_patience
                )
                # finalize with current population
                final_best = max(population, key=lambda c: c.fitness or -1e9)
                # Build next generation skipped; break the loop
                # Prepare last-summary stats for log at the end (after loop)
                population = sorted(population, key=lambda c: c.fitness or -1e9, reverse=True)
                # cut progress_log to actual gens run (already correct)
                break

            # Optional hill-climb
            if self.ga.enable_hill_climb and gen % self.ga.hill_climb_period == 0:
                pre = best.fitness or float("-inf")
                junctions = random.sample(range(self.num_overlaps), min(self.ga.hill_climb_junc, self.num_overlaps))
                for j in junctions:
                    cands = random.sample(self.valid_overlap_sets[j],
                                          min(self.ga.hill_climb_cands, len(self.valid_overlap_sets[j])))
                    for cand in cands:
                        trial = OverlapChromosome(best.overlaps[:])
                        trial.overlaps[j] = cand
                        if self.overlaps_disjoint(trial.overlaps):
                            fval = self._evaluate_fitness(trial)
                            if fval > best.fitness:
                                best = trial
                post = best.fitness or float("-inf")
                if post > pre:
                    logger.info("  hill-climb improved best: +%.3f → %.3f", post - pre, post)
                    prev_best = max(prev_best, post)
                    self._no_improve = 0
                    self._no_improve_total = 0

            # --- Build next generation ---
            base_inject = self.ga.random_injection_rate
            effective_inject = base_inject
            if gen <= self._boost_until_gen:
                effective_inject = min(0.95, base_inject + self.ga.stagnation_inject_add)

            sorted_pop = sorted(population, key=lambda c: c.fitness or -1e9, reverse=True)
            new_pop = sorted_pop[: self.ga.elitism_count]

            while len(new_pop) < self.ga.population_size:
                if random.random() < effective_inject:
                    ovrs = [random.choice(c) for c in self.valid_overlap_sets]
                    if self.overlaps_disjoint(ovrs):
                        new_pop.append(OverlapChromosome(ovrs))
                else:
                    p1, p2 = self.select_parents(population)
                    child = self.crossover(p1, p2)
                    self.mutate(child)
                    new_pop.append(child)

            population = new_pop

        # ---- Finalization & run summary ----
        final_best = max(population, key=lambda c: c.fitness or -1e9)
        last_fvals = [c.fitness if c.fitness is not None else float("-inf") for c in population]
        last_unique_sets = len({tuple(ov[0] for ov in c.overlaps) for c in population})
        last_uniq_pct = 100.0 * last_unique_sets / max(1, len(population))

        gens_run = len(progress_log)

        self.last_run_summary = RunSummary(
            generations_run=gens_run,
            early_stopped=early_stopped,
            early_stop_gen=early_stop_gen,
            best_fitness=final_best.fitness or float("-inf"),
            boosts=len(self._boost_events),
            restarts=self._restart_count,
            boost_gens_total=boost_gens_total,
            unique_pct_last=last_uniq_pct,
            boost_gens=self._boost_events[:],
            restart_gens=self._restart_events[:],
        )

        logger.info(
            "GA run summary: gens=%d, early_stop=%s%s, best=%.3f, boosts=%d (boost_gens=%d @%s), "
            "restarts=%d @%s, uniq_last=%.1f%%, last_pop_size=%d",
            gens_run,
            str(early_stopped),
            f"@Gen{early_stop_gen:02d}" if early_stopped and early_stop_gen is not None else "",
            self.last_run_summary.best_fitness,
            self.last_run_summary.boosts,
            self.last_run_summary.boost_gens_total,
            self.last_run_summary.boost_gens,
            self.last_run_summary.restarts,
            self.last_run_summary.restart_gens,
            last_uniq_pct,
            len(population),
        )

        return final_best, progress_log
