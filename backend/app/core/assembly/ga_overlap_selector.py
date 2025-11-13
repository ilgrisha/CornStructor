# File: backend/app/core/assembly/ga_overlap_selector.py
# Version: v1.8.0

"""
Genetic Algorithm to select orthogonal, specific overlaps between neighboring oligo bodies.

v1.8.0 — Light Tm shaping in fitness
- Slight reward for higher average overlap Tm.
- Very slight reward for tighter (smaller) Tm spread (std).

v1.7.9 — Junction-straddling enforcement + logging fix
- Enforce that every overlap **straddles** its junction: for junction between
  left body [left_s,left_e) and right body [right_s,right_e), a candidate (abs_s,abs_e)
  must satisfy: abs_s < left_e AND abs_e > right_s. This prevents the GA from
  pulling overlaps entirely inside a body, which could shrink a child's full length
  below min_children_size.
- Fix GA run-summary logger format: print event lists with %s (not %d).

v1.7.8 — Per-junction minimum overlap enforcement
- Optional `min_overlap_per_junction: List[int]` sets L_j ≥ min_j at each junction.

v1.7.6–1.7.7 — Minor fixes and exports
"""

from __future__ import annotations

__all__ = ["GAOverlapSelector", "NoOverlapCandidatesError"]

import logging
import multiprocessing
import random
random.seed(7)
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
    gc_fraction,          # returns GC **percent** [0..100]
    contains_any_motif,
    max_same_base_run,
    reverse_complement,
    compute_tm,           # PRIMER3 by default (falls back automatically)
)

logger = logging.getLogger(__name__)

# === Fitness weights ==========================================================
ORTHO_AVG_WEIGHT = 100.0
ORTHO_WORST_WEIGHT = 200.0
MIS3_PENALTY = 2.0
MISK_PENALTY = 1.0
OVERLAP_REWARD = 10.0
ALLOWED_MOTIF_SCALE = 1.0

# Prefer shorter overlaps slightly (within constraints)
OVERLAP_LENGTH_PENALTY = 0.8

# NEW (v1.8.0): gentle Tm shaping
# - Average overlap Tm gets a slight bonus.
# - Tm std gets a very slight penalty (so lower std => higher score).
TM_AVG_BONUS = 0.25
TM_STD_PENALTY = 0.05

# === GA defaults ==============================================================
POPULATION_SIZE = 200
NUM_GENERATIONS = 3
MUTATION_RATE_INITIAL = 0.35
CROSSOVER_RATE = 0.6
ELITISM_COUNT = 1
RANDOM_INJECTION_RATE = 0.05
TOURNAMENT_SIZE = 3

# Local search
ENABLE_HILL_CLIMB = True
HILL_CLIMB_PERIOD = 5
HILL_CLIMB_JUNC = 2
HILL_CLIMB_CANDS = 10

# Very-hard penalties for full-length violations
VERY_HARD_PENALTY = 1e9
HARD_SCALE_PENALTY = 1e6


# ---------------------------------------------------------------------
# Dataclasses / Errors
# ---------------------------------------------------------------------

@dataclass
class GAParameters:
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
    workers: int = 0
    batch_chunk_size: int = 128

    # Anti-stagnation
    stagnation_patience: int = 5
    stagnation_improve_epsilon: float = 1e-6
    stagnation_mutation_mult: float = 2.0
    mutation_rate_cap: float = 0.9
    stagnation_inject_add: float = 0.15
    stagnation_boost_gens: int = 2
    enable_restart: bool = True
    stagnation_restart_fraction: float = 0.5
    stagnation_max_restarts: int = 3
    diversity_min_unique_pct: float = 20.0

    # Early-stop
    enable_early_stop: bool = True
    early_stop_multiplier: float = 3.0


@dataclass
class RunSummary:
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
    """list of overlaps per junction. Each overlap is (sequence, abs_start, abs_end)."""
    def __init__(self, overlaps: List[Tuple[str, int, int]]):
        self.overlaps = overlaps
        self.fitness: Optional[float] = None


# ---------------------------------------------------------------------
# Core GA
# ---------------------------------------------------------------------

class GAOverlapSelector:
    """Optimize overlaps at each junction via a genetic algorithm."""
    last_run_summary: Optional[RunSummary] = None  # filled after evolve()

    def __init__(
        self,
        full_sequence: str,
        oligo_positions: List[Tuple[int, int]],
        oligo_seqs: List[str],
        overlap_min_size: int,
        overlap_max_size: int,
        tm_min: Optional[float] = None,
        tm_max: Optional[float] = None,
        gc_min: Optional[float] = None,
        gc_max: Optional[float] = None,
        run_min: Optional[int] = None,
        run_max: Optional[int] = None,
        disallowed_motifs: Optional[List[str]] = None,
        allowed_motifs: Optional[Dict[str, float]] = None,
        max_gap_allowed: int = 0,
        enforce_span: bool = False,
        tm_method: str = "PRIMER3",
        tm_params: Optional[Dict[str, float]] = None,
        ga_params: Optional[object] = None,
        # Child FULL-length constraints
        min_child_full_len: Optional[int] = None,
        max_child_full_len: Optional[int] = None,
        # Per-junction min overlap lengths
        min_overlap_per_junction: Optional[List[int]] = None,
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

        # GA knobs
        self.ga: GAParameters = _coerce_ga_params(ga_params)
        self.mutation_rate = self.ga.mutation_rate_initial

        # Workers
        auto_workers = max(1, int(multiprocessing.cpu_count() * 0.9))
        self.workers = self.ga.workers or auto_workers
        self.batch_opts = BatchOptions(workers=self.workers, chunk_size=self.ga.batch_chunk_size)

        # Anti-stagnation runtime state
        self._no_improve = 0
        self._no_improve_total = 0
        self._restart_count = 0
        self._boost_until_gen = 0
        self._boost_events: List[int] = []
        self._restart_events: List[int] = []

        # Per-generation progress details
        self.progress_detail: List[Dict[str, float]] = []

        # Child FULL-length constraints
        self.child_len_min = min_child_full_len
        self.child_len_max = max_child_full_len

        # Per-junction min overlaps
        if min_overlap_per_junction is not None and len(min_overlap_per_junction) != self.num_overlaps:
            raise ValueError("min_overlap_per_junction length must equal number of junctions")
        self.min_overlap_per_junction = list(min_overlap_per_junction) if min_overlap_per_junction else None

        # Best evaluated chromosome across the run (global best)
        self._best_overall: Optional[OverlapChromosome] = None

        logger.info(
            "GAOverlapSelector init: workers=%d, chunk=%d, pop=%d, gens=%d, tm=%s, child_len=[%s..%s]",
            self.workers, self.ga.batch_chunk_size, self.ga.population_size,
            self.ga.num_generations, self.tm_method,
            str(self.child_len_min), str(self.child_len_max),
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
        """
        Generate and filter candidates around each junction.
        Enforce:
          • per-junction minimum length (if provided),
          • overlap must **straddle** the junction (abs_s < left_e AND abs_e > right_s),
          • Tm/GC/run/motif filters.
        """
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
                "motif_disallowed": 0, "span_fail": 0,
            }

            # Per-junction minimum L
            local_min = self.overlap_min_size
            if self.min_overlap_per_junction is not None:
                local_min = max(local_min, int(self.min_overlap_per_junction[j]))

            cands: List[Tuple[str, int, int]] = []
            for L in range(local_min, self.overlap_max_size + 1):
                if L <= 0 or L > len(window):
                    continue
                for k in range(0, len(window) - L + 1):
                    abs_s = ws + k
                    abs_e = abs_s + L

                    # Must straddle the junction (bridge across the cut)
                    if not (abs_s < left_e and abs_e > right_s):
                        reasons["span_fail"] += 1
                        continue

                    seq = self.full_sequence[abs_s:abs_e]

                    # Tm check
                    if (self.tm_min is not None) or (self.tm_max is not None):
                        t = compute_tm(seq, method=self.tm_method, **self.tm_params)
                        if (self.tm_min is not None) and (t < self.tm_min):
                            reasons["tm_low"] += 1; continue
                        if (self.tm_max is not None) and (t > self.tm_max):
                            reasons["tm_high"] += 1; continue

                    # GC%
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
        # Greedy seed
        greedy = [max(cands, key=lambda x: len(x[0])) for cands in self.valid_overlap_sets]
        pop.append(OverlapChromosome(greedy))
        # Random fill
        while len(pop) < self.ga.population_size:
            ovrs = [random.choice(cands) for cands in self.valid_overlap_sets]
            if self.overlaps_disjoint(ovrs):
                pop.append(OverlapChromosome(ovrs))
        return pop

    def _build_full_oligos(self, chromo: OverlapChromosome) -> List[str]:
        """
        Expand bodies with overlaps to obtain FULL oligo sequences (body + overlaps on both sides).
        """
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

    # ----------------------------- Fitness ------------------------------------

    def _evaluate_fitness(self, chromo: OverlapChromosome) -> float:
        full_oligos = self._build_full_oligos(chromo)
        overlap_seqs = [ov[0] for ov in chromo.overlaps]

        # Very-hard penalties: child full-length range violations
        penalty = 0.0
        if self.child_len_min is not None and self.child_len_max is not None:
            lo, hi = self.child_len_min, self.child_len_max
            for L in (len(s) for s in full_oligos):
                if L < lo:
                    penalty -= VERY_HARD_PENALTY + (lo - L) * HARD_SCALE_PENALTY
                elif L > hi:
                    penalty -= VERY_HARD_PENALTY + (L - hi) * HARD_SCALE_PENALTY

        # Orthogonality of overlaps
        avg_nd, min_nd = pairwise_norm_edit_stats(overlap_seqs, options=self.batch_opts)
        score = ORTHO_AVG_WEIGHT * avg_nd + ORTHO_WORST_WEIGHT * min_nd

        # 3' mismatch & rc k-mer hit proxies
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
        # Nudge toward shorter overlaps (helps increase body sizes and reduce oligo count)
        score -= OVERLAP_LENGTH_PENALTY * sum(len(s) for s in overlap_seqs)

        # Rewards for allowed motifs (if configured)
        if self.allowed_motifs:
            motif_reward = 0.0
            for ov in overlap_seqs:
                for m, w in self.allowed_motifs.items():
                    if not m:
                        continue
                    motif_reward += w * count_motif_occurrences(ov, m)
            score += self.ga.allowed_motif_scale * motif_reward

        # --- NEW (v1.8.0): gentle Tm shaping over overlaps ---
        if overlap_seqs:
            tms = [compute_tm(s, method=self.tm_method, **self.tm_params) for s in overlap_seqs]
            avg_tm = mean(tms)
            std_tm = pstdev(tms) if len(tms) > 1 else 0.0
            # Slight reward for higher average Tm
            score += TM_AVG_BONUS * avg_tm
            # Very slight reward for tighter Tm spread (subtract std)
            score -= TM_STD_PENALTY * std_tm

        score += penalty  # apply penalties last

        chromo.fitness = score

        # Track best overall (clone overlaps so later mutation/crossover won't mutate it)
        if (self._best_overall is None) or (score > (self._best_overall.fitness or float("-inf"))):
            self._best_overall = OverlapChromosome(chromo.overlaps[:])
            self._best_overall.fitness = score

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
        Run the GA and return the best chromosome (GLOBAL best across the entire run)
        and a per-generation best-fitness log. Also stores a `RunSummary` in
        `self.last_run_summary`. Detailed per-generation stats are in `self.progress_detail`.
        """
        population = self.initial_population()
        progress_log: List[float] = []
        self.progress_detail = []
        prev_best = float("-inf")
        gens = max(1, int(self.ga.num_generations))
        self.mutation_rate = self.ga.mutation_rate_initial

        boosts = 0
        boost_gens_total = 0
        early_stopped = False
        early_stop_gen: Optional[int] = None

        early_patience = int(max(1, round(self.ga.stagnation_patience * self.ga.early_stop_multiplier)))

        for gen in range(1, gens + 1):
            # --- Heartbeat: start-of-generation ---
            logger.info(
                "Gen%02d: evaluating %d chromosomes across %d junctions (workers=%d)…",
                gen, len(population), self.num_overlaps, self.workers
            )

            # Evaluate fitness (threads) with mid-gen progress milestones
            total = len(population)
            milestones = {max(1, int(total * p)) for p in (0.10, 0.25, 0.50, 0.75, 1.00)}
            completed = 0

            with ThreadPoolExecutor(max_workers=self.workers) as exe:
                futures = {exe.submit(self._evaluate_fitness, c): c for c in population}
                for fut in as_completed(futures):
                    futures[fut].fitness = fut.result()
                    completed += 1
                    if completed in milestones:
                        logger.info("  Gen%02d progress: %d/%d evaluated (%.0f%%)", gen, completed, total, 100.0 * completed / total)

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
            # per-gen detail; include cum_best (monotonic best-so-far)
            cum_prev = self.progress_detail[-1].get("cum_best", float("-inf")) if self.progress_detail else float("-inf")
            cum_best = max(cum_prev, float(best_f))
            self.progress_detail.append({"gen": float(gen), "best": float(best_f), "mean": float(mean_f), "std": float(std_f), "cum_best": float(cum_best)})

            logger.info(
                "Gen%02d: size=%d uniq=%d (%.1f%%) best=%.3f mean=%.3f median=%.3f worst=%.3f std=%.3f mut=%.2f Δbest=%.3f",
                gen, pop_size, unique_sets, uniq_pct, best_f, mean_f, median_f, worst_f, std_f,
                self.mutation_rate, impr,
            )

            # Stagnation detection & handling
            if best_f > prev_best + self.ga.stagnation_improve_epsilon:
                self._no_improve = 0
                self._no_improve_total = 0
                prev_best = best_f
                self.mutation_rate = max(self.ga.mutation_rate_initial, self.mutation_rate * 0.95)
            else:
                self._no_improve += 1
                self._no_improve_total += 1
                diversity_collapsed = (uniq_pct < self.ga.diversity_min_unique_pct)
                if self._no_improve >= self.ga.stagnation_patience or diversity_collapsed:
                    self.mutation_rate = min(self.ga.mutation_rate_cap, self.mutation_rate * self.ga.stagnation_mutation_mult)
                    self._boost_until_gen = max(self._boost_until_gen, gen + self.ga.stagnation_boost_gens)
                    boosts += 1
                    self._boost_events.append(gen)
                    logger.info(
                        "  stagnation detected (no_improve=%d, uniq=%.1f%%) → boost: mut=%.2f; boost_until_gen=%d",
                        self._no_improve, uniq_pct, self.mutation_rate, self._boost_until_gen,
                    )
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
                        logger.info("  partial restart: kept=%d, rebuilt_to=%d, restarts=%d", keep, len(population), self._restart_count)
                        self._no_improve = 0  # reset recent patience

            if gen <= self._boost_until_gen:
                boost_gens_total += 1

            # Early-stop
            early_patience = int(max(1, round(self.ga.stagnation_patience * self.ga.early_stop_multiplier)))
            if self.ga.enable_early_stop and self._no_improve_total >= early_patience:
                early_stopped = True
                early_stop_gen = gen
                logger.info(
                    "Early-stop triggered at Gen%02d (no improvement for %d consecutive gens ≥ patience %d).",
                    gen, self._no_improve_total, early_patience
                )
                population = sorted(population, key=lambda c: c.fitness or -1e9, reverse=True)
                break

            # Optional hill-climb on the current best
            if self.ga.enable_hill_climb and gen % self.ga.hill_climb_period == 0:
                pre = best.fitness or float("-inf")
                junctions = random.sample(range(self.num_overlaps), min(self.ga.hill_climb_junc, self.num_overlaps))
                for j in junctions:
                    cands = random.sample(self.valid_overlap_sets[j], min(self.ga.hill_climb_cands, len(self.valid_overlap_sets[j])))
                    for cand in cands:
                        trial = OverlapChromosome(best.overlaps[:])
                        trial.overlaps[j] = cand
                        if self.overlaps_disjoint(trial.overlaps):
                            fval = self._evaluate_fitness(trial)  # updates global best if improved
                            if fval > (best.fitness or float("-inf")):
                                best = trial
                post = best.fitness or float("-inf")
                if post > pre:
                    logger.info("  hill-climb improved best: +%.3f → %.3f", post - pre, post)
                    prev_best = max(prev_best, post)
                    self._no_improve = 0
                    self._no_improve_total = 0
                    # Ensure the improved hill-climb result survives into the next generation
                    worst_idx = min(
                        range(len(population)),
                        key=lambda i: population[i].fitness if population[i].fitness is not None else float("inf")
                    )
                    population[worst_idx] = best

            # Build next generation
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
        final_best = self._best_overall or (max(population, key=lambda c: c.fitness or -1e9) if population else None)
        if final_best is None:
            raise RuntimeError("GAOverlapSelector: no final best chromosome available")

        last_unique_sets = len({tuple(ov[0] for ov in population[i].overlaps) for i in range(len(population))}) if population else 0
        last_uniq_pct = 100.0 * last_unique_sets / max(1, len(population))
        gens_run = len(progress_log)

        self.last_run_summary = RunSummary(
            generations_run=gens_run,
            early_stopped=early_stopped,
            early_stop_gen=early_stop_gen if early_stopped else None,
            best_fitness=final_best.fitness or float("-inf"),
            boosts=len(self._boost_events),
            restarts=self._restart_count,
            boost_gens_total=boost_gens_total,
            unique_pct_last=last_uniq_pct,
            boost_gens=self._boost_events[:],
            restart_gens=self._restart_events[:],
        )

        # LOGGING FIX: use %s for lists
        logger.info(
            "GA run summary: gens=%d, early_stop=%s%s, best=%.3f, boosts_total=%d (boost_gens=%d @%s), "
            "restarts=%d @%s, uniq_last=%.1f%%, last_pop_size=%d",
            gens_run,
            str(early_stopped),
            f"@Gen{early_stop_gen:02d}" if early_stopped and early_stop_gen is not None else "",
            self.last_run_summary.best_fitness,
            self.last_run_summary.boost_gens_total,
            len(self.last_run_summary.boost_gens),
            self.last_run_summary.boost_gens,
            self.last_run_summary.restarts,
            self.last_run_summary.restart_gens,
            last_uniq_pct,
            len(population),
        )

        return final_best, progress_log
