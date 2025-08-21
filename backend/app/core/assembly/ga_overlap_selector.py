# File: backend/app/core/assembly/ga_overlap_selector.py
# Version: v1.4.0

"""
Genetic Algorithm to select orthogonal, specific overlaps between neighboring oligo bodies.

What's new in v1.4.0 (logging & diagnostics):
- Logs resolved performance knobs at initialization:
    INFO: workers, batch_chunk_size, population_size, num_generations, tm_method
- Rich per-generation statistics (INFO):
    * population size
    * number and percentage of unique chromosomes
    * best / mean / median / worst fitness
    * population standard deviation of fitness
    * current mutation rate
    * best-improvement delta vs. previous generation

What's new in v1.3.2 (previous):
- Back-compat: accept any GA-like params object (e.g., GAParamsView from config),
  copying known fields and defaulting missing ones (incl. `workers`, `batch_chunk_size`).
- Keeps v1.3.x perf improvements: threaded evaluation + batch helpers.
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

# Module logger (inherits CLI/root config)
logger = logging.getLogger(__name__)

# === Default GA weights/knobs (can be overridden via GAParameters) ===

# Fitness function weights
ORTHO_AVG_WEIGHT = 100.0   # weight for average orthogonality
ORTHO_WORST_WEIGHT = 200.0 # weight for worst-case orthogonality
MIS3_PENALTY = 2.0         # per 3′ mis-annealing event
MISK_PENALTY = 1.0         # per k-mer mis-annealing hit
OVERLAP_REWARD = 10.0      # reward per overlap in solution
ALLOWED_MOTIF_SCALE = 1.0  # scales contribution of allowed motifs (sum w*count)

# Population parameters
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
    """Container for GA hyperparameters, fitness scales, and perf options."""
    population_size: int = POPULATION_SIZE
    num_generations: int = NUM_GENERATIONS
    mutation_rate_initial: float = MUTATION_RATE_INITIAL
    crossover_rate: float = CROSOVER_RATE if False else CROSSOVER_RATE  # static analyzers safe
    elitism_count: int = ELITISM_COUNT
    random_injection_rate: float = RANDOM_INJECTION_RATE
    tournament_size: int = TOURNAMENT_SIZE
    enable_hill_climb: bool = ENABLE_HILL_CLIMB
    hill_climb_period: int = HILL_CLIMB_PERIOD
    hill_climb_junc: int = HILL_CLIMB_JUNC
    hill_climb_cands: int = HILL_CLIMB_CANDS
    allowed_motif_scale: float = ALLOWED_MOTIF_SCALE
    # Performance tuning for threaded evaluation/batching
    workers: int = 0              # 0 → auto (cpu_count*0.9)
    batch_chunk_size: int = 128   # chunk size for batch helpers


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
                # Keep default if type conversion fails
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
    Represents one candidate solution: a list of overlap segments for each junction.
    Each overlap is a tuple: (sequence, abs_start, abs_end).
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

    Args:
        full_sequence:   The original full DNA string.
        oligo_positions: List[(start,end)] for each oligo body at this level.
        oligo_seqs:      List of the body sequences themselves.
        overlap_min_size, overlap_max_size: allowed candidate lengths (bp).
        tm_min, tm_max:  Tm window in °C; None disables bound(s).
        gc_min, gc_max:  GC **percent** window; None disables.
        run_min, run_max: Allowed range for max homopolymer run; None disables.
        disallowed_motifs: exact motifs to forbid.
        allowed_motifs: {motif: weight} to *reward* in fitness (not a filter).
        max_gap_allowed / enforce_span: compatibility placeholders (span not enforced).
        tm_method: "PRIMER3" (default), "NN", or "Wallace".
        tm_params: dict forwarded to `compute_tm` (expects mM/nM per wrapper).
        ga_params: GAParameters or GA-like object (e.g., GAParamsView). Missing fields default.
    """

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
        ga_params: Optional[object] = None,   # <- allow foreign GA param objects
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

        # Tm settings
        self.tm_method = tm_method
        self.tm_params = tm_params or {}

        # GA knobs (coerced to our dataclass with defaults)
        self.ga: GAParameters = _coerce_ga_params(ga_params)
        self.mutation_rate = self.ga.mutation_rate_initial

        # Worker calc (0.9× CPU) unless GAParameters.workers provided
        auto_workers = max(1, int(multiprocessing.cpu_count() * 0.9))
        self.workers = self.ga.workers or auto_workers
        self.batch_opts = BatchOptions(workers=self.workers, chunk_size=self.ga.batch_chunk_size)

        # Log resolved runtime knobs (visible at CLI INFO level)
        logger.info(
            "GAOverlapSelector init: workers=%d, batch_chunk_size=%d, pop=%d, gens=%d, tm=%s",
            self.workers,
            self.ga.batch_chunk_size,
            self.ga.population_size,
            self.ga.num_generations,
            self.tm_method,
        )

        # Precompute valid candidates per junction (w/ reasons)
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
        For each junction, generate candidates in a bounded window and filter them.
        Raises NoOverlapCandidatesError if a junction ends up empty.
        """
        N = len(self.full_sequence)
        all_sets: List[List[Tuple[str, int, int]]] = []

        for j in range(self.num_overlaps):
            left_s, left_e = self.oligo_positions[j]
            right_s, right_e = self.oligo_positions[j + 1]

            # Candidate search window around the junction
            ws = max(0, left_e - self.overlap_max_size)
            we = min(N, right_s + self.overlap_max_size)
            window = self.full_sequence[ws:we]

            reasons: Dict[str, int] = {
                "length_window": 0,  # placeholder for consistency
                "tm_low": 0,
                "tm_high": 0,
                "gc_low": 0,
                "gc_high": 0,
                "run_low": 0,
                "run_high": 0,
                "motif_disallowed": 0,
            }

            cands: List[Tuple[str, int, int]] = []
            # Enumerate all substrings within length bounds inside the window
            for L in range(self.overlap_min_size, self.overlap_max_size + 1):
                if L <= 0 or L > len(window):
                    continue
                for k in range(0, len(window) - L + 1):
                    abs_s = ws + k
                    abs_e = abs_s + L
                    seq = self.full_sequence[abs_s:abs_e]

                    # ---- Filters (NO SPAN CHECK) ----
                    # Tm
                    if (self.tm_min is not None) or (self.tm_max is not None):
                        t = compute_tm(seq, method=self.tm_method, **self.tm_params)
                        if (self.tm_min is not None) and (t < self.tm_min):
                            reasons["tm_low"] += 1
                            continue
                        if (self.tm_max is not None) and (t > self.tm_max):
                            reasons["tm_high"] += 1
                            continue

                    # GC percent
                    if (self.gc_min is not None) or (self.gc_max is not None):
                        gcp = gc_fraction(seq)  # 0..100
                        if (self.gc_min is not None) and (gcp < self.gc_min):
                            reasons["gc_low"] += 1
                            continue
                        if (self.gc_max is not None) and (gcp > self.gc_max):
                            reasons["gc_high"] += 1
                            continue

                    # Homopolymer run bounds
                    if (self.run_min is not None) or (self.run_max is not None):
                        r = max_same_base_run(seq)
                        if (self.run_min is not None) and (r < self.run_min):
                            reasons["run_low"] += 1
                            continue
                        if (self.run_max is not None) and (r > self.run_max):
                            reasons["run_high"] += 1
                            continue

                    # Disallowed motifs
                    if self.disallowed_motifs and contains_any_motif(seq, self.disallowed_motifs):
                        reasons["motif_disallowed"] += 1
                        continue

                    # Passed all filters
                    cands.append((seq, abs_s, abs_e))

            if not cands:
                raise NoOverlapCandidatesError(
                    junction_index=j,
                    left_end=left_e,
                    right_start=right_s,
                    reasons=reasons,
                )

            all_sets.append(cands)

        return all_sets

    # ---------------- GA primitives ----------------
    # (unchanged below here except for threaded evaluation)

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
            total_hitsk += rc_kmer_hits_total(ov, full_oligos, k=10, options=self.batch_opts)
            total_hitsk += rc_kmer_hits_total(rc_ov, full_oligos, k=10, options=self.batch_opts)
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
                    child_ov[i] = random.choice(self.valid_overlap_sets[i])
                    child_ov[i + 1] = random.choice(self.valid_overlap_sets[i + 1])
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

    def evolve(self) -> Tuple[OverlapChromosome, List[float]]:
        population = self.initial_population()
        progress_log: List[float] = []
        prev_best = float("-inf")
        gens = max(1, int(self.ga.num_generations))
        self.mutation_rate = self.ga.mutation_rate_initial

        for gen in range(1, gens + 1):
            # Evaluate fitness in parallel (threads; C-accelerated hot paths release the GIL)
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
            prev_best = max(prev_best, best_f)

            # Log a concise, information-rich line
            logger.info(
                "Gen%02d: size=%d uniq=%d (%.1f%%) best=%.3f mean=%.3f median=%.3f worst=%.3f std=%.3f mut=%.2f Δbest=%.3f",
                gen, pop_size, unique_sets, uniq_pct, best_f, mean_f, median_f, worst_f, std_f,
                self.mutation_rate, impr,
            )

            # Optional hill-climb
            if self.ga.enable_hill_climb and gen % self.ga.hill_climb_period == 0:
                pre_hc_best = best_f
                junctions = random.sample(range(self.num_overlaps), min(self.ga.hill_climb_junc, self.num_overlaps))
                for j in junctions:
                    cands = random.sample(
                        self.valid_overlap_sets[j],
                        min(self.ga.hill_climb_cands, len(self.valid_overlap_sets[j])),
                    )
                    for cand in cands:
                        trial = OverlapChromosome(best.overlaps[:])
                        trial.overlaps[j] = cand
                        if self.overlaps_disjoint(trial.overlaps):
                            fval = self._evaluate_fitness(trial)
                            if fval > best.fitness:
                                best = trial
                                best_f = fval
                if best_f > pre_hc_best:
                    logger.info("  hill-climb improved best: +%.3f → %.3f", best_f - pre_hc_best, best_f)

            # Form next generation: elitism + offspring + random injection
            sorted_pop = sorted(population, key=lambda c: c.fitness or -1e9, reverse=True)
            new_pop = sorted_pop[: self.ga.elitism_count]
            while len(new_pop) < self.ga.population_size:
                if random.random() < self.ga.random_injection_rate:
                    ovrs = [random.choice(c) for c in self.valid_overlap_sets]
                    if self.overlaps_disjoint(ovrs):
                        new_pop.append(OverlapChromosome(ovrs))
                else:
                    p1, p2 = self.select_parents(population)
                    child = self.crossover(p1, p2)
                    self.mutate(child)
                    new_pop.append(child)
            population = new_pop

        final_best = max(population, key=lambda c: c.fitness or -1e9)
        return final_best, progress_log
