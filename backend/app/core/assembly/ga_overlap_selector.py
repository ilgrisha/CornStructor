# File: backend/app/core/assembly/ga_overlap_selector.py
# Version: v1.1.0

"""
Genetic Algorithm to select orthogonal, specific overlaps between neighboring oligo bodies.

New in v1.1.0:
- Per-level 'max_gap_allowed'; enforce overlap to span the junction within an allowed gap
  window; if 'enforce_span' is True then require across-junction coverage with ±max_gap.
- Biopython-based Tm (NN or Wallace) via fitness_utils.compute_tm (delegates to Bio.SeqUtils).
- When a junction loses all candidates after filters, raise NoOverlapCandidatesError and
  store detailed reason tallies (gc_low/gc_high, tm_low/tm_high, run_low/run_high,
  disallowed_motif, span_violation). Caller can adjust division points and retry.
- GA knobs are externally configurable through GAParameters (from GlobalConfig).

If any constraint is None/empty, it is ignored.
"""

import random
import multiprocessing
from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict
from concurrent.futures import ProcessPoolExecutor, as_completed

from .fitness_utils import (
    average_normalized_distance,
    min_normalized_distance,
    count_3prime_mismatch_events,
    count_rc_kmer_hits,
    reverse_complement,
    gc_fraction,
    compute_tm,
    max_same_base_run,
    contains_any_motif,
    count_motif_occurrences,
)


# --------------------
# Exceptions & params
# --------------------

class NoOverlapCandidatesError(RuntimeError):
    """
    Raised when constraint filters eliminate all candidates at a junction.

    Attributes:
        junction_index: index of the junction between child j and j+1.
        left_end:       absolute end coordinate of the left body.
        right_start:    absolute start coordinate of the right body.
        reason_counts:  mapping reason → count during filtering.
    """
    def __init__(self, junction_index: int, left_end: int, right_start: int, reason_counts: Dict[str, int]):
        super().__init__("No valid overlap candidates for junction "
                         f"{junction_index} at {left_end}|{right_start}")
        self.junction_index = junction_index
        self.left_end = left_end
        self.right_start = right_start
        self.reason_counts = reason_counts

    def formatted_reasons(self) -> str:
        return ", ".join(f"{k}:{v}" for k, v in sorted(self.reason_counts.items(), key=lambda x: (-x[1], x[0])))


@dataclass
class GAParameters:
    """Container for GA knobs (tunable via globals)."""
    population_size: int = 200
    num_generations: int = 3
    mutation_rate_initial: float = 0.35
    crossover_rate: float = 0.6
    elitism_count: int = 1
    random_injection_rate: float = 0.05
    tournament_size: int = 3
    enable_hill_climb: bool = True
    hill_climb_period: int = 5
    hill_climb_junc: int = 2
    hill_climb_cands: int = 10
    allowed_motif_scale: float = 1.0


class OverlapChromosome:
    """
    Represents one candidate solution: a list of overlap segments for each junction.
    Each overlap is a tuple: (sequence, abs_start, abs_end).
    """
    def __init__(self, overlaps: List[Tuple[str, int, int]]):
        self.overlaps = overlaps
        self.fitness: Optional[float] = None


class GAOverlapSelector:
    """
    Core GA class for optimizing overlaps.
    """

    def __init__(
        self,
        full_sequence:     str,
        oligo_positions:   List[Tuple[int, int]],
        oligo_seqs:        List[str],
        overlap_min_size:  int,
        overlap_max_size:  int,
        *,
        tm_min: Optional[float] = None,
        tm_max: Optional[float] = None,
        gc_min: Optional[float] = None,
        gc_max: Optional[float] = None,
        run_min: Optional[int] = None,
        run_max: Optional[int] = None,
        disallowed_motifs: Optional[List[str]] = None,
        allowed_motifs: Optional[Dict[str, float]] = None,
        max_gap_allowed: Optional[int] = 0,
        enforce_span: bool = True,
        tm_method: str = "NN",
        tm_params: Optional[Dict[str, float]] = None,
        ga_params: Optional[GAParameters] = None
    ):
        self.full_sequence      = full_sequence
        self.oligo_positions    = oligo_positions
        self.oligo_seqs         = oligo_seqs
        self.num_overlaps       = len(oligo_positions) - 1
        self.overlap_min_size   = overlap_min_size
        self.overlap_max_size   = overlap_max_size

        # Advanced constraints
        self.tm_min             = tm_min
        self.tm_max             = tm_max
        self.gc_min             = gc_min
        self.gc_max             = gc_max
        self.run_min            = run_min
        self.run_max            = run_max
        self.disallowed_motifs  = [m.upper() for m in (disallowed_motifs or [])]
        self.allowed_motifs     = {k.upper(): float(v) for k, v in (allowed_motifs or {}).items()}
        self.max_gap_allowed    = int(max_gap_allowed or 0)
        self.enforce_span       = bool(enforce_span)

        # Tm setup (Biopython-backed)
        self.tm_method          = tm_method
        self.tm_params          = tm_params or {}

        # GA knob setup
        self.params             = ga_params or GAParameters()

        # Parallelism
        self.cpu_count          = max(1, int(multiprocessing.cpu_count() * 0.9))
        self.mutation_rate      = self.params.mutation_rate_initial

        # Candidates
        self.valid_overlap_sets = self.extract_valid_overlap_candidates()

    # ---------- Tuning ----------

    def set_cpu_fraction(self, frac: float):
        """Adjust the CPU worker count fraction (0<frac<=1)."""
        frac = max(0.05, min(1.0, frac))
        self.cpu_count = max(1, int(multiprocessing.cpu_count() * frac))

    # ---------- Candidate extraction with filters ----------

    def extract_valid_overlap_candidates(self) -> List[List[Tuple[str, int, int]]]:
        """
        For each junction, collect all substrings in the search window with lengths
        in [overlap_min_size, overlap_max_size]. Filter by:
          - Disallowed motifs
          - GC fraction range
          - Tm range (Biopython)
          - Homopolymer run (max_same_base_run) range
          - Span across junction within ±max_gap_allowed if enforce_span=True

        If filtered set is empty at a junction, raise NoOverlapCandidatesError with
        reason tallies; callers may adjust division points and retry.
        """
        N = len(self.full_sequence)
        all_cands: List[List[Tuple[str,int,int]]] = []

        for i in range(self.num_overlaps):
            left_s, left_e = self.oligo_positions[i]
            right_s, right_e = self.oligo_positions[i+1]
            # Search window around junction
            ws = max(0, left_e - self.overlap_max_size)
            we = min(N, right_s + self.overlap_max_size)
            window = self.full_sequence[ws:we]

            raw_cands: List[Tuple[str, int, int]] = []
            for L in range(self.overlap_min_size, self.overlap_max_size + 1):
                for k in range(len(window) - L + 1):
                    abs_s = ws + k
                    abs_e = abs_s + L
                    seq   = self.full_sequence[abs_s:abs_e]
                    raw_cands.append((seq, abs_s, abs_e))

            # Filter with reason counting
            reasons: Dict[str, int] = {
                "disallowed_motif": 0,
                "gc_low": 0, "gc_high": 0,
                "tm_low": 0, "tm_high": 0,
                "run_low": 0, "run_high": 0,
                "span_violation": 0,
            }
            filtered: List[Tuple[str, int, int]] = []
            for seq, abs_s, abs_e in raw_cands:
                # Span requirement (overlap should cross the junction; allow small gaps)
                if self.enforce_span:
                    if abs_s > (left_e + self.max_gap_allowed) or abs_e < (right_s - self.max_gap_allowed):
                        reasons["span_violation"] += 1
                        continue

                # Disallowed motifs
                if self.disallowed_motifs and contains_any_motif(seq, self.disallowed_motifs):
                    reasons["disallowed_motif"] += 1
                    continue

                # GC fraction
                if self.gc_min is not None or self.gc_max is not None:
                    gcf = gc_fraction(seq)
                    if self.gc_min is not None and gcf < self.gc_min:
                        reasons["gc_low"] += 1
                        continue
                    if self.gc_max is not None and gcf > self.gc_max:
                        reasons["gc_high"] += 1
                        continue

                # Tm
                if self.tm_min is not None or self.tm_max is not None:
                    tm = compute_tm(seq, method=self.tm_method, **self.tm_params)
                    if self.tm_min is not None and tm < self.tm_min:
                        reasons["tm_low"] += 1
                        continue
                    if self.tm_max is not None and tm > self.tm_max:
                        reasons["tm_high"] += 1
                        continue

                # Same-base run bounds
                if self.run_min is not None or self.run_max is not None:
                    rmax = max_same_base_run(seq)
                    if self.run_min is not None and rmax < self.run_min:
                        reasons["run_low"] += 1
                        continue
                    if self.run_max is not None and rmax > self.run_max:
                        reasons["run_high"] += 1
                        continue

                filtered.append((seq, abs_s, abs_e))

            if not filtered:
                # Report and let caller react
                raise NoOverlapCandidatesError(
                    junction_index=i,
                    left_end=left_e,
                    right_start=right_s,
                    reason_counts={k: v for k, v in reasons.items() if v > 0}
                )

            all_cands.append(filtered)

        return all_cands

    # ---------- GA machinery ----------

    def overlaps_disjoint(self, overlaps: List[Tuple[str, int, int]]) -> bool:
        """Ensure genomic disjointness of adjacent overlaps."""
        for i in range(len(overlaps) - 1):
            _, _, e1 = overlaps[i]
            _, s2, _ = overlaps[i+1]
            if e1 > s2:
                return False
        return True

    def initial_population(self) -> List[OverlapChromosome]:
        """
        Seed with a greedy longest-overlap chromosome, plus random valid chromosomes.
        """
        pop: List[OverlapChromosome] = []

        greedy = [max(cands, key=lambda x: len(x[0])) for cands in self.valid_overlap_sets]
        pop.append(OverlapChromosome(greedy))

        while len(pop) < self.params.population_size:
            ovrs = [random.choice(cands) for cands in self.valid_overlap_sets]
            if self.overlaps_disjoint(ovrs):
                pop.append(OverlapChromosome(ovrs))

        return pop

    def _build_full_oligos(self, chromo: OverlapChromosome) -> List[str]:
        """
        Build full fragment sequences by adding overlaps to bodies and alternating strand.
        """
        bodies    = self.oligo_seqs
        positions = self.oligo_positions
        ovs       = chromo.overlaps

        full_oligos: List[str] = []
        for idx, ((s0, e0), _) in enumerate(zip(positions, bodies)):
            prev = ovs[idx-1] if idx > 0 else None
            nxt  = ovs[idx]   if idx < len(bodies)-1 else None

            start = prev[1] if prev else s0
            end   = nxt[2]  if nxt else e0

            frag = self.full_sequence[start:end]
            seq = frag if idx % 2 == 0 else reverse_complement(frag)
            full_oligos.append(seq)

        return full_oligos

    def _evaluate_fitness(self, chromo: OverlapChromosome) -> float:
        """
        Fitness = (orthogonality bonuses) - (misannealing penalties)
                  + (overlap count reward) + (allowed-motif bonus)
        """
        full_oligos  = self._build_full_oligos(chromo)
        overlap_seqs = [ov[0] for ov in chromo.overlaps]

        # Orthogonality
        avg_dist   = average_normalized_distance(overlap_seqs)
        worst_dist = min_normalized_distance(overlap_seqs)
        score  = 100 * avg_dist
        score += 200 * worst_dist

        # Mis-annealing penalties (both orientations)
        total_events3 = 0
        total_hitsk   = 0
        for ov in overlap_seqs:
            rc_ov = reverse_complement(ov)
            total_events3 += count_3prime_mismatch_events(ov, full_oligos)
            total_events3 += count_3prime_mismatch_events(rc_ov, full_oligos)
            total_hitsk   += count_rc_kmer_hits(ov, full_oligos)
            total_hitsk   += count_rc_kmer_hits(rc_ov, full_oligos)

        score -= 2 * total_events3
        score -= 1 * total_hitsk

        # Reward each overlap retained
        score += 10 * len(overlap_seqs)

        # Allowed-motif bonus
        if self.allowed_motifs:
            motif_bonus = 0.0
            for ov in overlap_seqs:
                s = ov.upper()
                for motif, weight in self.allowed_motifs.items():
                    if motif:
                        cnt = count_motif_occurrences(s, motif)
                        if cnt:
                            motif_bonus += self.params.allowed_motif_scale * weight * cnt
            score += motif_bonus

        chromo.fitness = score
        return score

    def tournament(self, population: List[OverlapChromosome]) -> OverlapChromosome:
        """Tournament selection."""
        k = max(2, min(self.params.tournament_size, len(population)))
        cohort = random.sample(population, k)
        return max(cohort, key=lambda c: c.fitness or -1e9)

    def select_parents(self, population: List[OverlapChromosome]) -> Tuple[OverlapChromosome, OverlapChromosome]:
        """Select two parents independently."""
        return self.tournament(population), self.tournament(population)

    def crossover(self, p1: OverlapChromosome, p2: OverlapChromosome) -> OverlapChromosome:
        """Single-point crossover with repair."""
        if self.num_overlaps < 2 or random.random() > self.params.crossover_rate:
            return OverlapChromosome(p1.overlaps[:])
        pt = random.randint(1, self.num_overlaps - 1)
        child_ov = p1.overlaps[:pt] + p2.overlaps[pt:]

        if not self.overlaps_disjoint(child_ov):
            for i in range(len(child_ov) - 1):
                _, _, e1 = child_ov[i]
                _, s2, _ = child_ov[i+1]
                if e1 > s2:
                    child_ov[i]   = random.choice(self.valid_overlap_sets[i])
                    child_ov[i+1] = random.choice(self.valid_overlap_sets[i+1])

        if not self.overlaps_disjoint(child_ov):
            return OverlapChromosome(p1.overlaps[:])

        return OverlapChromosome(child_ov)

    def mutate(self, chromo: OverlapChromosome):
        """Randomly mutate overlaps while preserving disjointness."""
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
        """Main GA loop with optional hill-climbing."""
        population = self.initial_population()
        progress_log: List[float] = []

        for gen in range(1, self.params.num_generations + 1):
            # Evaluate in parallel
            with ProcessPoolExecutor(max_workers=self.cpu_count) as exe:
                futures = {exe.submit(self._evaluate_fitness, c): c for c in population}
                for fut in as_completed(futures):
                    futures[fut].fitness = fut.result()

            # Diagnostics
            unique_sets = len({tuple(ov[0] for ov in c.overlaps) for c in population})
            best = max(population, key=lambda c: c.fitness or -1e9)
            print(f"Gen{gen}: uniq={unique_sets}, best={best.fitness:.1f}")
            progress_log.append(best.fitness or -1e9)

            # Hill-climbing
            if self.params.enable_hill_climb and gen % self.params.hill_climb_period == 0:
                juncs = random.sample(range(self.num_overlaps), min(self.params.hill_climb_junc, self.num_overlaps))
                for j in juncs:
                    cands = random.sample(self.valid_overlap_sets[j],
                                          min(self.params.hill_climb_cands, len(self.valid_overlap_sets[j])))
                    for cand in cands:
                        trial = OverlapChromosome(best.overlaps[:])
                        trial.overlaps[j] = cand
                        if self.overlaps_disjoint(trial.overlaps):
                            fval = self._evaluate_fitness(trial)
                            if fval > best.fitness:
                                best = trial

            # Next generation
            sorted_pop = sorted(population, key=lambda c: c.fitness or -1e9, reverse=True)
            new_pop = sorted_pop[: self.params.elitism_count]
            while len(new_pop) < self.params.population_size:
                if random.random() < self.params.random_injection_rate:
                    ovrs = [random.choice(c) for c in self.valid_overlap_sets]
                    if self.overlaps_disjoint(ovrs):
                        new_pop.append(OverlapChromosome(ovrs))
                else:
                    p1, p2 = self.select_parents(population)
                    child  = self.crossover(p1, p2)
                    self.mutate(child)
                    new_pop.append(child)

            population = new_pop

        final_best = max(population, key=lambda c: c.fitness or -1e9)
        return final_best, progress_log
