# File: backend/app/core/assembly/ga_overlap_selector.py
# Version: v0.9.2

"""
Genetic Algorithm to select orthogonal, specific overlaps between neighboring oligo bodies.
Supports per-level overlap size bounds via constructor arguments.
"""

import random
import multiprocessing
from typing import List, Tuple, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed

from .fitness_utils import (
    average_normalized_distance,
    min_normalized_distance,
    count_3prime_mismatch_events,
    count_rc_kmer_hits,
    reverse_complement
)

# === Configuration constants ===

# Overlap length constraints (bp)
OVERLAP_SIZE_MIN      = 20
OVERLAP_SIZE_MAX      = 35

# GA population parameters
POPULATION_SIZE       = 200
NUM_GENERATIONS       = 3
MUTATION_RATE_INITIAL = 0.35
CROSSOVER_RATE        = 0.6
ELITISM_COUNT         = 1
RANDOM_INJECTION_RATE = 0.05
TOURNAMENT_SIZE       = 3

# Fitness function weights
ORTHO_AVG_WEIGHT      = 100   # weight for average orthogonality
ORTHO_WORST_WEIGHT    = 200   # weight for worst-case orthogonality
MIS3_PENALTY          =   2   # per 3′ mis-annealing event
MISK_PENALTY          =   1   # per k-mer mis-annealing hit
OVERLAP_REWARD        =  10   # reward per overlap in solution

# Hill-climbing local search parameters
ENABLE_HILL_CLIMB     = True
HILL_CLIMB_PERIOD     = 5     # run local search every 5 generations
HILL_CLIMB_JUNC       = 2     # number of junctions to sample
HILL_CLIMB_CANDS      = 10    # candidates per sampled junction


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

    Attributes:
        full_sequence:   The original full DNA string.
        oligo_positions: List of (start,end) positions for each oligo body.
        oligo_seqs:      List of the body sequences themselves.
        num_overlaps:    Number of junctions = len(oligo_positions) - 1.
        valid_overlap_sets: Precomputed list of all possible overlap candidates at each junction.
        cpu_count:       Number of parallel workers to use.
        mutation_rate:   Current mutation rate (adaptive over time).
    """

    def __init__(self,
                 full_sequence:     str,
                 oligo_positions:   List[Tuple[int, int]],
                 oligo_seqs:        List[str],
                 overlap_min_size:  int = OVERLAP_SIZE_MIN,
                 overlap_max_size:  int = OVERLAP_SIZE_MAX):
        self.full_sequence      = full_sequence
        self.oligo_positions    = oligo_positions
        self.oligo_seqs         = oligo_seqs
        self.num_overlaps       = len(oligo_positions) - 1
        self.overlap_min_size   = overlap_min_size
        self.overlap_max_size   = overlap_max_size
        self.valid_overlap_sets = self.extract_valid_overlap_candidates()
        self.cpu_count          = max(1, int(multiprocessing.cpu_count() * 0.9))
        self.mutation_rate      = MUTATION_RATE_INITIAL

    def extract_valid_overlap_candidates(self) -> List[List[Tuple[str, int, int]]]:
        """
        For each junction between oligo bodies, collect all substrings of the full_sequence
        whose lengths range from OVERLAP_SIZE_MIN to OVERLAP_SIZE_MAX, along with their
        absolute start/end coordinates.
        """
        N = len(self.full_sequence)
        all_cands: List[List[Tuple[str,int,int]]] = []

        for i in range(self.num_overlaps):
            # define a search window around the junction
            _, end_i   = self.oligo_positions[i]
            start_j, _ = self.oligo_positions[i+1]
            ws = max(0, end_i - self.overlap_max_size)
            we = min(N, start_j + self.overlap_max_size)
            window = self.full_sequence[ws:we]

            cands: List[Tuple[str, int, int]] = []
            for L in range(self.overlap_min_size, self.overlap_max_size + 1):
                for k in range(len(window) - L + 1):
                    abs_s = ws + k
                    abs_e = abs_s + L
                    seq   = self.full_sequence[abs_s:abs_e]
                    cands.append((seq, abs_s, abs_e))
            all_cands.append(cands)

        return all_cands

    def overlaps_disjoint(self, overlaps: List[Tuple[str, int, int]]) -> bool:
        """
        Check that no two adjacent overlaps in genomic coordinates overlap:
          overlap[i].end <= overlap[i+1].start
        Prevents physically impossible configurations.
        """
        for i in range(len(overlaps) - 1):
            _, _, e1 = overlaps[i]
            _, s2, _ = overlaps[i+1]
            if e1 > s2:
                return False
        return True

    def initial_population(self) -> List[OverlapChromosome]:
        """
        Create an initial population with:
          - one greedy chromosome using the longest available overlaps at each junction
          - the rest sampled randomly, enforcing disjointness of overlap coordinates.
        """
        pop: List[OverlapChromosome] = []

        # Greedy seed: pick the longest overlap at each junction
        greedy = [max(cands, key=lambda x: len(x[0])) for cands in self.valid_overlap_sets]
        pop.append(OverlapChromosome(greedy))

        # Fill the rest randomly (ensuring no coordinate conflicts)
        while len(pop) < POPULATION_SIZE:
            ovrs = [random.choice(cands) for cands in self.valid_overlap_sets]
            if self.overlaps_disjoint(ovrs):
                pop.append(OverlapChromosome(ovrs))

        return pop

    def _build_full_oligos(self, chromo: OverlapChromosome) -> List[str]:
        """
        Given a chromosome of overlaps, reconstruct the full-length oligo sequences by:
          - extending each body with its upstream/downstream overlaps
          - alternating strand orientation: even indices → sense; odd → antisense.
        Returns the list of resulting oligo sequences.
        """
        bodies    = self.oligo_seqs
        positions = self.oligo_positions
        ovs       = chromo.overlaps

        full_oligos: List[str] = []
        for idx, ((s0, e0), _) in enumerate(zip(positions, bodies)):
            prev = ovs[idx-1] if idx > 0 else None
            nxt  = ovs[idx]   if idx < len(bodies)-1 else None

            # genomic start/end incorporate overlaps
            start = prev[1] if prev else s0
            end   = nxt[2]  if nxt else e0

            frag = self.full_sequence[start:end]
            # alternate sense/antisense
            seq = frag if idx % 2 == 0 else reverse_complement(frag)
            full_oligos.append(seq)

        return full_oligos

    def _evaluate_fitness(self, chromo: OverlapChromosome) -> float:
        """
        Fitness = (orthogonality bonus) - (misannealing penalties) + (overlap count reward)

        Orthogonality measures:
          - average normalized edit distance across all overlap pairs
          - minimal normalized edit distance (worst pair)

        Mis-annealing counts:
          - count_3prime_mismatch_events: number of 3′-end binding events
          - count_rc_kmer_hits: number of shared k-mers in reverse orientation

        Returns the computed fitness score.
        """
        # 1) reconstruct full oligos for misannealing checks
        full_oligos  = self._build_full_oligos(chromo)
        overlap_seqs = [ov[0] for ov in chromo.overlaps]

        # Compute orthogonality metrics
        avg_dist   = average_normalized_distance(overlap_seqs)
        worst_dist = min_normalized_distance(overlap_seqs)
        score  = ORTHO_AVG_WEIGHT   * avg_dist
        score += ORTHO_WORST_WEIGHT * worst_dist

        # Mis-annealing penalties (both orientations)
        total_events3 = 0
        total_hitsk   = 0
        for ov in overlap_seqs:
            rc_ov = reverse_complement(ov)
            # count 3' binding events in both orientations
            total_events3 += count_3prime_mismatch_events(ov, full_oligos)
            total_events3 += count_3prime_mismatch_events(rc_ov, full_oligos)
            # count k-mer hits in both orientations
            total_hitsk   += count_rc_kmer_hits(ov, full_oligos)
            total_hitsk   += count_rc_kmer_hits(rc_ov, full_oligos)

        score -= MIS3_PENALTY * total_events3
        score -= MISK_PENALTY * total_hitsk

        # Reward each overlap retained
        score += OVERLAP_REWARD * len(overlap_seqs)

        chromo.fitness = score
        return score

    def tournament(self, population: List[OverlapChromosome]) -> OverlapChromosome:
        """
        Tournament selection: randomly sample TOURNAMENT_SIZE candidates
        and return the one with highest fitness.
        """
        cohort = random.sample(population, min(TOURNAMENT_SIZE, len(population)))
        return max(cohort, key=lambda c: c.fitness or -1e9)

    def select_parents(self, population: List[OverlapChromosome]) \
            -> Tuple[OverlapChromosome, OverlapChromosome]:
        """
        Select two parents independently via tournament selection.
        """
        return self.tournament(population), self.tournament(population)

    def crossover(self, p1: OverlapChromosome, p2: OverlapChromosome) -> OverlapChromosome:
        """
        Single-point crossover between two parents. With probability (1 - CROSSOVER_RATE)
        we simply clone parent1. Otherwise, pick a cut point and combine.
        Then repair any coordinate conflicts by re-sampling conflicting junctions.
        """
        # if there's no real crossover possible (<=1 junction) or we skip by rate, just clone parent1
        if self.num_overlaps < 2 or random.random() > CROSSOVER_RATE:
            return OverlapChromosome(p1.overlaps[:])
        # safe: num_overlaps >= 2, so (1, num_overlaps-1) is non-empty
        pt = random.randint(1, self.num_overlaps - 1)
        child_ov = p1.overlaps[:pt] + p2.overlaps[pt:]

        # Repair disjointness violations
        if not self.overlaps_disjoint(child_ov):
            for i in range(len(child_ov) - 1):
                _, _, e1 = child_ov[i]
                _, s2, _ = child_ov[i+1]
                if e1 > s2:
                    # re-draw those two overlaps at random
                    child_ov[i]   = random.choice(self.valid_overlap_sets[i])
                    child_ov[i+1] = random.choice(self.valid_overlap_sets[i+1])

        # If still invalid, fallback to parent1 clone
        if not self.overlaps_disjoint(child_ov):
            return OverlapChromosome(p1.overlaps[:])

        return OverlapChromosome(child_ov)

    def mutate(self, chromo: OverlapChromosome):
        """
        For each junction, with probability mutation_rate, attempt to
        replace that overlap by a random candidate. Retry up to 10 times
        to preserve coordinate disjointness; otherwise revert.
        """
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
        """
        Main GA loop:
          1) Initialize population.
          2) For each generation:
             - Parallel fitness evaluation
             - Record best fitness & population diversity
             - Optionally perform hill-climbing local search
             - Elitism + tournament + injection to form next generation
          3) Return best solution and fitness log.
        """
        # 1) seed population
        population = self.initial_population()
        progress_log: List[float] = []

        for gen in range(1, NUM_GENERATIONS + 1):
            # 2a) evaluate fitness in parallel
            with ProcessPoolExecutor(max_workers=self.cpu_count) as exe:
                futures = {exe.submit(self._evaluate_fitness, c): c
                           for c in population}
                for fut in as_completed(futures):
                    futures[fut].fitness = fut.result()

            # 2b) diagnostics
            unique_sets = len({tuple(ov[0] for ov in c.overlaps)
                               for c in population})
            best = max(population, key=lambda c: c.fitness or -1e9)
            print(f"Gen{gen}: uniq={unique_sets}, best={best.fitness:.1f}")
            progress_log.append(best.fitness or -1e9)

            # 2c) periodic hill-climbing to refine best solution
            if ENABLE_HILL_CLIMB and gen % HILL_CLIMB_PERIOD == 0:
                junctions = random.sample(range(self.num_overlaps),
                                          min(HILL_CLIMB_JUNC, self.num_overlaps))
                for j in junctions:
                    # sample a few candidates at this junction
                    cands = random.sample(self.valid_overlap_sets[j],
                                          min(HILL_CLIMB_CANDS,
                                              len(self.valid_overlap_sets[j])))
                    for cand in cands:
                        trial = OverlapChromosome(best.overlaps[:])
                        trial.overlaps[j] = cand
                        if self.overlaps_disjoint(trial.overlaps):
                            fval = self._evaluate_fitness(trial)
                            if fval > best.fitness:
                                best = trial

            # 2d) form next generation
            sorted_pop = sorted(population,
                                key=lambda c: c.fitness or -1e9,
                                reverse=True)
            new_pop = sorted_pop[:ELITISM_COUNT]  # carry over elites

            while len(new_pop) < POPULATION_SIZE:
                # random injection
                if random.random() < RANDOM_INJECTION_RATE:
                    ovrs = [random.choice(c)
                            for c in self.valid_overlap_sets]
                    if self.overlaps_disjoint(ovrs):
                        new_pop.append(OverlapChromosome(ovrs))
                else:
                    # crossover + mutation offspring
                    p1, p2 = self.select_parents(population)
                    child  = self.crossover(p1, p2)
                    self.mutate(child)
                    new_pop.append(child)

            population = new_pop

        # 3) return final best
        final_best = max(population, key=lambda c: c.fitness or -1e9)
        return final_best, progress_log
