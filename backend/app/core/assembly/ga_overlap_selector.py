# File: backend/app/core/assembly/ga_overlap_selector.py
# Version: v0.8.2

import random
import multiprocessing
from typing import List, Tuple, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed

from .fitness_utils import is_orthogonal, has_bad_3prime_binding_fast, rc_kmer_misannealing, reverse_complement

# === Configuration ===
OVERLAP_SIZE_MIN  = 15
OVERLAP_SIZE_MAX  = 25
POPULATION_SIZE   = 100
NUM_GENERATIONS   = 10
MUTATION_RATE     = 0.2
CROSSOVER_RATE    = 0.6
# =====================


class OverlapChromosome:
    """
    overlaps: List of tuples (seq, abs_start, abs_end)
    """
    def __init__(self, overlaps: List[Tuple[str,int,int]]):
        self.overlaps = overlaps
        self.fitness: Optional[float] = None


class GAOverlapSelector:
    """
    Genetic Algorithm to select exact-overlap segments between oligo bodies,
    then evaluates full-length oligos (overlap + body + overlap) for fitness.
    """

    def __init__(self,
                 full_sequence:    str,
                 oligo_positions:  List[Tuple[int,int]],
                 oligo_seqs:       List[str]):
        self.full_sequence   = full_sequence
        self.oligo_positions = oligo_positions
        self.oligo_seqs      = oligo_seqs
        self.num_overlaps    = len(oligo_positions) - 1
        self.valid_overlap_sets = self.extract_valid_overlap_candidates()
        self.cpu_count = max(1, int(multiprocessing.cpu_count() * 0.9))

    def extract_valid_overlap_candidates(self) -> List[List[Tuple[str,int,int]]]:
        """
        For each junction, collect all substrings of lengths MIN..MAX,
        along with their absolute start/end in full_sequence.
        """
        candidates: List[List[Tuple[str,int,int]]] = []
        N = len(self.full_sequence)

        for i in range(self.num_overlaps):
            _, end_i = self.oligo_positions[i]
            start_j, _ = self.oligo_positions[i+1]
            ws = max(0, end_i - OVERLAP_SIZE_MAX)
            we = min(N, start_j + OVERLAP_SIZE_MAX)
            window = self.full_sequence[ws:we]

            ovl_cands: List[Tuple[str,int,int]] = []
            for L in range(OVERLAP_SIZE_MIN, OVERLAP_SIZE_MAX + 1):
                for k in range(0, len(window) - L + 1):
                    abs_s = ws + k
                    abs_e = abs_s + L
                    seq = self.full_sequence[abs_s:abs_e]
                    ovl_cands.append((seq, abs_s, abs_e))
            candidates.append(ovl_cands)

        return candidates

    def initial_population(self) -> List[OverlapChromosome]:
        return [
            OverlapChromosome([random.choice(cands) for cands in self.valid_overlap_sets])
            for _ in range(POPULATION_SIZE)
        ]

    def _build_full_oligos(self, chromo: OverlapChromosome) -> List[str]:
        """
        Build each full-length oligo (including upstream/downstream overlaps),
        assigning alternating strand orientation (+/-).
        """
        bodies    = self.oligo_seqs
        positions = self.oligo_positions
        ovs       = chromo.overlaps  # List of (seq, abs_start, abs_end)

        full_oligos: List[str] = []
        for idx, ((s0, e0), body) in enumerate(zip(positions, bodies)):
            # pick upstream/downstream overlap coords
            prev_ov = ovs[idx-1] if idx > 0 else None
            next_ov = ovs[idx]   if idx < len(bodies)-1 else None

            # determine fragment bounds
            start = prev_ov[1] if prev_ov else s0
            end   = next_ov[2] if next_ov else e0

            # extract full sense sequence
            frag = self.full_sequence[start:end]

            # alternate strand
            if idx % 2 == 0:
                full_oligos.append(frag)
            else:
                # reverse-complement for antisense
                rc = frag.translate(str.maketrans("ACGT", "TGCA"))[::-1]
                full_oligos.append(rc)

        return full_oligos

    def _evaluate_fitness(self, chromo: OverlapChromosome) -> float:
        """
        1) Build full oligos
        2) Score each overlap AND its reverse complement against those oligos
        """
        full_oligos  = self._build_full_oligos(chromo)
        overlap_seqs = [ov[0] for ov in chromo.overlaps]

        score = 0
        # 1. Orthogonality among overlaps
        if not is_orthogonal(overlap_seqs):
            score -= 100

        # 2. Misannealing: test both ov and rc_ov
        for ov in overlap_seqs:
            rc_ov = reverse_complement(ov)
            # if either orientation can misanneal, penalize
            if has_bad_3prime_binding_fast(ov, full_oligos) or \
               has_bad_3prime_binding_fast(rc_ov, full_oligos):
                score -= 15
            if rc_kmer_misannealing(ov, full_oligos) or \
               rc_kmer_misannealing(rc_ov, full_oligos):
                score -= 20

        # Reward per valid overlap
        score += len(overlap_seqs) * 10

        chromo.fitness = score
        return score

    def select_parents(self, population: List[OverlapChromosome]) \
            -> Tuple[OverlapChromosome, OverlapChromosome]:
        sorted_pop = sorted(population,
                            key=lambda c: c.fitness or -1e9,
                            reverse=True)
        return sorted_pop[0], sorted_pop[1]

    def crossover(self, p1: OverlapChromosome, p2: OverlapChromosome) -> OverlapChromosome:
        if random.random() > CROSSOVER_RATE:
            return OverlapChromosome(p1.overlaps[:])
        point = random.randint(1, self.num_overlaps - 1)
        return OverlapChromosome(p1.overlaps[:point] + p2.overlaps[point:])

    def mutate(self, chromo: OverlapChromosome):
        for i in range(self.num_overlaps):
            if random.random() < MUTATION_RATE:
                chromo.overlaps[i] = random.choice(self.valid_overlap_sets[i])

    def evolve(self) -> Tuple[OverlapChromosome, List[float]]:
        population   = self.initial_population()
        progress_log: List[float] = []

        for gen in range(NUM_GENERATIONS):
            # Parallel fitness evaluation
            with ProcessPoolExecutor(max_workers=self.cpu_count) as exe:
                futures = { exe.submit(self._evaluate_fitness, c): c for c in population }
                for fut in as_completed(futures):
                    futures[fut].fitness = fut.result()

            best = max(population, key=lambda c: c.fitness or -1e9)
            progress_log.append(best.fitness or -1e9)

            # Breed next generation
            new_pop: List[OverlapChromosome] = []
            for _ in range(POPULATION_SIZE):
                p1, p2 = self.select_parents(population)
                child  = self.crossover(p1, p2)
                self.mutate(child)
                new_pop.append(child)
            population = new_pop

        best_final = max(population, key=lambda c: c.fitness or -1e9)
        return best_final, progress_log
