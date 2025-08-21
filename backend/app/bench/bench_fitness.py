# File: backend/app/bench/bench_fitness.py
# Version: v0.1.0

"""
Micro-benchmarks for fitness_utils hot paths.

Usage:
    PYTHONPATH=$(pwd) backend/.venv/bin/python backend/app/bench/bench_fitness.py

Or:
    make bench

What it measures:
- Normalized edit distances over N overlaps (pairwise)
- RC k-mer hits across M oligos
- Batch Tm computations (primer3 if available)

Tune N/M/K and sequence lengths to your real workloads.
"""

from __future__ import annotations

import random
import string
import time
from statistics import mean

from backend.app.core.assembly import fitness_utils as fu


def rand_dna(n: int) -> str:
    return "".join(random.choice("ACGT") for _ in range(n))


def timeit(fn, repeat=5):
    times = []
    for _ in range(repeat):
        t0 = time.perf_counter()
        fn()
        times.append(time.perf_counter() - t0)
    return min(times), mean(times)


def bench_edit_distance(n_overlaps=200, overlap_len=30):
    overlaps = [rand_dna(overlap_len) for _ in range(n_overlaps)]
    def run():
        _ = fu.average_normalized_distance(overlaps)
        _ = fu.min_normalized_distance(overlaps)
    return timeit(run)


def bench_kmer_hits(k=10, overlap_len=40, m_oligos=400, oligo_len=60):
    overlap = rand_dna(overlap_len)
    oligos = [rand_dna(oligo_len) for _ in range(m_oligos)]
    def run():
        _ = fu.count_rc_kmer_hits(overlap, oligos, k=k)
    return timeit(run)


def bench_tm(n_primers=500, primer_len=22, method="PRIMER3"):
    primers = [rand_dna(primer_len) for _ in range(n_primers)]
    def run():
        for s in primers:
            _ = fu.compute_tm(s, method=method)
    return timeit(run, repeat=3)


def main():
    random.seed(1337)
    print("== fitness_utils micro-bench ==")
    best, avg = bench_edit_distance()
    print(f"pairwise distances: best {best:.3f}s, avg {avg:.3f}s")
    best, avg = bench_kmer_hits()
    print(f"RC k-mer hits:      best {best:.3f}s, avg {avg:.3f}s")
    try:
        best, avg = bench_tm(method="PRIMER3")
        print(f"Tm (primer3):       best {best:.3f}s, avg {avg:.3f}s")
    except Exception as e:
        print(f"Tm (primer3) skipped: {e}")
        best, avg = bench_tm(method="NN")
        print(f"Tm (NN):            best {best:.3f}s, avg {avg:.3f}s")


if __name__ == "__main__":
    main()
