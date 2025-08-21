# File: backend/tests/test_batch_fitness.py
# Version: v0.1.0
"""
Tests for batch_fitness parallel helpers.
"""

from __future__ import annotations

import math

from backend.app.core.assembly import batch_fitness as bf
from backend.app.core.assembly import fitness_utils as fu


def test_pairwise_stats_matches_scalar():
    overlaps = ["ACGTAC", "ACGTTC", "TTTTTT", "GGGGGG"]
    avg_s, min_s = bf.pairwise_norm_edit_stats(overlaps, options=bf.BatchOptions(workers=1))
    # scalar reference = compute directly
    ref_avg = fu.average_normalized_distance(overlaps)
    ref_min = fu.min_normalized_distance(overlaps)
    assert math.isclose(avg_s, ref_avg, rel_tol=1e-9, abs_tol=1e-12)
    assert math.isclose(min_s, ref_min, rel_tol=1e-9, abs_tol=1e-12)

    # parallel path (just ensure equal results)
    avg_p, min_p = bf.pairwise_norm_edit_stats(overlaps, options=bf.BatchOptions(workers=4, chunk_size=2))
    assert math.isclose(avg_p, ref_avg, rel_tol=1e-9, abs_tol=1e-12)
    assert math.isclose(min_p, ref_min, rel_tol=1e-9, abs_tol=1e-12)


def test_rc_kmer_hits_total():
    ov = "ACGTACGTACGT"
    oligos = ["TTTTACGTGGGG", "ACGTACGTACGT", "CCCCACGTAAAA"]
    total = bf.rc_kmer_hits_total(ov, oligos, k=4, options=bf.BatchOptions(workers=2, chunk_size=1))
    # Should be >= scalar (equal, since it's just sum of per-oligo)
    total_ref = sum(fu.count_rc_kmer_hits(ov, [x], k=4) for x in oligos)
    assert total == total_ref


def test_tm_batch():
    seqs = ["ACGTACGTACGT", "TTTTTTTTTTTT", "GGGGGGGGGGGG"]
    tms = bf.compute_tm_batch(seqs, options=bf.BatchOptions(workers=2, chunk_size=1))
    assert len(tms) == len(seqs)
    assert all(isinstance(x, float) for x in tms)
