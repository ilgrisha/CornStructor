# File: backend/tests/test_fitness_utils.py
# Version: v0.1.1

"""
Unit tests for C-accelerated fitness utilities.

v0.1.1:
- No code changes to assertions; add a warnings filter for primer3 just in case.
"""

import math
import warnings

import pytest

from backend.app.core.assembly.fitness_utils import (
    _PRIMER3_AVAILABLE,
    average_normalized_distance,
    min_normalized_distance,
    normalized_edit_distance,
    extract_kmer_hashes,
    count_rc_kmer_hits,
    count_3prime_mismatch_events,
    gc_fraction,
    max_same_base_run,
    contains_any_motif,
    count_motif_occurrences,
    compute_tm,
    reverse_complement,
)


def test_edit_distance_normalization_basic():
    assert normalized_edit_distance("AAAA", "AAAA") == 0.0
    assert math.isclose(normalized_edit_distance("AAAA", "AAAT"), 1 / 4)


def test_pairwise_metrics_simple():
    overlaps = ["ACGTAC", "ACGTTC", "TTTTTT"]
    avg = average_normalized_distance(overlaps)
    mmin = min_normalized_distance(overlaps)
    assert 0.0 < mmin <= avg <= 1.0


def test_kmer_hashes_and_hits():
    a = "ACGTACGTAC"
    b = reverse_complement(a)
    k = 4
    ha = extract_kmer_hashes(a, k)
    assert len(ha) > 0
    hits = count_rc_kmer_hits(a, [b], k=k)
    assert hits > 0


def test_3prime_misanneal_events():
    oligos = ["ACGTACGTACGT", "TTTTACGTGGGG", "CCCCCCCC"]
    overlap = "GGGGACGT"
    events = count_3prime_mismatch_events(overlap, oligos, three_prime_len=4, max_mismatches=0)
    assert events >= 1


def test_biophysical_helpers():
    s = "AAACCCGGGTTT"
    assert math.isclose(gc_fraction(s), 6 / 12 * 100.0)
    assert max_same_base_run("AAATTTGG") == 3
    assert contains_any_motif("ACGTACGT", ["TAC", "GGG"]) is True
    assert count_motif_occurrences("AAAA", "AA") == 3


@pytest.mark.parametrize("method", ["PRIMER3", "NN", "Wallace"])
def test_compute_tm_basic(method):
    s = "ACGTACGTACGT"
    if method == "PRIMER3" and not _PRIMER3_AVAILABLE:
        pytest.skip("primer3 not available in test env")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t = compute_tm(s, method=method)
    assert isinstance(t, float)
    assert t > 0
