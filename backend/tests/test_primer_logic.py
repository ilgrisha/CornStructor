# File: backend/tests/test_primer_logic.py
# Version: v0.1.0
"""
Unit tests for the primer subsystem:
- thermodynamics (gc%, revcomp, tm)
- constraints (homopolymer, 3' complementarity, dimer risk)
- offtarget counting
- designer end-to-end on a tiny target window
"""

from __future__ import annotations

import math
import pytest

from app.core.primer.thermodynamics import gc_percent, revcomp, tm_nearest_neighbor
from app.core.primer.constraints import longest_homopolymer, check_dimer_risk, three_prime_match_count
from app.core.primer.offtarget import count_offtargets
from app.core.primer.parameters import PrimerDesignParameters
from app.core.primer.designer import PrimerDesigner


def test_gc_and_revcomp():
    assert math.isclose(gc_percent("ATGC"), 50.0)
    assert revcomp("ATGC") == "GCAT"


def test_homopolymer():
    assert longest_homopolymer("ATTTTGC") == 4
    assert longest_homopolymer("AAAAA") == 5
    assert longest_homopolymer("ATGC") == 1


def test_dimer_risk_basic():
    consec, total = check_dimer_risk("ATGC", "GCAT")  # perfect revcomp
    assert consec >= 4
    assert total >= 4


def test_three_prime_match_count():
    # 3' window "TTT" vs target RC "AAA" should show 3 matches at zero shift
    assert three_prime_match_count("GGGTTT", "AAA", 3) >= 3


def test_offtargets():
    full, three = count_offtargets("ATGC", "TTTATGCGGG", 0, 3, 0)
    # ATGC appears once; 3' window "TGC" appears once
    assert full >= 1
    assert three >= 1


def test_designer_e2e_small_target():
    seq = "TTTAAACCCGGGATGCTAGCTAGCTAGGGTTTAAACCC"
    # target window contains "ATGCTAGC"
    start = seq.index("ATGCTAGC") + 1
    end = start + len("ATGCTAGC") - 1

    params = PrimerDesignParameters(
        primerLengthMin=18, primerLengthMax=18,
        primerTmMin=40.0, primerTmMax=70.0,
        primerGCMin=30.0, primerGCMax=70.0,
        primerHomopolymerMax=4,
        primerThreePrimeEndLength=5,
        primerThreePrimePrimerMatchMax=2,
        primerTargetMatchMax=2,
        primerTargetMatchNumberMax=2,
        primerSecondaryStructureDeltaGMin=-9.0,
        primerTmDifferenceMax=5.0,
    )
    d = PrimerDesigner()
    res = d.design(seq, start, end, params)
    assert res.forward_seq
    assert res.reverse_seq
    assert 30.0 <= res.forward_gc <= 70.0
    assert 40.0 <= res.forward_tm <= 80.0
