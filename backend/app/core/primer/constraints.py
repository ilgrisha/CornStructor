# File: backend/app/core/primer/constraints.py
# Version: v0.1.0
"""
Filtering constraints for primer candidates.

Includes:
- Homopolymer run limit
- 3' end complementarity checks (self and primer-primer)
- Self-/cross-dimer risk (consecutive & total matches)
"""

from __future__ import annotations

from typing import Tuple

from .thermodynamics import revcomp


def longest_homopolymer(seq: str) -> int:
    best = 0
    run = 0
    prev = ""
    for c in seq.upper():
        if c == prev:
            run += 1
        else:
            run = 1
            prev = c
        if run > best:
            best = run
    return best


def consecutive_complement_runs(a: str, b_rc: str) -> Tuple[int, int]:
    """
    Compute:
     - max consecutive complementary matches when a is aligned to RC of b without gaps
     - total complementary matches in the best ungapped alignment (max over shifts)

    Returns:
        (max_consecutive, max_total)
    """
    a = a.upper()
    b_rc = b_rc.upper()
    max_consec = 0
    max_total = 0

    # Align a against b_rc with all relative shifts (ungapped)
    for shift in range(-len(b_rc) + 1, len(a)):
        consec = 0
        total = 0
        for i in range(len(a)):
            j = i - shift
            if 0 <= j < len(b_rc):
                if a[i] == b_rc[j]:
                    total += 1
                    consec += 1
                    max_consec = max(max_consec, consec)
                else:
                    consec = 0
        max_total = max(max_total, total)

    return max_consec, max_total


def three_prime_match_count(primer: str, target: str, window_len: int) -> int:
    """
    Count number of complementary matches between primer 3' window and target (ungapped),
    maximizing over all ungapped shifts.
    """
    if not primer or window_len <= 0:
        return 0
    w = primer[-window_len:].upper()
    t_rc = revcomp(target.upper())
    # best total matches for w vs RC(target) across shifts (ungapped)
    _, tot = consecutive_complement_runs(w, t_rc)
    return tot


def check_dimer_risk(a: str, b: str) -> Tuple[int, int]:
    """
    Compute dimer risk between a and b (both 5'->3'):
     - returns (max_consecutive, max_total) complementary matches
    """
    return consecutive_complement_runs(a, revcomp(b))
