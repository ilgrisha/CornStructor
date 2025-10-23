# File: backend/app/core/primer/offtarget.py
# Version: v0.1.0
"""
Off-target attachment counters.

Approach:
- Ungapped sliding comparison of primer vs. target and target RC.
- Count local windows with at most (mismatch_threshold) mismatches for both:
    * full primer length (strict)
    * 3' window of specified length (len3)
- This approximates local alignment with strong gap-open penalties while being fast.
"""

from __future__ import annotations

from typing import Tuple

from .thermodynamics import revcomp


def _count_ungapped_matches(primer: str, subject: str, max_mismatches: int) -> int:
    p = primer.upper()
    s = subject.upper()
    n = len(p)
    if n == 0 or len(s) < 1:
        return 0
    count = 0
    for start in range(0, len(s) - n + 1):
        window = s[start : start + n]
        mism = sum(1 for i in range(n) if p[i] != window[i])
        if mism <= max_mismatches:
            count += 1
    return count


def count_offtargets(primer: str, target: str, max_mismatches_full: int, len3: int, max_mismatches_3p: int) -> Tuple[int, int]:
    """
    Count occurrences where the primer (or 3' window) can attach to target (or RC(target))
    with at most the given mismatch thresholds.

    Returns:
        (full_matches_count, three_prime_matches_count)
    """
    t = target.upper()
    trc = revcomp(t)
    p = primer.upper()
    p3 = p[-len3:] if len3 > 0 else ""

    full = _count_ungapped_matches(p, t, max_mismatches_full) + _count_ungapped_matches(p, trc, max_mismatches_full)
    three_p = 0
    if p3:
        three_p = _count_ungapped_matches(p3, t, max_mismatches_3p) + _count_ungapped_matches(p3, trc, max_mismatches_3p)
    return full, three_p
