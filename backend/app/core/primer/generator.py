# File: backend/app/core/primer/generator.py
# Version: v0.1.0
"""
Candidate primer generator.

- Forward candidates: substrings from left to right within [min_len, max_len]
- Reverse candidates: substrings from right to left within [min_len, max_len], then reverse-complemented
"""

from __future__ import annotations

from typing import Iterable, List, Tuple

from .thermodynamics import revcomp


def generate_forward_candidates(target: str, min_len: int, max_len: int) -> List[str]:
    s = target.upper()
    out: List[str] = []
    for L in range(min_len, max_len + 1):
        if L > len(s):
            break
        for i in range(0, len(s) - L + 1):
            out.append(s[i : i + L])
    return out


def generate_reverse_candidates(target: str, min_len: int, max_len: int) -> List[str]:
    s = target.upper()
    out: List[str] = []
    for L in range(min_len, max_len + 1):
        if L > len(s):
            break
        for i in range(0, len(s) - L + 1):
            seg = s[i : i + L]
            out.append(revcomp(seg))
    return out
