# File: backend/app/core/primer/alignments.py
# Version: v0.1.0
"""
ASCII alignment renderers for UI visualization.

Stub strings that keep the API stable; real implementation will use local alignment
and self-/cross-dimer scans in a later step.
"""

from __future__ import annotations


class AlignmentRenderer:
    def forward_vs_target(self, forward: str, target: str) -> str:
        return "|||||   ||||  ||||"  # placeholder

    def reverse_vs_target_rc(self, reverse: str, target_rc: str) -> str:
        return "|||  ||||||  |||||"  # placeholder

    def self_dimer(self, primer: str) -> str:
        return "|||||   ||||"  # placeholder

    def cross_dimer(self, forward: str, reverse_rc: str) -> str:
        return "||| |||  |||"  # placeholder
