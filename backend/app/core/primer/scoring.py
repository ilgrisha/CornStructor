# File: backend/app/core/primer/scoring.py
# Version: v0.1.0
"""
Composite scoring for primer pairs.

Lower score is better. Components:
- |Tm - targetTm| for each primer + |ΔTm between primers|
- GC deviation from midpoint of allowed range
- Self-/cross-dimer penalties (consecutive and total)
- 3' end complementarity penalties
- Off-target penalties (full + 3' window counts)
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class ScoreContext:
    # Desired/limits for penalties
    tm_min: float
    tm_max: float
    gc_min: float
    gc_max: float
    wTm: float
    wGC: float
    wDimer: float
    w3p: float
    wOff: float
    pair_dt_max: float


def _clamp(a: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, a))


def _gc_mid(gc_min: float, gc_max: float) -> float:
    return (gc_min + gc_max) / 2.0


def score_pair(
    ctx: ScoreContext,
    f_tm: float,
    f_gc: float,
    r_tm: float,
    r_gc: float,
    self_dimer_max_consec: int,
    self_dimer_total: int,
    cross_dimer_max_consec: int,
    cross_dimer_total: int,
    three_p_self: int,
    three_p_cross: int,
    off_full: int,
    off_3p: int,
) -> float:
    # Tm penalties: distance from allowed band + ΔTm
    f_tm_pen = 0.0 if ctx.tm_min <= f_tm <= ctx.tm_max else min(abs(f_tm - ctx.tm_min), abs(f_tm - ctx.tm_max))
    r_tm_pen = 0.0 if ctx.tm_min <= r_tm <= ctx.tm_max else min(abs(r_tm - ctx.tm_min), abs(r_tm - ctx.tm_max))
    dtm = abs(f_tm - r_tm)
    dtm_pen = max(0.0, dtm - ctx.pair_dt_max)

    # GC deviation from midpoint
    gmid = _gc_mid(ctx.gc_min, ctx.gc_max)
    f_gc_pen = abs(f_gc - gmid)
    r_gc_pen = abs(r_gc - gmid)

    # Dimer penalties (weight consecutive higher)
    dimer_pen = (2.0 * (self_dimer_max_consec + cross_dimer_max_consec)) + 0.5 * (self_dimer_total + cross_dimer_total)

    # 3' complementarity penalties
    three_p_pen = 1.5 * three_p_self + 2.0 * three_p_cross

    # Off-target
    off_pen = 1.0 * off_full + 2.0 * off_3p

    total = (
        ctx.wTm * (f_tm_pen + r_tm_pen + dtm_pen)
        + ctx.wGC * (f_gc_pen + r_gc_pen)
        + ctx.wDimer * dimer_pen
        + ctx.w3p * three_p_pen
        + ctx.wOff * off_pen
    )
    return total
