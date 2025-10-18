# File: backend/app/core/assembly/endaware_planner.py
# Version: v0.3.2
"""
End-aware body planner with pre-scan overlap length estimation.

v0.3.2
- NEW: honor cfg.enforce_even_children_count. If the minimal feasible n is odd
  and the flag is True, prefer n-1 if feasible; otherwise use n+1.

v0.3.1
- Return (positions, hat_final) where hat_final are per-junction minimum overlap
  lengths that passed local filters in a pre-scan.

Returns:
- positions: list[(s,e)] of **body** (pre-overlap) spans that exactly tile [start:end)
- hat_final: list[int] predicted **shortest workable** overlap lengths per junction
"""

from __future__ import annotations

import logging
import math
from statistics import median
from typing import List, Tuple, Optional

from backend.app.config.config_level import LevelConfig
from .fitness_utils import (
    compute_tm,
    gc_fraction,          # returns %GC [0..100]
    max_same_base_run,
    contains_any_motif,
)

log = logging.getLogger(__name__)


def _conservative_internal_bounds(cfg: LevelConfig) -> Tuple[int, int]:
    lb = max(1, cfg.min_children_size - 2 * cfg.overlap_min_size)
    ub = max(1, cfg.max_children_size - 2 * cfg.overlap_max_size)
    if lb > ub:
        lb = ub = max(1, min(lb, ub))
    return lb, ub


def _equal_spaced_bounds(start: int, end: int, n: int) -> List[Tuple[int, int]]:
    L = end - start
    return [(start + math.floor(i * L / n), start + math.floor((i + 1) * L / n)) for i in range(n)]


def _estimate_overlap_lengths_at_junctions(
    full_seq: str,
    positions: List[Tuple[int, int]],
    cfg: LevelConfig,
    tm_method: str,
    tm_params: dict,
) -> List[int]:
    N = len(full_seq)
    hat_o: List[int] = []

    for j in range(len(positions) - 1):
        left_s, left_e = positions[j]
        right_s, right_e = positions[j + 1]
        found: Optional[int] = None

        for L in range(cfg.overlap_min_size, cfg.overlap_max_size + 1):
            ok = False
            # Left-flanking
            a_s = max(0, left_e - L)
            a_e = left_e
            if a_e - a_s == L:
                seqA = full_seq[a_s:a_e]
                ok = _passes_filters(seqA, cfg, tm_method, tm_params)

            # Right-flanking
            if not ok:
                b_s = right_s
                b_e = min(N, right_s + L)
                if b_e - b_s == L:
                    seqB = full_seq[b_s:b_e]
                    ok = _passes_filters(seqB, cfg, tm_method, tm_params)

            if ok:
                found = L
                break

        if found is None:
            found = cfg.overlap_max_size

        hat_o.append(found)

    return hat_o


def _passes_filters(seq: str, cfg: LevelConfig, tm_method: str, tm_params: dict) -> bool:
    if cfg.overlap_tm_min is not None or cfg.overlap_tm_max is not None:
        t = compute_tm(seq, method=tm_method, **tm_params)
        if cfg.overlap_tm_min is not None and t < cfg.overlap_tm_min:
            return False
        if cfg.overlap_tm_max is not None and t > cfg.overlap_tm_max:
            return False

    if cfg.overlap_gc_min is not None or cfg.overlap_gc_max is not None:
        gcp = gc_fraction(seq)
        if cfg.overlap_gc_min is not None and gcp < cfg.overlap_gc_min:
            return False
        if cfg.overlap_gc_max is not None and gcp > cfg.overlap_gc_max:
            return False

    if cfg.overlap_run_min is not None or cfg.overlap_run_max is not None:
        r = max_same_base_run(seq)
        if cfg.overlap_run_min is not None and r < cfg.overlap_run_min:
            return False
        if cfg.overlap_run_max is not None and r > cfg.overlap_run_max:
            return False

    if cfg.overlap_disallowed_motifs and contains_any_motif(seq, cfg.overlap_disallowed_motifs):
        return False

    return True


def _bounds_end_with_o(cfg: LevelConfig, o_len: int) -> Tuple[int, int]:
    lb = max(1, cfg.min_children_size - o_len)
    ub = max(1, cfg.max_children_size - o_len)
    if lb > ub:
        lb = ub = max(1, min(lb, ub))
    return lb, ub


def _bounds_internal_with_o(cfg: LevelConfig, o_left: int, o_right: int) -> Tuple[int, int]:
    lb = max(1, cfg.min_children_size - (o_left + o_right))
    ub = max(1, cfg.max_children_size - (o_left + o_right))
    if lb > ub:
        lb = ub = max(1, min(lb, ub))
    return lb, ub


def _feasible_sum_range_endaware(LBe: int, UBe: int, LBis: List[int], UBis: List[int]) -> Tuple[int, int]:
    if not LBis:
        return (LBe, UBe)
    sum_min = LBe + LBe + sum(LBis)
    sum_max = UBe + UBe + sum(UBis)
    return (sum_min, sum_max)


def _distribute_with_bounds(
    L: int,
    LBe: int,
    UBe: int,
    LBis: List[int],
    UBis: List[int],
) -> List[int]:
    if not LBis:  # n == 1
        x = max(LBe, min(UBe, L))
        return [x]

    n = len(LBis) + 2
    bodies = [LBe] + LBis[:] + [LBe]
    caps   = [(UBe - LBe)] + [UBis[i] - LBis[i] for i in range(len(LBis))] + [(UBe - LBe)]

    total = sum(bodies)
    residual = L - total

    for idx in (0, n - 1):
        if residual <= 0:
            break
        add = min(caps[idx], residual)
        bodies[idx] += add
        caps[idx]   -= add
        residual    -= add

    i = 1
    while residual > 0 and len(LBis) > 0:
        add = min(caps[i], residual)
        if add > 0:
            bodies[i] += add
            caps[i]   -= add
            residual  -= add
        i += 1
        if i >= n - 1:
            i = 1
        if all(c == 0 for c in caps):
            break

    cur = sum(bodies)
    if cur != L:
        delta = L - cur
        direction = 1 if delta > 0 else -1
        pool = ([0, n - 1] + list(range(1, n - 1))) if direction > 0 else (list(range(1, n - 1)) + [0, n - 1])
        for idx in pool:
            if direction > 0:
                cap = caps[idx]
                add = min(cap, abs(delta))
                bodies[idx] += add
                delta -= add
            else:
                lb = LBe if idx in (0, n - 1) else LBis[idx - 1]
                take = min(bodies[idx] - lb, abs(delta))
                bodies[idx] -= take
                delta += take
            if delta == 0:
                break

    return bodies


def plan_bodies_endaware_with_prescan(
    full_seq: str,
    start: int,
    end: int,
    cfg: LevelConfig,
    tm_method: str,
    tm_params: dict,
) -> Tuple[List[Tuple[int, int]], List[int]]:
    """
    Return positions and per-junction **minimum** overlaps (hat_final) that passed
    local filters in a pre-scan. GA should be constrained to L_j >= hat_final[j].
    Also honors cfg.enforce_even_children_count by adjusting n to an even value
    (prefer n-1 if feasible; otherwise n+1).
    """
    L = end - start

    # (1) conservative initial split
    LBi_cons, UBi_cons = _conservative_internal_bounds(cfg)
    n0 = max(1, math.ceil(L / max(1, UBi_cons)))
    positions0 = _equal_spaced_bounds(start, end, n0)

    # (2) pre-scan predicted overlap lengths
    hat_o = _estimate_overlap_lengths_at_junctions(full_seq, positions0, cfg, tm_method, tm_params)
    o_eff = median(hat_o) if hat_o else cfg.overlap_max_size
    o_first = hat_o[0] if hat_o else o_eff
    o_last  = hat_o[-1] if hat_o else o_eff

    # (3) feasible n with end-aware bounds
    LBe, UBe = _bounds_end_with_o(cfg, int(o_first))
    LBi, UBi = _bounds_internal_with_o(cfg, int(o_eff), int(o_eff))

    def feasible_range_for_n(n: int) -> Tuple[int, int]:
        if n == 1:
            return (LBe, UBe)
        LBis = [LBi] * (n - 2)
        UBis = [UBi] * (n - 2)
        return _feasible_sum_range_endaware(LBe, UBe, LBis, UBis)

    # minimal n
    n = 1
    while True:
        smin, smax = feasible_range_for_n(n)
        if smin <= L <= smax:
            break
        n += 1
        if n > max(1, n0 + 10):
            n = n0
            break

    # Honor "enforce_even_children_count" by nudging n to the nearest feasible even
    chosen_n = n
    if cfg.enforce_even_children_count and (n % 2 == 1):
        # Prefer fewer oligos if possible
        if n > 1:
            smin_m1, smax_m1 = feasible_range_for_n(n - 1)
            if smin_m1 <= L <= smax_m1:
                chosen_n = n - 1
        if chosen_n == n:  # try n+1 if n-1 not feasible
            smin_p1, smax_p1 = feasible_range_for_n(n + 1)
            if smin_p1 <= L <= smax_p1:
                chosen_n = n + 1

    # (4) distribute exact sizes → positions
    n_use = chosen_n
    LBis = [LBi] * max(0, n_use - 2)
    UBis = [UBi] * max(0, n_use - 2)
    sizes = _distribute_with_bounds(L, LBe, UBe, LBis, UBis)
    pos = []
    cur = start
    for sz in sizes:
        pos.append((cur, cur + sz))
        cur += sz

    # recompute per-junction for final cuts (the **minimums** GA must honor)
    hat_final = _estimate_overlap_lengths_at_junctions(full_seq, pos, cfg, tm_method, tm_params)

    log.info(
        "End-aware planner: L=%d, n0=%d -> n_min=%d, n_applied=%d%s, o_eff≈%d (first=%d,last=%d), bounds_end=[%d..%d], bounds_int=[%d..%d]",
        L, n0, n, n_use,
        " (even-enforced)" if (cfg.enforce_even_children_count and (n_use % 2 == 0)) else "",
        int(o_eff), int(o_first), int(o_last), LBe, UBe, LBi, UBi
    )

    return pos, hat_final
