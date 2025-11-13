# File: backend/app/config/config_level.py
# Version: v0.7.0
"""
Per-level configuration loader.

v0.7.0
- NEW optional boolean per level:
    * enforce_even_children_count (bool, default False): if True, the planner
      will enforce an even number of children (oligos). When the minimal
      feasible count is odd, it prefers n-1 if feasible; otherwise uses n+1.

v0.6.0
- NEW optional field per level:
    * force_body_count_deduction (int ≥ 0, default 0): after computing the minimal
      feasible number of bodies (children) for the fragment, the planner will
      attempt to reduce that count by up to this amount while preserving
      feasibility (size bounds with estimated overlaps). If 0, no deduction.
- Docstring updated.

v0.5.1
- Backward-compatible loader: accepts either a JSON array OR an object with "levels":[...].
- Child-count settings removed in favor of child-size constraints.
- NEW required fields per level:
    * min_children_size / max_children_size  (FULL oligo length, i.e., body + overlaps)
- Optional:
    * body_max_safety_buffer  (keeps pre-overlap bodies under max_children_size with slack)
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional


@dataclass(frozen=True)
class LevelConfig:
    level: int

    # Parent fragment size constraints (validate fragments at this level)
    fragment_min_size: int
    fragment_max_size: int

    # FULL child (oligo) size constraints (post-overlap)
    min_children_size: int
    max_children_size: int

    # Overlap length constraints
    overlap_min_size: int
    overlap_max_size: int

    # Optional overlap quality constraints (kept for compatibility)
    overlap_tm_min: Optional[float] = None
    overlap_tm_max: Optional[float] = None
    overlap_gc_min: Optional[float] = None  # percent [0..100]
    overlap_gc_max: Optional[float] = None  # percent [0..100]
    overlap_run_min: Optional[int] = None
    overlap_run_max: Optional[int] = None
    overlap_disallowed_motifs: Optional[List[str]] = None
    overlap_allowed_motifs: Optional[Dict[str, float]] = None

    # Legacy compat (assembler passes enforce_span=False anyway)
    max_gap_allowed: int = 0

    # Optional slack to keep pre-overlap bodies below max
    body_max_safety_buffer: Optional[int] = None

    # Attempt to force-reduce the number of bodies by up to this amount (0 = no deduction)
    force_body_count_deduction: int = 0

    # NEW: enforce even number of children/oligos (False by default)
    enforce_even_children_count: bool = False


def _opt_float(v: Any) -> Optional[float]:
    if v is None:
        return None
    try:
        return float(v)
    except Exception:
        return None


def _opt_int(v: Any) -> Optional[int]:
    if v is None:
        return None
    try:
        return int(v)
    except Exception:
        return None


def _load_levels_list(data: Any) -> List[dict]:
    """Accept either a list[...] or {"levels":[...]} root structure."""
    if isinstance(data, list):
        return data
    if isinstance(data, dict) and "levels" in data and isinstance(data["levels"], list):
        return data["levels"]
    raise ValueError("levels.json must be a JSON array or an object with a 'levels' array.")


def load_levels_config(path: Path | str) -> Dict[int, LevelConfig]:
    """
    Load levels config into {level: LevelConfig}.
    """
    p = Path(path)
    raw = json.loads(p.read_text(encoding="utf-8"))
    items = _load_levels_list(raw)

    out: Dict[int, LevelConfig] = {}
    for entry in items:
        try:
            level = int(entry["level"])
            fragment_min_size = int(entry["fragment_min_size"])
            fragment_max_size = int(entry["fragment_max_size"])
            overlap_min_size = int(entry["overlap_min_size"])
            overlap_max_size = int(entry["overlap_max_size"])
            # required
            min_children_size = int(entry["min_children_size"])
            max_children_size = int(entry["max_children_size"])
        except KeyError as ke:
            raise KeyError(f"Missing required key in level entry: {ke!s}") from None

        if fragment_min_size <= 0 or fragment_max_size <= 0 or fragment_min_size > fragment_max_size:
            raise ValueError(f"[L{level}] invalid fragment bounds: {fragment_min_size}..{fragment_max_size}")
        if overlap_min_size <= 0 or overlap_max_size <= 0 or overlap_min_size > overlap_max_size:
            raise ValueError(f"[L{level}] invalid overlap bounds: {overlap_min_size}..{overlap_max_size}")
        if min_children_size <= 0 or max_children_size <= 0 or min_children_size > max_children_size:
            raise ValueError(f"[L{level}] invalid child (full oligo) bounds: {min_children_size}..{max_children_size}")

        tm_min = _opt_float(entry.get("overlap_tm_min"))
        tm_max = _opt_float(entry.get("overlap_tm_max"))
        gc_min = _opt_float(entry.get("overlap_gc_min"))
        gc_max = _opt_float(entry.get("overlap_gc_max"))
        run_min = _opt_int(entry.get("overlap_run_min"))
        run_max = _opt_int(entry.get("overlap_run_max"))
        max_gap_allowed = int(entry.get("max_gap_allowed", 0))

        disallowed = entry.get("overlap_disallowed_motifs") or []
        allowed = entry.get("overlap_allowed_motifs") or {}

        body_buf = entry.get("body_max_safety_buffer")
        body_buf = None if body_buf is None else int(body_buf)
        if body_buf is not None and body_buf < 0:
            raise ValueError(f"[L{level}] body_max_safety_buffer must be >= 0")

        force_deduct = int(entry.get("force_body_count_deduction", 0))
        if force_deduct < 0:
            raise ValueError(f"[L{level}] force_body_count_deduction must be >= 0")

        enforce_even = bool(entry.get("enforce_even_children_count", False))

        cfg = LevelConfig(
            level=level,
            fragment_min_size=fragment_min_size,
            fragment_max_size=fragment_max_size,
            min_children_size=min_children_size,
            max_children_size=max_children_size,
            overlap_min_size=overlap_min_size,
            overlap_max_size=overlap_max_size,
            overlap_tm_min=tm_min,
            overlap_tm_max=tm_max,
            overlap_gc_min=gc_min,
            overlap_gc_max=gc_max,
            overlap_run_min=run_min,
            overlap_run_max=run_max,
            overlap_disallowed_motifs=[str(m).upper() for m in disallowed],
            overlap_allowed_motifs={str(k).upper(): float(v) for k, v in allowed.items()},
            max_gap_allowed=max_gap_allowed,
            body_max_safety_buffer=body_buf,
            force_body_count_deduction=force_deduct,
            enforce_even_children_count=enforce_even,
        )
        out[level] = cfg

    return out
