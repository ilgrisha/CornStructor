# File: backend/app/config/config_level.py
# Version: v0.4.0
"""
Per-level configuration loader.

v0.4.0
- Moved from backend.app.core.config -> backend.app.config.config_level.
- No functional changes; preserved schema including `max_gap_allowed`.
- Added stricter typing and input normalization for allowed motifs.

This module reads `levels.json` (an array of level entries) and produces a
dictionary of LevelConfig keyed by the integer level number.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Any


@dataclass(frozen=True)
class LevelConfig:
    """
    Configuration for a single assembly level.

    Attributes:
        level: Integer tree level (0 = root).
        fragment_min_size / fragment_max_size: Target size range (bp) for fragments at this level.
        min_children / max_children: Allowed number of children fragments per parent at this level.
        overlap_min_size / overlap_max_size: Candidate overlap length bounds (bp).
        overlap_tm_min / overlap_tm_max: Optional acceptable Tm range for overlaps (°C).
        overlap_gc_min / overlap_gc_max: Optional acceptable GC% range for overlaps.
        overlap_run_min / overlap_run_max: Optional bounds on homopolymer runs (e.g. max run of same base).
        overlap_disallowed_motifs: List of substrings that must NOT appear in overlaps.
        overlap_allowed_motifs: Mapping of preferred motifs -> weight (>0, default 1.0). Empty means no preference.
        max_gap_allowed: Maximum number of bp an overlap may miss the junction by (>=0). 0 means ideally spans.
    """
    level: int

    fragment_min_size: int
    fragment_max_size: int
    min_children: int
    max_children: int

    overlap_min_size: int
    overlap_max_size: int

    overlap_tm_min: Optional[float] = None
    overlap_tm_max: Optional[float] = None
    overlap_gc_min: Optional[float] = None
    overlap_gc_max: Optional[float] = None
    overlap_run_min: Optional[int] = None
    overlap_run_max: Optional[int] = None

    overlap_disallowed_motifs: List[str] = field(default_factory=list)
    overlap_allowed_motifs: Dict[str, float] = field(default_factory=dict)

    max_gap_allowed: int = 0


def _normalize_allowed_motifs(value: Optional[Any]) -> Dict[str, float]:
    """
    Normalize the overlap_allowed_motifs input value, accepting:
      - dict[str, number] -> dict[str, float]
      - list[str] -> dict[str, 1.0] for each element
      - None / missing -> {}
    """
    if value is None:
        return {}
    if isinstance(value, dict):
        out: Dict[str, float] = {}
        for k, v in value.items():
            if not isinstance(k, str):
                continue
            try:
                out[k] = float(v)
            except Exception:
                out[k] = 1.0
        return out
    if isinstance(value, list):
        return {str(k): 1.0 for k in value if isinstance(k, str)}
    try:
        return {str(value): 1.0}
    except Exception:
        return {}


def load_levels_config(path: Path | str) -> Dict[int, LevelConfig]:
    """
    Load an array of per-level configuration entries from JSON.
    Returns a mapping of level -> LevelConfig.

    Args:
        path: Path to levels.json.
    """
    p = Path(path)
    data = json.loads(p.read_text())

    if not isinstance(data, list):
        raise ValueError("levels.json root must be a JSON array")

    levels: Dict[int, LevelConfig] = {}
    for entry in data:
        try:
            level = int(entry["level"])
            fragment_min_size = int(entry["fragment_min_size"])
            fragment_max_size = int(entry["fragment_max_size"])
            min_children = int(entry["min_children"])
            max_children = int(entry["max_children"])
            overlap_min_size = int(entry["overlap_min_size"])
            overlap_max_size = int(entry["overlap_max_size"])
        except KeyError as ke:
            raise KeyError(f"Missing required key in level entry: {ke!s}") from None

        cfg = LevelConfig(
            level=level,
            fragment_min_size=fragment_min_size,
            fragment_max_size=fragment_max_size,
            min_children=min_children,
            max_children=max_children,
            overlap_min_size=overlap_min_size,
            overlap_max_size=overlap_max_size,
            overlap_tm_min=_opt_float(entry.get("overlap_tm_min")),
            overlap_tm_max=_opt_float(entry.get("overlap_tm_max")),
            overlap_gc_min=_opt_float(entry.get("overlap_gc_min")),
            overlap_gc_max=_opt_float(entry.get("overlap_gc_max")),
            overlap_run_min=_opt_int(entry.get("overlap_run_min")),
            overlap_run_max=_opt_int(entry.get("overlap_run_max")),
            overlap_disallowed_motifs=list(entry.get("overlap_disallowed_motifs", []) or []),
            overlap_allowed_motifs=_normalize_allowed_motifs(entry.get("overlap_allowed_motifs")),
            max_gap_allowed=int(entry.get("max_gap_allowed", 0)),
        )
        levels[cfg.level] = cfg

    return levels


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
