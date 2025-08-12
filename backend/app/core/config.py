# File: backend/app/core/config.py
# Version: v0.3.0

"""
Per-level configuration loader.

New in v0.3.0:
- Add per-level 'max_gap_allowed' (int, default 0) to control how far an overlap
  may miss spanning the junction.
- Advanced overlap constraints (Tm/GC/run/motifs) kept from prior version.
"""

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional


@dataclass
class LevelConfig:
    """
    Configuration for one level of the hierarchical assembly tree.
    """
    level: int
    fragment_min_size:   int
    fragment_max_size:   int
    min_children:        int
    max_children:        int
    overlap_min_size:    int
    overlap_max_size:    int

    # Optional advanced overlap constraints
    overlap_tm_min:      Optional[float] = None
    overlap_tm_max:      Optional[float] = None
    overlap_gc_min:      Optional[float] = None
    overlap_gc_max:      Optional[float] = None
    overlap_run_min:     Optional[int]   = None
    overlap_run_max:     Optional[int]   = None

    # Motif lists
    overlap_disallowed_motifs: List[str]        = field(default_factory=list)
    overlap_allowed_motifs:    Dict[str, float] = field(default_factory=dict)

    # New: allowed gap around the junction for an overlap to be considered spanning
    max_gap_allowed:     int = 0


def _normalize_allowed_motifs(raw_val) -> Dict[str, float]:
    """Accept dict motif->weight OR list of {'motif','weight'}; normalize to dict."""
    if raw_val is None:
        return {}
    if isinstance(raw_val, dict):
        return {str(k): float(v) for k, v in raw_val.items()}
    if isinstance(raw_val, list):
        out: Dict[str, float] = {}
        for item in raw_val:
            try:
                m = str(item.get("motif"))
                w = float(item.get("weight"))
                if m:
                    out[m] = w
            except Exception:
                continue
        return out
    return {}


def load_levels_config(path: Path) -> Dict[int, LevelConfig]:
    """
    Load a JSON list of level configs into a dict[level → LevelConfig].
    """
    raw = json.loads(path.read_text(encoding="utf-8"))
    levels: Dict[int, LevelConfig] = {}
    for entry in raw:
        cfg = LevelConfig(
            level               = entry["level"],
            fragment_min_size   = entry["fragment_min_size"],
            fragment_max_size   = entry["fragment_max_size"],
            min_children        = entry["min_children"],
            max_children        = entry["max_children"],
            overlap_min_size    = entry["overlap_min_size"],
            overlap_max_size    = entry["overlap_max_size"],
            overlap_tm_min      = entry.get("overlap_tm_min"),
            overlap_tm_max      = entry.get("overlap_tm_max"),
            overlap_gc_min      = entry.get("overlap_gc_min"),
            overlap_gc_max      = entry.get("overlap_gc_max"),
            overlap_run_min     = entry.get("overlap_run_min"),
            overlap_run_max     = entry.get("overlap_run_max"),
            overlap_disallowed_motifs = entry.get("overlap_disallowed_motifs", []),
            overlap_allowed_motifs    = _normalize_allowed_motifs(entry.get("overlap_allowed_motifs")),
            max_gap_allowed     = int(entry.get("max_gap_allowed", 0))
        )
        levels[cfg.level] = cfg
    return levels
