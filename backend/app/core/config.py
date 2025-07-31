# File: backend/app/core/config.py
# Version: v0.1.0

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List


@dataclass
class LevelConfig:
    """
    Configuration for one level of the hierarchical assembly tree.
    """
    level: int
    fragment_min_size:   int   # minimum bp length of a fragment at this level
    fragment_max_size:   int   # maximum bp length of a fragment at this level
    min_children:        int   # minimum number of sub‐fragments (children)
    max_children:        int   # maximum number of sub‐fragments
    overlap_min_size:    int   # GA: minimum overlap length at this level
    overlap_max_size:    int   # GA: maximum overlap length at this level


def load_levels_config(path: Path) -> Dict[int, LevelConfig]:
    """
    Load a JSON list of level configs into a dict[level → LevelConfig].

    JSON format should be:
    [
      {
        "level": 0,
        "fragment_min_size": 200,
        "fragment_max_size": 500,
        "min_children": 2,
        "max_children": 4,
        "overlap_min_size": 20,
        "overlap_max_size": 35
      },
      {
        "level": 1,
        "fragment_min_size": 80,
        "fragment_max_size": 200,
        "min_children": 3,
        "max_children": 6,
        "overlap_min_size": 15,
        "overlap_max_size": 25
      },
      …
    ]
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
            overlap_max_size    = entry["overlap_max_size"]
        )
        levels[cfg.level] = cfg
    # Optional: warn if levels are non‐contiguous or missing 0
    return levels
