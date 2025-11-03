# File: backend/app/config/config_primers.py
# Version: v0.2.0
"""
Primer parameters configuration loader/saver (strict, no legacy compatibility).

- Reads defaults from: backend/app/config/primers_param_default.json
- Reads/writes current from: backend/app/config/primers_param.json
- Validates payloads with PrimerDesignParameters (Pydantic) from core/primer/parameters.py

Usage:
    from backend.app.config.config_primers import load_current_params, save_current_params

Notes
-----
- This loader expects the new JSON schema (camelCase keys), for example:

  {
    "primerLengthMin": 18,
    "primerLengthMax": 25,
    "primerTmMin": 55.0,
    "primerTmMax": 62.0,
    "primerGCMin": 40.0,
    "primerGCMax": 65.0,
    "primerHomopolymerMax": 4,
    "primerThreePrimeEndLength": 5,
    "primerThreePrimePrimerMatchMax": 2,
    "primerTargetMatchMax": 5,
    "primerTargetMatchNumberMax": 2,
    "primerSecondaryStructureDeltaGMin": -9.0,
    "primerTmDifferenceMax": 3.0,
    "weights": { "wTm": 1.0, "wGC": 0.5, "wDimer": 1.0, "w3p": 1.25, "wOff": 1.0 },
    "randomSeed": null,
    "templateExact": true
  }

Thread-safety:
- Uses atomic writes (tmp + replace) to avoid partial/dirty writes.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Tuple

from backend.app.core.primer.parameters import PrimerDesignParameters

# Resolve config directory relative to this file
_THIS_DIR = Path(__file__).resolve().parent
CONFIG_DIR = _THIS_DIR
DEFAULT_FILE = CONFIG_DIR / "primers_param_default.json"
CURRENT_FILE = CONFIG_DIR / "primers_param.json"


def _read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def _atomic_write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
        f.write("\n")
    os.replace(tmp, path)


def load_default_params() -> PrimerDesignParameters:
    """Load default primer design parameters from primers_param_default.json."""
    payload = _read_json(DEFAULT_FILE)
    return PrimerDesignParameters.model_validate(payload or {})


def load_current_params(fallback_to_default: bool = True) -> PrimerDesignParameters:
    """
    Load current (editable) primer design parameters.
    If file missing/invalid and fallback is True, return defaults.
    """
    payload = _read_json(CURRENT_FILE)
    if not payload and fallback_to_default:
        return load_default_params()
    return PrimerDesignParameters.model_validate(payload or {})


def save_current_params(params: PrimerDesignParameters) -> None:
    """Persist current parameters to primers_param.json (atomic write)."""
    _atomic_write_json(CURRENT_FILE, params.model_dump())


def ensure_current_exists() -> Tuple[bool, PrimerDesignParameters]:
    """
    Ensure primers_param.json exists; if not, initialize from defaults.
    Returns (created, params).
    """
    if CURRENT_FILE.exists():
        return False, load_current_params()
    defaults = load_default_params()
    save_current_params(defaults)
    return True, defaults
