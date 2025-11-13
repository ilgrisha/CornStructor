# File: backend/app/api/v1/params.py
# Version: v0.1.0
"""
Parameters API for construction tree defaults (Globals & Levels).

Endpoints:
- GET /params/defaults -> { globals: {...}, levels: [...] }

Notes:
- Reads from backend/app/config/globals_param_default.json and
  backend/app/config/levels_param_default.json.
- Schema is *locked*: keys are dictated by the JSON files; clients should only edit values.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List

from fastapi import APIRouter, HTTPException

router = APIRouter(prefix="/params", tags=["params"])

CFG_DIR = Path("backend/app/config")
GLOBALS_PATH = CFG_DIR / "globals_param_default.json"
LEVELS_PATH  = CFG_DIR / "levels_param_default.json"


@router.get("/defaults")
def get_defaults() -> Dict[str, Any]:
    """Return default Globals and Levels as loaded from JSON files."""
    try:
        g = json.loads(GLOBALS_PATH.read_text(encoding="utf-8"))
        lv = json.loads(LEVELS_PATH.read_text(encoding="utf-8"))
    except FileNotFoundError as e:
        raise HTTPException(status_code=500, detail=f"Defaults file missing: {e}")
    except json.JSONDecodeError as e:
        raise HTTPException(status_code=500, detail=f"Invalid defaults JSON: {e}")

    if not isinstance(g, dict) or not isinstance(lv, list):
        raise HTTPException(status_code=500, detail="Defaults must be object (globals) and array (levels)")

    # Do not allow arbitrary keys in clients; backend defines schema via these files.
    return {"globals": g, "levels": lv}
