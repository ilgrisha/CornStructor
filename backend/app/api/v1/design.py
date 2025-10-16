# File: backend/app/api/v1/design.py
# Version: v0.4.1
"""
Design API:
- POST /design/start -> { job_id, jobId }
- GET  /design/{job_id}/logs -> SSE stream

v0.4.1:
- FIX: Do not pass unsupported `params=` kwarg to JobManager.start_with_existing_id().
- Keep behavior from v0.4.0:
  * Merge/lock incoming tree parameters (Globals & Levels) against backend defaults.
  * Write merged "current" parameters to backend/app/config/globals.json and levels.json
    so the construction pipeline uses the latest edits.
  * Persist merged params on the Design record.
"""
from __future__ import annotations

import asyncio
import json
from pathlib import Path
from typing import Any, AsyncGenerator, Dict, List, Tuple

from fastapi import APIRouter, HTTPException, Depends
from fastapi.responses import StreamingResponse
from pydantic import BaseModel, Field
from sqlalchemy.orm import Session

from backend.app.core.jobs.runner import job_manager
from backend.app.db.session import get_db
from backend.app.services.design_store import create_design_with_run

router = APIRouter(prefix="/design", tags=["design"])


def _sse(data: str) -> bytes:
    return f"data: {data}\n\n".encode("utf-8")


# --------- Models ---------
class DesignStartRequest(BaseModel):
    """Payload to kick off a run."""
    sequence: str = Field(..., min_length=10, description="DNA sequence (A/C/G/T/N)")
    params: Dict[str, Any] | None = None
    toggles: Dict[str, bool] | None = None


class DesignStartResponse(BaseModel):
    """Return both snake_case and camelCase for compatibility with existing UIs."""
    job_id: str
    jobId: str


# --------- Config helpers (locked schema + IO) ---------
CFG_DIR = Path("backend/app/config")
DEFAULT_GLOBALS = CFG_DIR / "globals_param_default.json"
DEFAULT_LEVELS = CFG_DIR / "levels_param_default.json"
CURRENT_GLOBALS = CFG_DIR / "globals.json"   # <- actively used by pipeline
CURRENT_LEVELS = CFG_DIR / "levels.json"     # <- actively used by pipeline


def _load_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        raise HTTPException(status_code=500, detail=f"Required config file missing: {path}")
    except json.JSONDecodeError as e:
        raise HTTPException(status_code=500, detail=f"Invalid JSON in {path.name}: {e}")


def _write_json(path: Path, data: Any) -> None:
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")


def _lock_object(template: Dict[str, Any], incoming: Dict[str, Any] | None) -> Dict[str, Any]:
    """
    Enforce locked schema for a dict object based on a template:
    - Only keys from template are allowed, in the same order.
    - Values come from incoming if present; otherwise template defaults are used.
    """
    out: Dict[str, Any] = {}
    incoming = incoming or {}
    for k in template.keys():
        out[k] = incoming[k] if k in incoming else template[k]
    return out


def _lock_levels(template_levels: List[Dict[str, Any]], incoming_levels: List[Dict[str, Any]] | None) -> List[Dict[str, Any]]:
    """
    Enforce locked schema for the Levels array.
    - Each level row is locked against the template row shape (use first template object).
    - Index reflow of the `level` field is applied if present.
    """
    if not template_levels:
        return []
    tmpl = template_levels[0]
    rows_in = incoming_levels or template_levels
    locked: List[Dict[str, Any]] = []
    for i, row in enumerate(rows_in):
        locked_row = _lock_object(tmpl, row or {})
        if "level" in locked_row:
            locked_row["level"] = i
        locked.append(locked_row)
    return locked


def _merge_and_lock_tree(incoming_params: Dict[str, Any] | None) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    Returns (merged_globals, merged_levels_obj) with locked schema.
    - Reads defaults from DEFAULT_GLOBALS/DEFAULT_LEVELS
    - Applies incoming overrides from params.tree.globals / params.tree.levels
    - Writes nothing here; just returns merged structures.
    """
    defaults_globals = _load_json(DEFAULT_GLOBALS)
    defaults_levels = _load_json(DEFAULT_LEVELS)

    if not isinstance(defaults_globals, dict) or not isinstance(defaults_levels, list):
        raise HTTPException(status_code=500, detail="Defaults must be object (globals) and array (levels)")

    tree = (incoming_params or {}).get("tree") or {}
    inc_globals = tree.get("globals")
    inc_levels = tree.get("levels")

    merged_globals = _lock_object(defaults_globals, inc_globals if isinstance(inc_globals, dict) else None)
    merged_levels_list = _lock_levels(defaults_levels, inc_levels if isinstance(inc_levels, list) else None)

    # Return both the merged globals and a dict containing levels (for symmetry/persistence)
    return merged_globals, {"levels": merged_levels_list}


@router.post("/start", response_model=DesignStartResponse)
async def start_design(req: DesignStartRequest, db: Session = Depends(get_db)):
    """
    Start a construction tree design job:

    - Merge/lock any provided tree.globals and tree.levels against backend defaults.
    - Write the merged to CURRENT_GLOBALS and CURRENT_LEVELS (so pipeline uses latest).
    - Persist merged params on the Design (snake+camel casing of job id returned).
    - Trigger the async job runner.
    """
    # Reserve a job id so we can atomically persist Run + Design
    job_id = job_manager.reserve_job_id()

    # Merge/lock tree params against defaults
    merged_globals, merged_levels_obj = _merge_and_lock_tree(req.params or {})

    # Build the final params object we persist:
    final_params: Dict[str, Any] = dict(req.params or {})
    tree = dict(final_params.get("tree") or {})
    tree["globals"] = merged_globals
    tree["levels"] = merged_levels_obj["levels"]
    final_params["tree"] = tree

    # Persist Run + Design with final (merged) params
    params_json = json.dumps(final_params)
    create_design_with_run(db, job_id=job_id, sequence=req.sequence, params_json=params_json)

    # Write "current" parameters to the active config files used by the pipeline
    _write_json(CURRENT_GLOBALS, merged_globals)
    _write_json(CURRENT_LEVELS, tree["levels"])

    # Kick off the async pipeline (runner no longer receives `params` kwarg)
    await job_manager.start_with_existing_id(job_id=job_id, fasta_text=req.sequence)

    return DesignStartResponse(job_id=job_id, jobId=job_id)


@router.get("/{job_id}/logs")
async def stream_logs(job_id: str):
    job = job_manager.get(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Unknown job_id")

    async def gen() -> AsyncGenerator[bytes, None]:
        while True:
            item = await job.queue.get()
            if item == "__EOF__":
                yield _sse("Log stream completed")
                break
            yield _sse(str(item))
            await asyncio.sleep(0.01)

    return StreamingResponse(gen(), media_type="text/event-stream")
