# File: backend/app/api/v1/reports_dynamic.py
# Version: v1.1.0
"""
Dynamic reports router (no disk persistence required).

Serves /reports/{job_id}/... by rendering artifacts on the fly from database
(Design.sequence, Design.params_json, Design.tree_json, Design.ga_progress_json).

New in v1.1.0
-------------
- Add endpoints to serve stored config snapshots directly:
    GET /reports/{job_id}/levels.snapshot.json
    GET /reports/{job_id}/globals.snapshot.json

These read snapshots from Design.params_json if present, otherwise fall back to
the current config files on disk.
"""
from __future__ import annotations

import json
from pathlib import Path
from fastapi import APIRouter, Depends, HTTPException
from fastapi.responses import FileResponse, JSONResponse
from sqlalchemy.orm import Session

from backend.app.db.session import get_db
from backend.app.services.design_store import get_design_by_run
from backend.app.core.reporting.on_the_fly import ensure_rendered, RENDERED_INDEX
from backend.app.core.config import settings

router = APIRouter(prefix="/reports", tags=["reports (dynamic)"])


def _ensure_and_read(job_id: str, relpath: str, db: Session) -> Path:
    design = get_design_by_run(db, job_id)
    if not design:
        raise HTTPException(status_code=404, detail="Design not found for run")
    outdir = ensure_rendered(job_id=job_id, design=design)
    target = outdir / relpath
    if not target.is_file():
        raise HTTPException(status_code=404, detail="Not Found")
    return target


@router.get("/{job_id}/index.html")
def report_index(job_id: str, db: Session = Depends(get_db)):
    path = _ensure_and_read(job_id, RENDERED_INDEX, db)
    return FileResponse(str(path), media_type="text/html; charset=utf-8")


@router.get("/{job_id}/analysis.html")
def report_analysis_html(job_id: str, db: Session = Depends(get_db)):
    path = _ensure_and_read(job_id, "analysis.html", db)
    return FileResponse(str(path), media_type="text/html; charset=utf-8")


@router.get("/{job_id}/tree.html")
def report_tree_html(job_id: str, db: Session = Depends(get_db)):
    path = _ensure_and_read(job_id, "tree.html", db)
    return FileResponse(str(path), media_type="text/html; charset=utf-8")


@router.get("/{job_id}/fragments.fasta")
def report_fragments(job_id: str, db: Session = Depends(get_db)):
    path = _ensure_and_read(job_id, "fragments.fasta", db)
    return FileResponse(str(path), media_type="text/plain; charset=utf-8")


@router.get("/{job_id}/ga_progress.html")
def report_ga_progress_html(job_id: str, db: Session = Depends(get_db)):
    path = _ensure_and_read(job_id, "ga_progress.html", db)
    return FileResponse(str(path), media_type="text/html; charset=utf-8")


@router.get("/{job_id}/ga_progress.json")
def report_ga_progress_json(job_id: str, db: Session = Depends(get_db)):
    path = _ensure_and_read(job_id, "ga_progress.json", db)
    return FileResponse(str(path), media_type="application/json")


@router.get("/{job_id}/tree.json")
def report_tree_json(job_id: str, db: Session = Depends(get_db)):
    path = _ensure_and_read(job_id, "tree.json", db)
    return FileResponse(str(path), media_type="application/json")


@router.get("/{job_id}/clusters/{name}")
def report_cluster_html(job_id: str, name: str, db: Session = Depends(get_db)):
    # Basic path safety
    if "/" in name or "\\" in name:
        raise HTTPException(status_code=400, detail="Invalid path")
    path = _ensure_and_read(job_id, f"clusters/{name}", db)
    return FileResponse(str(path), media_type="text/html; charset=utf-8")


# ---------------- NEW: expose stored snapshots directly ---------------- #

def _load_snapshot_from_design(design_params_json: str | None, key: str) -> dict | None:
    if not design_params_json:
        return None
    try:
        data = json.loads(design_params_json)
        snap = data.get(key)
        if isinstance(snap, dict):
            return snap
        # allow array or other JSON types if your config uses them
        if snap is not None:
            return snap
    except Exception:
        pass
    return None


def _load_json_file(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


@router.get("/{job_id}/levels.snapshot.json")
def levels_snapshot(job_id: str, db: Session = Depends(get_db)):
    """
    Return the exact levels config used for this run (snapshot).
    Falls back to current settings.LEVELS_PATH if snapshot was not saved.
    """
    design = get_design_by_run(db, job_id)
    if not design:
        raise HTTPException(status_code=404, detail="Design not found for run")

    snap = _load_snapshot_from_design(design.params_json, "levels")
    if snap is None:
        # fallback to current file on disk
        snap = _load_json_file(settings.LEVELS_PATH)

    return JSONResponse(snap)


@router.get("/{job_id}/globals.snapshot.json")
def globals_snapshot(job_id: str, db: Session = Depends(get_db)):
    """
    Return the exact global config used for this run (snapshot).
    Falls back to current settings.GLOBALS_PATH if snapshot was not saved.
    """
    design = get_design_by_run(db, job_id)
    if not design:
        raise HTTPException(status_code=404, detail="Design not found for run")

    snap = _load_snapshot_from_design(design.params_json, "globals")
    if snap is None:
        # fallback to current file on disk
        snap = _load_json_file(Path(settings.GLOBALS_PATH))

    return JSONResponse(snap)
