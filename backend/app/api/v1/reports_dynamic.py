# File: backend/app/api/v1/reports_dynamic.py
# Version: v1.0.1
"""
Dynamic reports router (no disk persistence required).

Serves /reports/{job_id}/... by rendering artifacts on the fly from database
(Design.sequence, Design.params_json, Design.tree_json, Design.ga_progress_json).

Endpoints (examples):
- GET /reports/{job_id}/index.html
- GET /reports/{job_id}/analysis.html
- GET /reports/{job_id}/ga_progress.html
- GET /reports/{job_id}/ga_progress.json
- GET /reports/{job_id}/tree.json
- GET /reports/{job_id}/tree.html
- GET /reports/{job_id}/fragments.fasta
- GET /reports/{job_id}/clusters/<file>.html
"""
from __future__ import annotations

from pathlib import Path
from fastapi import APIRouter, Depends, HTTPException
from fastapi.responses import FileResponse
from sqlalchemy.orm import Session

from backend.app.db.session import get_db
from backend.app.services.design_store import get_design_by_run
from backend.app.core.reporting.on_the_fly import ensure_rendered, RENDERED_INDEX

# NOTE: prefix='/reports' ensures these routes live at the root /reports/*,
# even when this router is included from outside the /api aggregator.
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
