# File: backend/app/api/v1/reports_dynamic.py
# Version: v1.3.0
"""
Dynamic reports router (no disk persistence required).

v1.3.0:
- Add routes for input.fasta and assembly.gb (GenBank).
- Keep explicit routes for globals.json and levels.json (downloadable).
"""
from __future__ import annotations

from pathlib import Path
from fastapi import APIRouter, Depends, HTTPException
from fastapi.responses import FileResponse
from sqlalchemy.orm import Session

from backend.app.db.session import get_db
from backend.app.services.design_store import get_design_by_run
from backend.app.core.reporting.on_the_fly import ensure_rendered, RENDERED_INDEX

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


@router.get("/{job_id}/input.fasta")
def report_input_fasta(job_id: str, db: Session = Depends(get_db)):
    path = _ensure_and_read(job_id, "input.fasta", db)
    return FileResponse(str(path), media_type="text/plain; charset=utf-8")


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


@router.get("/{job_id}/oligos.fasta")
def report_oligos(job_id: str, db: Session = Depends(get_db)):
    path = _ensure_and_read(job_id, "oligos.fasta", db)
    return FileResponse(str(path), media_type="text/plain; charset=utf-8")


@router.get("/{job_id}/assembly.gb")
def report_genbank(job_id: str, db: Session = Depends(get_db)):
    path = _ensure_and_read(job_id, "assembly.gb", db)
    return FileResponse(str(path), media_type="application/octet-stream", filename="assembly.gb")


@router.get("/{job_id}/fragments.csv")
def report_fragments_csv(job_id: str, db: Session = Depends(get_db)):
    path = _ensure_and_read(job_id, "fragments.csv", db)
    return FileResponse(str(path), media_type="text/csv; charset=utf-8", filename="fragments.csv")


@router.get("/{job_id}/oligos.csv")
def report_oligos_csv(job_id: str, db: Session = Depends(get_db)):
    path = _ensure_and_read(job_id, "oligos.csv", db)
    return FileResponse(str(path), media_type="text/csv; charset=utf-8", filename="oligos.csv")


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


@router.get("/{job_id}/globals.json")
def globals_json(job_id: str, db: Session = Depends(get_db)):
    path = _ensure_and_read(job_id, "globals.json", db)
    return FileResponse(str(path), media_type="application/json")


@router.get("/{job_id}/levels.json")
def levels_json(job_id: str, db: Session = Depends(get_db)):
    path = _ensure_and_read(job_id, "levels.json", db)
    return FileResponse(str(path), media_type="application/json")


@router.get("/{job_id}/clusters/{name}")
def report_cluster_html(job_id: str, name: str, db: Session = Depends(get_db)):
    # Basic path safety
    if "/" in name or "\\" in name:
        raise HTTPException(status_code=400, detail="Invalid path")
    path = _ensure_and_read(job_id, f"clusters/{name}", db)
    return FileResponse(str(path), media_type="text/html; charset=utf-8")
