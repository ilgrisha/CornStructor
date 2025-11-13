# File: backend/app/api/v1/reports_dynamic.py
# Version: v1.4.0
"""
Dynamic reports router (no disk persistence required).

v1.4.0:
- Add route to download a ZIP bundle containing *all* rendered report artifacts:
    GET /reports/{job_id}/bundle.zip

v1.3.0:
- Add routes for input.fasta and assembly.gb (GenBank).
- Keep explicit routes for globals.json and levels.json (downloadable).
"""
from __future__ import annotations

import io
import zipfile
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


def _ensure_and_get_outdir(job_id: str, db: Session) -> Path:
    """Render (if needed) and return the run's output directory."""
    design = get_design_by_run(db, job_id)
    if not design:
        raise HTTPException(status_code=404, detail="Design not found for run")
    return ensure_rendered(job_id=job_id, design=design)


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


# -----------------------
# NEW: ZIP bundle download
# -----------------------

@router.get("/{job_id}/bundle.zip")
def report_bundle_zip(job_id: str, db: Session = Depends(get_db)):
    """
    Create and return a ZIP of the entire rendered report directory.
    This is generated on-the-fly each time to guarantee freshness.
    """
    outdir = _ensure_and_get_outdir(job_id, db)
    zip_path = outdir / "bundle.zip"

    # Build zip in-memory then write to disk atomically to avoid partial reads
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for p in outdir.rglob("*"):
            if not p.is_file():
                continue
            # Exclude the zip file itself if present
            if p.name == "bundle.zip":
                continue
            # Keep paths inside the zip relative to the outdir root
            arcname = p.relative_to(outdir).as_posix()
            zf.write(p, arcname)
    buf.seek(0)
    zip_path.write_bytes(buf.read())

    # Offer a nice filename in the download prompt
    filename = f"{job_id}_report_bundle.zip"
    return FileResponse(str(zip_path), media_type="application/zip", filename=filename)
