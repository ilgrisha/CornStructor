# File: backend/app/api/v1/live_reports.py
# Version: v0.1.0
"""
Live (on-the-fly) report endpoints using a temporary directory.

Path B (temp-only render, minimal changes):
- We reuse existing HTML exporters where possible, writing to a TemporaryDirectory
  and returning the file bytes as the HTTP response.
- We avoid re-running the assembler. We leverage persisted DB fields:
  - Design.tree_json (already stored)
  - Design.ga_progress_json (added here)
  - (sequence/params_json available for future enhancements)

Endpoints
---------
GET /live/by-run/{job_id}/index.html
    Minimal index page linking to GA progress (live) and JSON artifacts.

GET /live/by-run/{job_id}/ga_progress.html
    Renders GA progress HTML using `export_ga_progress_html` and returns it.

GET /live/by-run/{job_id}/ga_progress.json
GET /live/by-run/{job_id}/tree.json
    Return the persisted JSON payloads directly from DB.
"""
from __future__ import annotations

import json
import tempfile
from fastapi import APIRouter, Depends, HTTPException, Response
from sqlalchemy.orm import Session

from backend.app.db.session import get_db
from backend.app.services.design_store import get_design_by_run
from backend.app.core.visualization.ga_progress_html_report import export_ga_progress_html

router = APIRouter(prefix="/live", tags=["live-reports"])


@router.get("/by-run/{job_id}/index.html", response_class=Response)
def live_index(job_id: str, db: Session = Depends(get_db)) -> Response:
    design = get_design_by_run(db, job_id)
    if not design:
        raise HTTPException(status_code=404, detail="Design not found for run")
    # Build a tiny HTML index that references live endpoints
    has_gaprog = bool(design.ga_progress_json)
    parts = []
    parts.append(f"<!doctype html><meta charset='utf-8'><title>Live Report — {job_id}</title>")
    parts.append("<style>body{font-family:system-ui,-apple-system,Segoe UI,Roboto,sans-serif;margin:2rem} h1{margin:0 0 .5rem 0} ul{margin:.25rem 0 1rem 1.2rem}</style>")
    parts.append(f"<h1>Live Report — {job_id}</h1>")
    parts.append("<section><h2>Artifacts</h2><ul>")
    parts.append(f"<li><a href='/api/live/by-run/{job_id}/tree.json' target='_blank' rel='noopener'>tree.json</a></li>")
    if has_gaprog:
        parts.append(f"<li><a href='/api/live/by-run/{job_id}/ga_progress.json' target='_blank' rel='noopener'>ga_progress.json</a></li>")
        parts.append(f"<li><a href='/api/live/by-run/{job_id}/ga_progress.html' target='_blank' rel='noopener'>GA progress (HTML)</a></li>")
    else:
        parts.append("<li><em>GA progress not available for this run.</em></li>")
    parts.append("</ul></section>")
    parts.append("<div class='back'><a href='/'>&larr; Back to CornStructor</a></div>")
    html = "\n".join(parts)
    return Response(content=html, media_type="text/html; charset=utf-8")


@router.get("/by-run/{job_id}/ga_progress.html", response_class=Response)
def live_ga_progress_html(job_id: str, db: Session = Depends(get_db)) -> Response:
    design = get_design_by_run(db, job_id)
    if not design or not design.ga_progress_json:
        raise HTTPException(status_code=404, detail="GA progress not found for run")
    try:
        data = json.loads(design.ga_progress_json)
    except Exception:
        raise HTTPException(status_code=500, detail="Invalid GA progress JSON in database")
    # Use a temp dir and exporter to produce HTML, then return it
    with tempfile.TemporaryDirectory() as tmp:
        from pathlib import Path
        out = Path(tmp) / "ga_progress.html"
        export_ga_progress_html(data, out)
        html = out.read_text(encoding="utf-8")
    return Response(content=html, media_type="text/html; charset=utf-8")


@router.get("/by-run/{job_id}/ga_progress.json", response_class=Response)
def live_ga_progress_json(job_id: str, db: Session = Depends(get_db)) -> Response:
    design = get_design_by_run(db, job_id)
    if not design or not design.ga_progress_json:
        raise HTTPException(status_code=404, detail="GA progress not found for run")
    return Response(content=design.ga_progress_json, media_type="application/json")


@router.get("/by-run/{job_id}/tree.json", response_class=Response)
def live_tree_json(job_id: str, db: Session = Depends(get_db)) -> Response:
    design = get_design_by_run(db, job_id)
    if not design or not design.tree_json:
        raise HTTPException(status_code=404, detail="Tree JSON not found for run")
    return Response(content=design.tree_json, media_type="application/json")
