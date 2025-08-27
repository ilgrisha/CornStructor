# File: backend/app/api/v1/runs.py
# Version: v0.1.0
"""
Runs API

Endpoints:
- GET  /runs            -> list runs (q, limit, offset)
- GET  /runs/{job_id}   -> get a run
- DELETE /runs/{job_id} -> delete a run (does not remove report files)
"""
from __future__ import annotations

from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session

from backend.app.db.session import get_db
from backend.app.db.schemas.run import RunItem, RunDetail, RunListResponse
from backend.app.services.run_store import list_runs as svc_list_runs, get_run as svc_get_run, delete_run as svc_delete_run

router = APIRouter(prefix="/runs", tags=["runs"])

@router.get("", response_model=RunListResponse)
def list_runs(q: str | None = Query(default=None, description="Filter by job_id or note"),
              limit: int = Query(default=50, ge=1, le=200),
              offset: int = Query(default=0, ge=0),
              db: Session = Depends(get_db)):
    total, items = svc_list_runs(db, q=q, limit=limit, offset=offset)
    return RunListResponse(
        total=total,
        items=[RunItem(
            job_id=r.job_id,
            status=r.status.value,
            created_at=r.created_at,
            updated_at=r.updated_at,
            report_url=r.report_url,
            exit_code=r.exit_code,
            sequence_len=r.sequence_len,
            note=r.note,
        ) for r in items],
    )

@router.get("/{job_id}", response_model=RunDetail)
def get_run(job_id: str, db: Session = Depends(get_db)):
    r = svc_get_run(db, job_id)
    if not r:
        raise HTTPException(status_code=404, detail="Run not found")
    return RunDetail(
        job_id=r.job_id,
        status=r.status.value,
        created_at=r.created_at,
        updated_at=r.updated_at,
        report_url=r.report_url,
        exit_code=r.exit_code,
        sequence_len=r.sequence_len,
        note=r.note,
        params_json=r.params_json,
    )

@router.delete("/{job_id}", status_code=204)
def delete_run(job_id: str, db: Session = Depends(get_db)):
    ok = svc_delete_run(db, job_id)
    if not ok:
        raise HTTPException(status_code=404, detail="Run not found")
    return
