# File: backend/app/services/run_store.py
# Version: v0.1.0
"""
Run persistence helpers (service layer).

These functions encapsulate the DB logic so callers (routers, job runner) don't
need to import SQLAlchemy session management details.
"""
from __future__ import annotations

from typing import Optional

from sqlalchemy import select, delete, func
from sqlalchemy.orm import Session

from backend.app.db.models import Run, RunStatus

def record_run_start(db: Session, *, job_id: str, sequence_len: Optional[int] = None, params_json: Optional[str] = None, note: Optional[str] = None) -> Run:
    """Create or update a run row to RUNNING."""
    run = db.execute(select(Run).where(Run.job_id == job_id)).scalar_one_or_none()
    if run is None:
        run = Run(job_id=job_id)
        db.add(run)
    run.status = RunStatus.RUNNING
    run.sequence_len = sequence_len
    run.params_json = params_json
    run.note = note
    run.report_url = None
    run.exit_code = None
    run.touch()
    db.commit()
    db.refresh(run)
    return run

def record_run_completion(db: Session, *, job_id: str, report_url: Optional[str], exit_code: int) -> Run | None:
    """Mark a run as COMPLETED or FAILED and set report_url/exit_code."""
    run = db.execute(select(Run).where(Run.job_id == job_id)).scalar_one_or_none()
    if run is None:
        return None
    run.status = RunStatus.COMPLETED if exit_code == 0 else RunStatus.FAILED
    run.report_url = report_url
    run.exit_code = exit_code
    run.touch()
    db.commit()
    db.refresh(run)
    return run

def list_runs(db: Session, *, q: Optional[str] = None, limit: int = 50, offset: int = 0) -> tuple[int, list[Run]]:
    """Return (total, items) filtered by optional free-text q against job_id or note."""
    stmt = select(Run).order_by(Run.created_at.desc()).limit(limit).offset(offset)
    count_stmt = select(func.count()).select_from(Run)
    if q:
        like = f"%{q}%"
        stmt = stmt.where((Run.job_id.ilike(like)) | (Run.note.ilike(like)))
        count_stmt = count_stmt.where((Run.job_id.ilike(like)) | (Run.note.ilike(like)))
    total = db.execute(count_stmt).scalar_one()
    items = list(db.execute(stmt).scalars())
    return total, items

def get_run(db: Session, job_id: str) -> Run | None:
    return db.execute(select(Run).where(Run.job_id == job_id)).scalar_one_or_none()

def delete_run(db: Session, job_id: str) -> bool:
    res = db.execute(delete(Run).where(Run.job_id == job_id))
    db.commit()
    return res.rowcount > 0
