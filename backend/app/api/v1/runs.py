# File: backend/app/api/v1/runs.py
# Version: v0.5.0
"""
Runs API.

Updates in v0.5.0
-----------------
- Correct total counting for SQLAlchemy 2.0 (COUNT(*)).
- Return `updated_at` in the list so the UI can react to recent changes.
- Explicit ordering by created_at DESC, id DESC for stable pagination.
- Add `Cache-Control: no-store` to avoid any intermediary caching.

Endpoints
---------
GET    /runs                 List runs (with pagination)
DELETE /runs/{job_id}        Delete a run, its linked design, and (optionally) its artifacts
"""
from __future__ import annotations

from datetime import datetime
from typing import Optional, List

from fastapi import APIRouter, Depends, HTTPException, Query, Response, status
from pydantic import BaseModel
from sqlalchemy import select, desc, func, or_
from sqlalchemy.orm import Session

from backend.app.db.session import get_db
from backend.app.db.models import Run
from backend.app.services.run_ops import delete_run_and_design

router = APIRouter(prefix="/runs", tags=["runs"])


# ---------- Schemas ----------
class RunItem(BaseModel):
    job_id: str
    status: str
    sequence_len: Optional[int] = None
    report_url: Optional[str] = None
    created_at: datetime
    updated_at: datetime
    note: Optional[str] = None

    class Config:
        from_attributes = True


class RunListResponse(BaseModel):
    items: List[RunItem]
    total: int


# ---------- Endpoints ----------
@router.get("", response_model=RunListResponse)
def list_runs(
    response: Response,
    db: Session = Depends(get_db),
    limit: int = Query(50, ge=1, le=500),
    offset: int = Query(0, ge=0),
    q: Optional[str] = Query(None, description="Filter by job_id or description substring"),
):
    """
    Return paginated runs ordered by newest first.

    Notes:
    - Uses COUNT(*) for total (SQLAlchemy 2.0).
    - Ordered by created_at DESC, id DESC for stability across inserts.
    """
    # Base statements
    stmt_base = select(Run)
    count_base = select(func.count()).select_from(Run)

    if q:
        like = f"%{q}%"
        predicate = or_(Run.job_id.ilike(like), Run.note.ilike(like))
        stmt_base = stmt_base.where(predicate)
        count_base = count_base.where(predicate)

    total: int = int(db.scalar(count_base) or 0)

    rows = (
        db.execute(
            stmt_base.order_by(desc(Run.created_at), desc(Run.id)).limit(limit).offset(offset)
        )
        .scalars()
        .all()
    )

    items = [RunItem.model_validate(r) for r in rows]

    # Avoid caching (the UI expects fresh lists)
    response.headers["Cache-Control"] = "no-store"
    return RunListResponse(items=items, total=total)


@router.delete("/{job_id}", status_code=status.HTTP_204_NO_CONTENT)
def delete_run(
    job_id: str,
    db: Session = Depends(get_db),
    delete_reports: bool = Query(True, description="Also delete /reports/<job_id> folder"),
):
    ok = delete_run_and_design(db, job_id=job_id, delete_reports=delete_reports)
    if not ok:
        raise HTTPException(status_code=404, detail="Run not found")
    return Response(status_code=status.HTTP_204_NO_CONTENT)
