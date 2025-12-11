# File: backend/app/api/v1/designs.py
# Version: v0.1.0
"""Design read APIs.

Endpoints:
- GET /designs/by-run/{job_id} -> design (sequence, params, tree_json)
- GET /designs/{design_id}     -> design by id
"""
from __future__ import annotations

from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session

from backend.app.db.session import get_db
from backend.app.db.schemas.design import DesignDetail, DesignByRunResponse
from backend.app.services.design_store import get_design_by_run
from backend.app.db.models import Design, Run

router = APIRouter(prefix="/designs", tags=["designs"])

@router.get("/by-run/{job_id}", response_model=DesignByRunResponse)
def read_design_by_run(job_id: str, db: Session = Depends(get_db)):
    design = get_design_by_run(db, job_id)
    if not design:
        raise HTTPException(status_code=404, detail="Design not found for run")
    # Join to get job_id back to UI
    run = db.query(Run).filter(Run.design_id == design.id).first()
    return {
        **DesignDetail.model_validate(design).model_dump(),
        "job_id": run.job_id if run else "",
        "run_note": run.note if run else None,
    }

@router.get("/{design_id}", response_model=DesignDetail)
def read_design(design_id: int, db: Session = Depends(get_db)):
    design = db.get(Design, design_id)
    if not design:
        raise HTTPException(status_code=404, detail="Design not found")
    return design
