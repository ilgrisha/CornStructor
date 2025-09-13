# File: backend/app/services/design_store.py
# Version: v0.2.1
"""Service-layer helpers for Design persistence and run linkage.

Functions
---------
create_design_with_run(db, job_id, sequence, params_json) -> (Run, Design)
    Create a RUNNING Run and an associated Design in one transaction.

attach_tree_to_design(db, job_id, tree_json) -> bool
    Persist the serialized tree JSON into the Design associated with a Run.

get_design_by_run(db, job_id) -> Optional[Design]
    Fetch the Design associated with a given Run (by job_id).

mark_run_completed(db, job_id, report_url, exit_code) -> None
    Update Run status (COMPLETED/FAILED), report_url and exit_code.
"""
from __future__ import annotations

from typing import Optional, Tuple

from sqlalchemy import select
from sqlalchemy.orm import Session

from backend.app.db.models import Design, Run, RunStatus


def create_design_with_run(
    db: Session,
    *,
    job_id: str,
    sequence: str,
    params_json: Optional[str],
) -> Tuple[Run, Design]:
    """Create a RUNNING Run and an associated Design in a single transaction."""
    run = Run(job_id=job_id, status=RunStatus.RUNNING)
    run.sequence_len = len(sequence)
    run.params_json = params_json

    design = Design(
        sequence=sequence,
        sequence_len=len(sequence),
        params_json=params_json,
    )
    run.design = design

    db.add(run)
    db.flush()  # assign PKs
    db.commit()
    db.refresh(run)
    db.refresh(design)
    return run, design


def attach_tree_to_design(db: Session, *, job_id: str, tree_json: str) -> bool:
    """Attach serialized tree JSON to the Design associated with a Run."""
    run = db.execute(select(Run).where(Run.job_id == job_id)).scalar_one_or_none()
    if not run or not run.design_id:
        return False
    design = db.get(Design, run.design_id)
    if design is None:
        return False
    design.tree_json = tree_json
    db.add(design)
    db.commit()
    return True


def get_design_by_run(db: Session, job_id: str) -> Optional[Design]:
    """Return the Design linked to a Run by job_id, or None."""
    run = db.execute(select(Run).where(Run.job_id == job_id)).scalar_one_or_none()
    if not run or not run.design_id:
        return None
    return db.get(Design, run.design_id)


def mark_run_completed(
    db: Session,
    *,
    job_id: str,
    report_url: Optional[str],
    exit_code: int,
) -> None:
    """Set Run.status, report_url, and exit_code based on pipeline result."""
    run = db.execute(select(Run).where(Run.job_id == job_id)).scalar_one_or_none()
    if not run:
        return
    run.exit_code = exit_code
    run.report_url = report_url
    run.status = RunStatus.COMPLETED if exit_code == 0 else RunStatus.FAILED
    db.add(run)
    db.commit()
