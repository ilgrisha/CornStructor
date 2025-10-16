# File: backend/app/services/run_ops.py
# Version: v0.1.0
"""
Run administration operations (delete, cleanup).

This module provides a transactional helper to delete a Run by job_id and,
optionally, its linked Design, as well as the report artifacts folder under
OUTPUT_DIR/<job_id>.

We delete both the Run and its Design to make it easy for users to run
experiments repeatedly without accumulating stale rows or artifacts.
"""
from __future__ import annotations

import shutil
from pathlib import Path
from typing import Optional

from sqlalchemy import select
from sqlalchemy.orm import Session

from backend.app.core.config import settings
from backend.app.db.models import Run, Design


def delete_run_and_design(
    db: Session,
    *,
    job_id: str,
    delete_reports: bool = True,
) -> bool:
    """Delete a run (by job_id), its linked design (if any), and optional artifacts.

    Returns
    -------
    bool
        True if a Run with this job_id existed and was deleted; False if not found.
    """
    run: Optional[Run] = db.execute(select(Run).where(Run.job_id == job_id)).scalar_one_or_none()
    if not run:
        return False

    # Hold a reference to linked design (1:1) before deletion
    design: Optional[Design] = run.design

    # Stage deletions
    db.delete(run)
    if design is not None:
        db.delete(design)

    db.commit()

    # Best-effort artifacts cleanup (does not affect DB tx)
    if delete_reports:
        outdir = Path(settings.OUTPUT_DIR) / job_id
        try:
            shutil.rmtree(outdir, ignore_errors=True)
        except Exception:
            # Ignore filesystem errors; we consider DB deletion successful
            pass

    return True
