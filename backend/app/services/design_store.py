# File: backend/app/services/design_store.py
# Version: v0.4.1
"""Service-layer helpers for Design persistence and run linkage.

v0.4.1:
- Clarify that `params_json` stores snapshots under keys 'globals' and 'levels'.
"""
from __future__ import annotations

import json
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
        params_json=params_json,  # may later be enriched with 'globals'/'levels' snapshots
    )
    run.design = design

    db.add(run)
    db.flush()
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


def save_ga_progress(db: Session, *, job_id: str, ga_progress_json: str) -> bool:
    """Persist GA progress JSON into the Design linked to a Run."""
    run = db.execute(select(Run).where(Run.job_id == job_id)).scalar_one_or_none()
    if not run or not run.design_id:
        return False
    design = db.get(Design, run.design_id)
    if design is None:
        return False
    design.ga_progress_json = ga_progress_json
    db.add(design)
    db.commit()
    return True


def save_config_snapshot(
    db: Session, *, job_id: str, globals_json: str, levels_json: str
) -> bool:
    """
    Merge snapshots of globals/levels into Design.params_json for deterministic
    re-rendering later.

    Stored structure:
        {
          "globals": <dict | JSON>,
          "levels": <dict | JSON>,
          "analysisParams": <existing analysis params if any (dict/JSON)>
        }
    """
    run = db.execute(select(Run).where(Run.job_id == job_id)).scalar_one_or_none()
    if not run or not run.design_id:
        return False
    design = db.get(Design, run.design_id)
    if not design:
        return False

    # base
    try:
        current = json.loads(design.params_json) if design.params_json else {}
    except Exception:
        current = {}

    analysis_params = current.get("analysisParams", current if isinstance(current, dict) else {})

    try:
        gl = json.loads(globals_json)
    except Exception:
        gl = {}
    try:
        lv = json.loads(levels_json)
    except Exception:
        lv = {}

    merged = {
        "globals": gl,
        "levels": lv,
        "analysisParams": analysis_params,
    }
    design.params_json = json.dumps(merged)
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
