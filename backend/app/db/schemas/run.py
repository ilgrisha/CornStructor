# File: backend/app/db/schemas/run.py
# Version: v0.1.0
"""
Pydantic schemas for Run resources.
"""
from __future__ import annotations

from datetime import datetime
from pydantic import BaseModel, Field
from typing import Optional, Literal

class RunBase(BaseModel):
    job_id: str = Field(..., description="External ID and reports folder name.")
    status: Literal["running", "completed", "failed"]
    created_at: datetime
    updated_at: datetime
    report_url: Optional[str] = None
    exit_code: Optional[int] = None
    sequence_len: Optional[int] = None
    note: Optional[str] = None

class RunItem(RunBase):
    """Compact list item representation."""
    pass

class RunDetail(RunBase):
    """Detailed representation; includes optional params_json."""
    params_json: Optional[str] = None

class RunCreateRequest(BaseModel):
    """Manual creation (rare). Normally created by the job runner."""
    job_id: str
    sequence_len: Optional[int] = None
    params_json: Optional[str] = None
    note: Optional[str] = None

class RunListResponse(BaseModel):
    total: int
    items: list[RunItem]
