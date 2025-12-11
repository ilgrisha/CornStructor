# File: backend/app/db/schemas/design.py
# Version: v0.2.0
"""Pydantic schemas for Design resources (includes GA progress)."""
from __future__ import annotations

from datetime import datetime
from pydantic import BaseModel, Field
from typing import Optional

class DesignBase(BaseModel):
    id: int
    created_at: datetime
    updated_at: datetime
    sequence: str = Field(..., description="DNA sequence (A/C/G/T/N)")
    sequence_len: int
    params_json: Optional[str] = Field(None, description="Creation parameters JSON (globals + levels)")
    tree_json: Optional[str] = Field(None, description="Serialized construction tree JSON")
    ga_progress_json: Optional[str] = Field(None, description="GA progress timeline JSON")
    name: Optional[str] = None

    class Config:
        from_attributes = True

class DesignCreateRequest(BaseModel):
    sequence: str
    params_json: Optional[str] = None
    name: Optional[str] = None

class DesignDetail(DesignBase):
    pass

class DesignByRunResponse(DesignBase):
    job_id: str
    run_note: Optional[str] = None
