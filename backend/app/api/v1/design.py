# File: backend/app/api/v1/design.py
# Version: v0.3.0
"""
Design API:
- POST /design/start -> { job_id, jobId }
- GET  /design/{job_id}/logs -> SSE stream

v0.3.0:
- On start, create a Design linked to the Run capturing sequence and params.
"""
from __future__ import annotations

import asyncio
from typing import Any, AsyncGenerator

from fastapi import APIRouter, HTTPException, Depends
from fastapi.responses import StreamingResponse
from pydantic import BaseModel, Field
from sqlalchemy.orm import Session
import json

from backend.app.core.jobs.runner import job_manager
from backend.app.db.session import get_db
from backend.app.services.design_store import create_design_with_run

router = APIRouter(prefix="/design", tags=["design"])

def _sse(data: str) -> bytes:
    return f"data: {data}\n\n".encode("utf-8")

class DesignStartRequest(BaseModel):
    """Payload to kick off a run."""
    sequence: str = Field(..., min_length=10, description="DNA sequence (A/C/G/T/N)")
    params: dict[str, Any] | None = None
    toggles: dict[str, bool] | None = None

class DesignStartResponse(BaseModel):
    """Return both snake_case and camelCase for compatibility with existing UIs."""
    job_id: str
    jobId: str

@router.post("/start", response_model=DesignStartResponse)
async def start_design(req: DesignStartRequest, db: Session = Depends(get_db)):
    # Reserve a job id so we can atomically persist Run + Design
    job_id = job_manager.reserve_job_id()

    # Persist Run + Design
    params_json = json.dumps(req.params) if req.params is not None else None
    create_design_with_run(db, job_id=job_id, sequence=req.sequence, params_json=params_json)

    # Kick off the async pipeline
    await job_manager.start_with_existing_id(job_id=job_id, fasta_text=req.sequence)

    return DesignStartResponse(job_id=job_id, jobId=job_id)

@router.get("/{job_id}/logs")
async def stream_logs(job_id: str):
    job = job_manager.get(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Unknown job_id")

    async def gen() -> AsyncGenerator[bytes, None]:
        while True:
            item = await job.queue.get()
            if item == "__EOF__":
                yield _sse("Log stream completed")
                break
            yield _sse(str(item))
            await asyncio.sleep(0.01)

    return StreamingResponse(gen(), media_type="text/event-stream")
