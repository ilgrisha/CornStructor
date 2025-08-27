# File: backend/app/api/v1/design.py
# Version: v0.2.1
"""
Design API:
- POST /design/start -> { job_id, jobId }
- GET  /design/{job_id}/logs -> SSE stream
"""
from __future__ import annotations

import asyncio
from typing import Any, AsyncGenerator

from fastapi import APIRouter, HTTPException
from fastapi.responses import StreamingResponse
from pydantic import BaseModel, Field

from backend.app.core.jobs.runner import job_manager

router = APIRouter(prefix="/design", tags=["design"])


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
async def start_design(req: DesignStartRequest) -> DesignStartResponse:
    bad = [c for c in req.sequence.upper() if c not in {"A", "C", "G", "T", "N"}]
    if bad:
        raise HTTPException(status_code=400, detail=f"Invalid characters in sequence: {sorted(set(bad))}")
    job_id = await job_manager.start(req.sequence)
    # Return both styles so either `job_id` or `jobId` works on the frontend
    return DesignStartResponse(job_id=job_id, jobId=job_id)


def _sse(data: str) -> bytes:
    return f"data: {data}\n\n".encode("utf-8")


@router.get("/{job_id}/logs")
async def stream_logs(job_id: str) -> StreamingResponse:
    """Attach to a job's log queue and stream as text/event-stream."""
    if job_id == "undefined" or not job_id:
        raise HTTPException(status_code=400, detail="Missing or invalid job_id")
    job = job_manager.get(job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Unknown job_id")

    async def gen() -> AsyncGenerator[bytes, None]:
        while True:
            item = await job.queue.get()
            if item == "__EOF__":
                yield _sse("Log stream completed")
                break
            yield _sse(item)
            await asyncio.sleep(0.01)

    return StreamingResponse(gen(), media_type="text/event-stream")
