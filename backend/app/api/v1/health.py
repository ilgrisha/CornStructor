# File: backend/app/api/v1/health.py
# Version: v0.1.0
"""
Simple healthcheck router.
"""
from __future__ import annotations
from fastapi import APIRouter

router = APIRouter(tags=["health"])


@router.get("/health")
def health() -> dict[str, str]:
    """Return a minimal health payload."""
    return {"status": "ok"}
