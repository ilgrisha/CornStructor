# File: backend/app/api/v1/api.py
# Version: v0.5.0
"""
v1 API aggregator.

Routers included:
- health
- design
- designs (persisted design records)
- runs
- secondary_structure (ViennaRNA-backed stems analysis)
- live (on-the-fly reports)
"""
from __future__ import annotations

from fastapi import APIRouter

from . import health as health_router
from . import design as design_router
from . import designs as designs_router
from . import runs as runs_router
from . import secondary_structure as stems_router
from . import live_reports as live_router

api_router = APIRouter()

api_router.include_router(health_router.router)
api_router.include_router(design_router.router)
api_router.include_router(designs_router.router)
api_router.include_router(runs_router.router)
api_router.include_router(stems_router.router)
api_router.include_router(live_router.router)
