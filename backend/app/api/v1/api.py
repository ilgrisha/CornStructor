# File: backend/app/api/v1/api.py
# Version: v0.4.0
"""
v1 API aggregator.

Collects and exposes all versioned v1 routers through a single `api_router`.
`main.py` mounts this once with the global API prefix from settings.

Routers included:
- health
- design
- designs (persisted design records)
- runs
- secondary_structure (ViennaRNA-backed stems analysis)
"""
from __future__ import annotations

from fastapi import APIRouter

# Import individual routers from the v1 package
from . import health as health_router
from . import design as design_router
from . import designs as designs_router
from . import runs as runs_router
from . import secondary_structure as stems_router

# Unified router for /api/v1
api_router = APIRouter()

# NOTE: No prefix here—`main.py` applies the versioned/api prefix once.
api_router.include_router(health_router.router)
api_router.include_router(design_router.router)
api_router.include_router(designs_router.router)
api_router.include_router(runs_router.router)
api_router.include_router(stems_router.router)
