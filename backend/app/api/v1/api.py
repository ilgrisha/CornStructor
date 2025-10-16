# File: backend/app/api/v1/api.py
# Version: v0.7.1
"""
v1 API aggregator.

Routers included under /api:
- health
- design
- designs (persisted design records)
- runs
- secondary_structure (ViennaRNA-backed stems analysis)
- live (on-the-fly JSON/HTML helpers)
- params (defaults for globals/levels)

Additionally, we expose `public_router` that mounts the dynamic reports
router at /reports/* (root, not under /api), so main.py doesn't need to
import reports_dynamic directly.
"""
from __future__ import annotations

from fastapi import APIRouter

from . import health as health_router
from . import design as design_router
from . import designs as designs_router
from . import runs as runs_router
from . import secondary_structure as stems_router
from . import live_reports as live_router

# Import the dynamic reports router ONLY here (not from main.py)
from . import reports_dynamic as reports_dynamic_router
from . import params as params_router

# All v1 JSON APIs live under /api via api_router
api_router = APIRouter()
api_router.include_router(health_router.router)
api_router.include_router(design_router.router)
api_router.include_router(designs_router.router)
api_router.include_router(runs_router.router)
api_router.include_router(stems_router.router)
api_router.include_router(live_router.router)
api_router.include_router(params_router.router)

# Public (root-level) router to expose /reports/* without prefix /api.
# main.py will include this router at the app root.
public_router = APIRouter()
public_router.include_router(reports_dynamic_router.router)  # exposes /reports/*
