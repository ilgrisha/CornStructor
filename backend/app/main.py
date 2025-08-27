# File: backend/app/main.py
# Version: v0.3.0
"""
FastAPI application entrypoint.

Mounts:
- /api/health
- /api/design
- /api/runs
- /reports (static) -> serves OUTPUT_DIR for dev parity (Nginx serves this in prod)

CORS enabled per settings.
OpenAPI/Docs exposed at /api/openapi.json and /docs.

Run (local):
  uvicorn backend.app.main:app --reload --port 8000
"""
from __future__ import annotations

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles

from backend.app.core.config import settings
from backend.app.api.v1 import health as health_router
from backend.app.api.v1 import design as design_router
from backend.app.api.v1 import runs as runs_router

# DB init
from backend.app.db.session import engine
from backend.app.db.models import Base

# Create tables at startup if they don't exist (simple bootstrap; Alembic recommended later)
Base.metadata.create_all(bind=engine)

app = FastAPI(title=settings.APP_NAME, version=settings.APP_VERSION)

# CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.cors_origins_list,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Routers
app.include_router(health_router.router, prefix=settings.API_PREFIX)
app.include_router(design_router.router, prefix=settings.API_PREFIX)
app.include_router(runs_router.router, prefix=settings.API_PREFIX)

# Dev/Parity: also serve reports directly from the backend (Nginx handles this in prod)
app.mount(
    "/reports",
    StaticFiles(directory=str(settings.OUTPUT_DIR), html=True),
    name="reports",
)
