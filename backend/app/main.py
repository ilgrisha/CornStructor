# File: backend/app/main.py
# Version: v0.5.0
"""
FastAPI application entrypoint (aggregator style).

Mounts:
- /api/* via a single aggregated v1 router (see backend/app/api/v1/api.py)
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
from backend.app.api.v1.api import api_router as api_v1_router

# DB init
from backend.app.db.session import engine
from backend.app.db.models import Base

# Create tables at startup if they don't exist (simple bootstrap; Alembic recommended later)
Base.metadata.create_all(bind=engine)

# App metadata is sourced from settings to keep it consistent across services
app = FastAPI(title=settings.APP_NAME, version=settings.APP_VERSION)

# CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.cors_origins_list,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Single include line: mount all v1 routes under the configured API prefix
# Example: if settings.API_PREFIX == "/api", routes become /api/health, /api/design, etc.
app.include_router(api_v1_router, prefix=settings.API_PREFIX)

# Dev/Parity: also serve reports directly from the backend (Nginx handles this in prod)
app.mount(
    "/reports",
    StaticFiles(directory=str(settings.OUTPUT_DIR), html=True),
    name="reports",
)
