# File: backend/app/main.py
# Version: v0.3.1
"""
FastAPI app entry.

- Keeps all route assembly in backend/app/api/v1/api.py.
- Mounts /api/* via `api_router`.
- Mounts /reports/* at root via `public_router` (no direct import of reports_dynamic).
- Optional SQLite auto-heal is guarded by SCHEMA_AUTOHEAL.
"""
from __future__ import annotations

import os
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from backend.app.api.v1.api import api_router, public_router
from backend.app.db.session import engine
from backend.app.db.maintenance import ensure_schema_sqlite

# ✅ Primers router
from backend.app.api.v1.primers.router import router as primers_router


app = FastAPI(title="CornStructor API", openapi_url="/api/openapi.json", docs_url="/api/docs")

# CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # tighten in prod
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ✅ Mount primers endpoints
app.include_router(primers_router)


@app.get("/healthz")
def healthz():
    return {"status": "ok"}


# APIs under /api
app.include_router(api_router, prefix="/api")

# Public dynamic reports at root /reports/*
app.include_router(public_router)


def _env_true(name: str, default: str = "false") -> bool:
    val = os.getenv(name, default).strip().lower()
    return val in ("1", "true", "yes")


@app.on_event("startup")
def _startup_autoheal() -> None:
    # Only try auto-heal if explicitly enabled AND using SQLite
    if engine.url.get_backend_name() == "sqlite" and _env_true("SCHEMA_AUTOHEAL", "false"):
        actions = ensure_schema_sqlite(engine)
        if actions:
            print("[schema-autoheal]", ", ".join(actions))
