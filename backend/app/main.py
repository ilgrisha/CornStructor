# File: backend/app/main.py
# Version: v0.2.3
"""
FastAPI app entry with optional SQLite auto-heal and robust reports serving.

- Mounts /reports to serve files from OUTPUT_DIR (StaticFiles)
- Adds a precise fallback route to serve files directly (FileResponse) if needed
- Startup logs confirm what directory is actually mounted
"""
from __future__ import annotations

import os
from pathlib import Path
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse
from starlette.staticfiles import StaticFiles

from backend.app.api.v1.api import api_router
from backend.app.db.session import engine
from backend.app.db.maintenance import ensure_schema_sqlite
from backend.app.core.config import settings

app = FastAPI(title="CornStructor API", openapi_url="/api/openapi.json", docs_url="/api/docs")

# CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # tighten in prod
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# API
app.include_router(api_router, prefix="/api")

# --- Reports static mount ---
reports_base = settings.REPORTS_PUBLIC_BASE.rstrip("/") or "/reports"
reports_dir = Path(settings.OUTPUT_DIR).resolve()

# Log what we’re mounting (helps diagnose mismatches)
print(f"[reports] mounting '{reports_base}' -> '{reports_dir}' (exists={reports_dir.exists()})")

# Serve static
app.mount(reports_base, StaticFiles(directory=str(reports_dir), html=False), name="reports")

# Explicit fallback: if StaticFiles misses, this route still serves any file
@app.get(f"{reports_base}" + "/{job_id}/{path:path}")
def _reports_fallback(job_id: str, path: str):
    candidate = reports_dir / job_id / path
    if candidate.is_file():
        return FileResponse(str(candidate))
    raise HTTPException(status_code=404, detail="Not Found")


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

# Quick diagnostics endpoint to verify what the backend sees
@app.get("/api/diag/reports")
def _diag_reports():
    try:
        jobs = []
        if reports_dir.exists():
            for p in sorted((reports_dir.glob("*/index.html"))):
                jobs.append(str(p.parent.name))
        return JSONResponse(
            {
                "reports_base": reports_base,
                "reports_dir": str(reports_dir),
                "exists": reports_dir.exists(),
                "job_ids_with_index": jobs,
            }
        )
    except Exception as ex:
        return JSONResponse({"error": repr(ex)}, status_code=500)
