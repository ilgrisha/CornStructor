# File: backend/app/core/config.py
# Version: v0.3.0
"""
Centralized application settings using Pydantic Settings.

Controls:
- App metadata and API prefix
- CORS origins
- Paths to levels/globals JSON
- Output dir for generated reports (shared volume)
- Public base path for reports (served by Nginx in prod; by FastAPI StaticFiles in dev)
- Database URL (SQLAlchemy)
"""
from __future__ import annotations

from pathlib import Path
from pydantic_settings import BaseSettings, SettingsConfigDict

class Settings(BaseSettings):
    # --- App ---
    API_PREFIX: str = "/api"
    APP_NAME: str = "CornStructor"
    APP_VERSION: str = "0.3.0"

    # --- CORS ---
    CORS_ORIGINS: str = "*"  # comma-separated or '*' for all

    # --- Data / pipeline config ---
    LEVELS_PATH: Path = Path("backend/app/config/levels.json")
    GLOBALS_PATH: Path = Path("backend/app/config/globals.json")

    # --- Output / reports ---
    OUTPUT_DIR: Path = Path("backend/data/out")
    REPORTS_PUBLIC_BASE: str = "/reports"  # served by backend in dev, Nginx in prod

    # --- DB ---
    DB_URL: str = "sqlite:///backend/app/data/cornstructor.db"

    # --- Server ---
    HOST: str = "0.0.0.0"
    PORT: int = 8000
    LOG_LEVEL: str = "info"

    model_config = SettingsConfigDict(env_file=".env", env_file_encoding="utf-8")

    @property
    def cors_origins_list(self) -> list[str]:
        raw = self.CORS_ORIGINS.strip()
        if raw == "*":
            return ["*"]
        return [o.strip() for o in raw.split(",") if o.strip()]

settings = Settings()
