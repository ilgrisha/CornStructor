# File: backend/app/core/config.py
# Version: v0.2.0
"""
Centralized application settings using Pydantic Settings.

Controls:
- App metadata and API prefix
- CORS origins
- Paths to levels/globals JSON
- Output dir for generated reports (shared volume)
- Public base path for reports (served by Nginx in prod; by FastAPI StaticFiles in dev)
"""
from __future__ import annotations

from pathlib import Path
from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    APP_NAME: str = "CornStructor API"
    APP_VERSION: str = "0.2.0"
    API_PREFIX: str = "/api"

    # CORS (comma-separated origins or "*")
    CORS_ORIGINS: str = "*"

    # Input config files (existing in repo)
    LEVELS_PATH: Path = Path("backend/app/config/levels.json")
    GLOBALS_PATH: Path = Path("backend/app/config/globals.json")

    # Output dir (shared volume)
    OUTPUT_DIR: Path = Path("/data/out")

    # Public URL base for reports
    REPORTS_PUBLIC_BASE: str = "/reports"

    # Uvicorn bind
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
