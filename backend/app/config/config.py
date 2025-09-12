# File: backend/app/core/config.py
# Version: v0.2.3
from __future__ import annotations

from pathlib import Path
from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    APP_NAME: str = "CornStructor API"
    APP_VERSION: str = "0.2.3"
    API_PREFIX: str = "/api"

    CORS_ORIGINS: str = "*"

    LEVELS_PATH: Path = Path("backend/app/config/levels.json")
    GLOBALS_PATH: Path = Path("backend/app/config/globals.json")

    OUTPUT_DIR: Path = Path("/data/out")
    REPORTS_PUBLIC_BASE: str = "/reports"

    HOST: str = "0.0.0.0"
    PORT: int = 8000
    LOG_LEVEL: str = "info"

    # Critical bits:
    # - extra="allow": unknown env vars (e.g., BACKEND_PORT) won't crash
    # - env_file=None: do NOT auto-load any .env inside the container
    model_config = SettingsConfigDict(
        extra="allow",
        env_file=None,
        env_file_encoding="utf-8",
    )

    @property
    def cors_origins_list(self) -> list[str]:
        raw = self.CORS_ORIGINS.strip()
        if raw == "*":
            return ["*"]
        return [o.strip() for o in raw.split(",") if o.strip()]


settings = Settings()
