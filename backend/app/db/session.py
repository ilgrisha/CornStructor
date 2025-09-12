# File: backend/app/db/session.py
# Version: v0.1.0
"""
SQLAlchemy engine and session factory.

- Uses SQLite by default (file path from SETTINGS.DB_URL or `sqlite:///./data/cornstructor.db`)
- Provides `get_db()` FastAPI dependency to manage session lifecycle.
- Creates the `data/` directory if using SQLite file URLs.

This is intentionally synchronous for simplicity. For our current write volume
(~a few rows per run), sync SQLAlchemy is perfectly adequate and simpler to test.
"""
from __future__ import annotations

from pathlib import Path
from typing import Generator

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, Session

from backend.app.core.config import settings

# Resolve DB URL (default to a local SQLite file inside backend/app/data)
DB_URL = settings.DB_URL
if DB_URL.startswith("sqlite:///"):
    db_path = DB_URL.replace("sqlite:///", "", 1)
    db_dir = Path(db_path).resolve().parent
    db_dir.mkdir(parents=True, exist_ok=True)

engine = create_engine(DB_URL, future=True)
SessionLocal = sessionmaker(bind=engine, autoflush=False, autocommit=False, future=True)

def get_db() -> Generator[Session, None, None]:
    """Yield a database session and guarantee closing it after use."""
    db: Session = SessionLocal()
    try:
        yield db
    finally:
        db.close()
