# File: backend/app/api/v1/primers/deps.py
# Version: v0.2.0
"""
Dependency providers for primer endpoints.

This version wires `get_db()` to CornStructor's existing DB session provider.
"""

from __future__ import annotations

from fastapi import Depends
from sqlalchemy.orm import Session

# ✅ Use the project's canonical DB session dependency
from backend.app.db.session import get_db  # type: ignore


def db_session(db: Session = Depends(get_db)) -> Session:
    """Return an active SQLAlchemy session."""
    return db
