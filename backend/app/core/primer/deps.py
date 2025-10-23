# File: backend/app/api/v1/primers/deps.py
# Version: v0.1.0
"""
Dependency providers for primer endpoints.
"""

from __future__ import annotations

from fastapi import Depends
from sqlalchemy.orm import Session

# Reuse existing db session factory in your project:
# from app.db.session import get_db
# For scaffold purposes, we declare a placeholder. Replace with your real get_db.
def get_db():
    """Yield a DB session (placeholder). Wire to your real session factory."""
    raise NotImplementedError("Wire get_db() to your existing session provider.")


def db_session(db: Session = Depends(get_db)) -> Session:
    """Return a SQLAlchemy session."""
    return db
