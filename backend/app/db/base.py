# File: backend/app/db/base.py
# Version: v0.1.0
"""
Declarative Base for CornStructor.

Other model modules should import Base from here:

    from app.db.base import Base

We do NOT import model modules here to avoid circular imports. Instead,
alembic/env.py will import model packages before autogenerate.
"""

from __future__ import annotations

from sqlalchemy.orm import DeclarativeBase


class Base(DeclarativeBase):
    pass
