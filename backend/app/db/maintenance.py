# File: backend/app/db/maintenance.py
# Version: v0.2.0
"""
SQLite schema maintenance helpers (dev-only, non-destructive).

- ensure_schema_sqlite(engine): creates only tables that are missing.
- Critically, this module imports `backend.app.db.models` (not just Base),
  so ALL ORM models (Design, Run, etc.) are registered in Base.metadata.

Usage:
  Set env var SCHEMA_AUTOHEAL=true and keep your DB backend as sqlite.
  On app startup, ensure_schema_sqlite(engine) will run and create any
  missing tables, logging each action for visibility.

Notes:
  * Safe to run multiple times; it never drops or alters existing tables.
  * Use Alembic migrations for staging/production changes.
"""
from __future__ import annotations

from typing import List

from sqlalchemy import inspect
from sqlalchemy.engine import Engine

# IMPORTANT: import the MODELS MODULE for side effects so that all model classes
# are defined and attached to THIS module's Base before we inspect metadata.
import backend.app.db.models as models  # noqa: F401


def ensure_schema_sqlite(engine: Engine) -> List[str]:
    """
    Create any missing tables declared on models.Base.metadata.

    Returns a list of human-readable action strings (e.g., "created table runs").
    """
    Base = models.Base  # same Declarative Base that defines Design, Run, etc.

    inspector = inspect(engine)
    existing = set(inspector.get_table_names())
    defined = set(Base.metadata.tables.keys())

    missing = sorted(defined - existing)
    actions: List[str] = []

    for name in missing:
        table = Base.metadata.tables[name]
        # checkfirst guards against races / repeated calls
        table.create(bind=engine, checkfirst=True)
        actions.append(f"created table {name}")

    if not actions:
        actions.append("all tables present")

    return actions
