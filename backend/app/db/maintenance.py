# File: backend/app/db/maintenance.py
# Version: v0.1.0
"""
SQLite schema auto-heal for local/dev runs.

This module inspects the current SQLite schema and:
- Creates `designs` table if missing.
- Adds `runs.design_id` column if missing.
- Adds a lightweight FK (SQLite won't enforce; Alembic migration adds proper FK elsewhere).

This is a pragmatic unblocker for environments without Alembic yet.
It is idempotent and safe to call on every startup.
"""
from __future__ import annotations

from typing import Set, Tuple, List, Any
from sqlalchemy import text
from sqlalchemy.engine import Engine


def _get_tables(engine: Engine) -> Set[str]:
    with engine.connect() as con:
        rows = con.execute(text("SELECT name FROM sqlite_master WHERE type='table';")).fetchall()
    return {r[0] for r in rows}


def _get_columns(engine: Engine, table: str) -> Set[str]:
    with engine.connect() as con:
        rows = con.execute(text(f"PRAGMA table_info('{table}')")).fetchall()
    return {r[1] for r in rows}  # r[1] is column name


def _create_designs_table(engine: Engine) -> None:
    # Mirrors backend/app/db/models.py
    stmt = """
    CREATE TABLE IF NOT EXISTS designs (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        created_at TEXT NOT NULL,
        updated_at TEXT NOT NULL,
        sequence TEXT NOT NULL,
        sequence_len INTEGER NOT NULL,
        params_json TEXT NULL,
        tree_json TEXT NULL,
        name VARCHAR(120) NULL
    );
    """
    with engine.begin() as con:
        con.execute(text(stmt))


def _add_runs_design_id(engine: Engine) -> None:
    # SQLite supports simple ALTER TABLE ADD COLUMN, no FK enforcement here.
    with engine.begin() as con:
        con.execute(text("ALTER TABLE runs ADD COLUMN design_id INTEGER NULL;"))


def ensure_schema_sqlite(engine: Engine) -> List[str]:
    """
    Ensures the minimal schema needed for design<->run linkage exists.

    Returns a list of actions performed for logging/diagnostics.
    """
    actions: List[str] = []
    tables = _get_tables(engine)

    # 1) designs table
    if "designs" not in tables:
        _create_designs_table(engine)
        actions.append("created table designs")

    # 2) runs.design_id column
    if "runs" in tables:
        cols = _get_columns(engine, "runs")
        if "design_id" not in cols:
            _add_runs_design_id(engine)
            actions.append("added column runs.design_id")
    else:
        # Nothing to do; runs table will be created by metadata on first run
        pass

    return actions
