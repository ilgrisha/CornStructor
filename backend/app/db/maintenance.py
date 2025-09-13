# File: backend/app/db/maintenance.py
# Version: v0.2.0
"""
SQLite schema auto-heal for local/dev runs.

Now also:
- Adds `designs.ga_progress_json` if missing.
"""
from __future__ import annotations

from typing import Set, List
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
    stmt = """
    CREATE TABLE IF NOT EXISTS designs (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        created_at TEXT NOT NULL,
        updated_at TEXT NOT NULL,
        sequence TEXT NOT NULL,
        sequence_len INTEGER NOT NULL,
        params_json TEXT NULL,
        tree_json TEXT NULL,
        ga_progress_json TEXT NULL,
        name VARCHAR(120) NULL
    );
    """
    with engine.begin() as con:
        con.execute(text(stmt))


def _add_runs_design_id(engine: Engine) -> None:
    with engine.begin() as con:
        con.execute(text("ALTER TABLE runs ADD COLUMN design_id INTEGER NULL;"))


def _add_designs_ga_progress(engine: Engine) -> None:
    with engine.begin() as con:
        con.execute(text("ALTER TABLE designs ADD COLUMN ga_progress_json TEXT NULL;"))


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
    else:
        cols = _get_columns(engine, "designs")
        if "ga_progress_json" not in cols:
            _add_designs_ga_progress(engine)
            actions.append("added column designs.ga_progress_json")

    # 2) runs.design_id column
    if "runs" in tables:
        cols = _get_columns(engine, "runs")
        if "design_id" not in cols:
            _add_runs_design_id(engine)
            actions.append("added column runs.design_id")

    return actions
