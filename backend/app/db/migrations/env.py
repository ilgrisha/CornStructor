# File: backend/app/db/migrations/env.py
# Version: v0.2.0
"""
Alembic environment for CornStructor.

- Ensures repo root is on sys.path so `import backend...` works regardless of CWD.
- Uses DATABASE_URL env var when present, otherwise falls back to alembic.ini.
- Registers ALL SQLAlchemy models on Base.metadata so `--autogenerate` can detect schema changes,
  including the new primers tables (primer_runs, primer_results).
"""
from __future__ import annotations

import os
import sys
from logging.config import fileConfig

from sqlalchemy import engine_from_config, pool
from alembic import context

# --------------------------------------------------------------------------------------
# Make the repository importable no matter where Alembic is invoked from
# This file lives at: backend/app/db/migrations/env.py
# repo_root = <this_file>/../../../../
# --------------------------------------------------------------------------------------
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(THIS_DIR, "../../../../"))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

# --------------------------------------------------------------------------------------
# Import the project's SQLAlchemy Base and models
# IMPORTANT: We must import any module that defines models so that Base.metadata is populated.
# --------------------------------------------------------------------------------------
from backend.app.db.models import Base  # noqa: E402  (aggregates core app models)

# 🔗 NEW: Import primers models so autogenerate sees primer_runs and primer_results
# (side-effect import; no symbols referenced directly)
import backend.app.api.v1.primers.models  # noqa: F401,E402

# --------------------------------------------------------------------------------------
# Alembic configuration
# --------------------------------------------------------------------------------------
config = context.config

# Configure logging from alembic.ini
if config.config_file_name is not None:
    fileConfig(config.config_file_name)

# Target metadata for autogenerate
target_metadata = Base.metadata


def run_migrations_offline() -> None:
    """
    Run migrations in 'offline' mode'.

    In this mode we configure the context with just a URL
    and not an Engine, though an Engine is acceptable here as well.
    By skipping the Engine creation we don't even need a DBAPI to be available.
    """
    url = os.getenv("DATABASE_URL") or config.get_main_option("sqlalchemy.url")
    if not url:
        raise RuntimeError(
            "DATABASE_URL is not set and sqlalchemy.url is empty in alembic.ini. "
            "Set a valid connection string (e.g., sqlite:///./cornstructor.db)."
        )

    context.configure(
        url=url,
        target_metadata=target_metadata,
        literal_binds=True,
        dialect_opts={"paramstyle": "named"},
        compare_type=True,
        compare_server_default=True,
        # include_schemas=True,  # enable if you manage multiple schemas
    )

    with context.begin_transaction():
        context.run_migrations()


def run_migrations_online() -> None:
    """
    Run migrations in 'online' mode'.

    In this scenario we need to create an Engine and associate a connection with the context.
    """
    # Prefer DATABASE_URL env var (keeps parity with app runtime), else alembic.ini
    db_url = os.getenv("DATABASE_URL")
    if db_url:
        config.set_main_option("sqlalchemy.url", db_url)

    # engine_from_config reads keys prefixed with "sqlalchemy." in alembic.ini
    connectable = engine_from_config(
        config.get_section(config.config_ini_section),
        prefix="sqlalchemy.",
        poolclass=pool.NullPool,
    )

    with connectable.connect() as connection:
        context.configure(
            connection=connection,
            target_metadata=target_metadata,
            compare_type=True,
            compare_server_default=True,
            # render_as_batch=True,  # enable for older SQLite migrations requiring batch ops
        )

        with context.begin_transaction():
            context.run_migrations()


if context.is_offline_mode():
    run_migrations_offline()
else:
    run_migrations_online()
