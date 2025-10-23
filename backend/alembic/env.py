# File: backend/alembic/env.py
# Version: v0.1.0
"""
Alembic environment for CornStructor backend.

- Reads DATABASE_URL from environment (recommended).
- Falls back to alembic.ini sqlalchemy.url if env var is missing.
- Uses app.db.base.Base.metadata as target_metadata.
- Ensures models are imported so autogenerate sees them (incl. primers tables).
"""

from __future__ import annotations

import os
from logging.config import fileConfig

from sqlalchemy import engine_from_config, pool
from alembic import context

# this is the Alembic Config object, which provides access to the values within the .ini file in use.
config = context.config

# Interpret the config file for Python logging.
if config.config_file_name is not None:
    fileConfig(config.config_file_name)

# --- Database URL wiring ---
DB_URL = os.getenv("DATABASE_URL") or config.get_main_option("sqlalchemy.url") or ""
if not DB_URL:
    raise RuntimeError(
        "DATABASE_URL is not set and sqlalchemy.url is empty in alembic.ini. "
        "Set DATABASE_URL (e.g., postgresql+psycopg2://user:pass@host/db)."
    )
config.set_main_option("sqlalchemy.url", DB_URL)

# --- Target metadata ---
from app.db.base import Base  # our declarative Base

# Import models so they are registered on Base.metadata for autogenerate
# Keep imports local to avoid circular dependencies at runtime
import app.api.v1.primers.models  # noqa: F401

target_metadata = Base.metadata


def run_migrations_offline():
    """Run migrations in 'offline' mode."""
    url = config.get_main_option("sqlalchemy.url")
    context.configure(
        url=url,
        target_metadata=target_metadata,
        literal_binds=True,
        dialect_opts={"paramstyle": "named"},
        include_object=None,
        compare_type=True,
        compare_server_default=True,
    )

    with context.begin_transaction():
        context.run_migrations()


def run_migrations_online():
    """Run migrations in 'online' mode."""
    connectable = engine_from_config(
        config.get_section(config.config_ini_section) or {},
        prefix="sqlalchemy.",
        poolclass=pool.NullPool,
    )

    with connectable.connect() as connection:
        context.configure(
            connection=connection,
            target_metadata=target_metadata,
            include_object=None,
            compare_type=True,
            compare_server_default=True,
        )

        with context.begin_transaction():
            context.run_migrations()


if context.is_offline_mode():
    run_migrations_offline()
else:
    run_migrations_online()
