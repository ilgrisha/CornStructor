# File: backend/app/db/migrations/versions/20250912_0002_add_ga_progress_to_designs.py
# Version: v0.1.0
"""
Add `ga_progress_json` column to designs (idempotent).

Revision ID: 0002_add_ga_progress_to_designs
Revises: 0001_add_designs_and_link_runs
Create Date: 2025-09-12
"""
from __future__ import annotations

from alembic import op
import sqlalchemy as sa

# revision identifiers
revision = "0002_add_ga_progress_to_designs"
down_revision = "0001_add_designs_and_link_runs"
branch_labels = None
depends_on = None


def upgrade() -> None:
    bind = op.get_bind()
    insp = sa.inspect(bind)
    if "designs" in insp.get_table_names():
        cols = {c["name"] for c in insp.get_columns("designs")}
        if "ga_progress_json" not in cols:
            with op.batch_alter_table("designs") as batch:
                batch.add_column(sa.Column("ga_progress_json", sa.Text(), nullable=True))


def downgrade() -> None:
    bind = op.get_bind()
    insp = sa.inspect(bind)
    if "designs" in insp.get_table_names():
        cols = {c["name"] for c in insp.get_columns("designs")}
        if "ga_progress_json" in cols:
            with op.batch_alter_table("designs") as batch:
                batch.drop_column("ga_progress_json")
