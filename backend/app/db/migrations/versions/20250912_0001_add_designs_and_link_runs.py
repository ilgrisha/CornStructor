# File: backend/app/db/migrations/versions/20250912_0001_add_designs_and_link_runs.py
# Version: v0.1.1
"""
Add `designs` table and link `runs.design_id` (one-to-one).

Idempotent for SQLite/local dev:
- Skips creating `designs` if it already exists.
- Skips adding `runs.design_id` / FK if the column already exists.

Revision ID: 0001_add_designs_and_link_runs
Revises: None
Create Date: 2025-09-12

Run:
  alembic upgrade head
"""
from __future__ import annotations

from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = "0001_add_designs_and_link_runs"
down_revision = None
branch_labels = None
depends_on = None


def upgrade() -> None:
    bind = op.get_bind()
    insp = sa.inspect(bind)

    # 1) Create designs table if missing
    tables = set(insp.get_table_names())
    if "designs" not in tables:
        op.create_table(
            "designs",
            sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
            sa.Column("created_at", sa.DateTime(), nullable=False),
            sa.Column("updated_at", sa.DateTime(), nullable=False),
            sa.Column("sequence", sa.Text(), nullable=False),
            sa.Column("sequence_len", sa.Integer(), nullable=False),
            sa.Column("params_json", sa.Text(), nullable=True),
            sa.Column("tree_json", sa.Text(), nullable=True),
            sa.Column("name", sa.String(length=120), nullable=True),
        )

    # 2) Add runs.design_id + FK if missing
    if "runs" in tables:
        run_cols = {c["name"] for c in insp.get_columns("runs")}
        if "design_id" not in run_cols:
            with op.batch_alter_table("runs") as batch:
                batch.add_column(sa.Column("design_id", sa.Integer(), nullable=True))
                # Create FK only if designs exists (it should by now)
                batch.create_foreign_key(
                    "fk_runs_design_id_designs",
                    "designs",
                    ["design_id"],
                    ["id"],
                    ondelete="SET NULL",
                )
        else:
            # Column exists; ensure FK exists (best-effort).
            # SQLite often reports no FKs; skip silently to keep idempotent behavior.
            pass


def downgrade() -> None:
    # Best-effort downgrade (drops FK & column if present; drops table if present)
    bind = op.get_bind()
    insp = sa.inspect(bind)
    tables = set(insp.get_table_names())

    if "runs" in tables:
        with op.batch_alter_table("runs") as batch:
            # Constraint name must match upgrade()
            try:
                batch.drop_constraint("fk_runs_design_id_designs", type_="foreignkey")
            except Exception:
                pass
            run_cols = {c["name"] for c in insp.get_columns("runs")}
            if "design_id" in run_cols:
                batch.drop_column("design_id")

    if "designs" in tables:
        op.drop_table("designs")
