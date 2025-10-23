# File: backend/alembic/versions/20251021_000001_add_primers_tables.py
# Version: v0.2.0
"""
Create tables: primer_runs, primer_results

Updated: set `down_revision` to the current head revision of your project.
If your latest migration id is different, update `down_revision` below accordingly.
"""

from __future__ import annotations

from alembic import op
import sqlalchemy as sa

# ✅ Set this to your current head revision BEFORE this migration.
# Example: down_revision = "a1b2c3d4e5f6"
down_revision = "SET_ME_TO_PREVIOUS_REVISION_ID"  # <-- UPDATE THIS VALUE

# revision identifiers, used by Alembic.
revision = "20251021_000001"
branch_labels = None
depends_on = None


def upgrade():
    op.create_table(
        "primer_runs",
        sa.Column("id", sa.String(), primary_key=True),
        sa.Column("created_at", sa.DateTime(timezone=True), server_default=sa.func.now(), nullable=False),
        sa.Column("sequence_digest", sa.String(), nullable=False, index=True),
        sa.Column("target_start", sa.Integer(), nullable=False),
        sa.Column("target_end", sa.Integer(), nullable=False),
        sa.Column("parameters_json", sa.dialects.postgresql.JSONB(), nullable=False),
        sa.Column("status", sa.String(), nullable=False, server_default="completed"),
    )
    op.create_table(
        "primer_results",
        sa.Column("id", sa.String(), primary_key=True),
        sa.Column("run_id", sa.String(), sa.ForeignKey("primer_runs.id", ondelete="CASCADE"), nullable=False),
        sa.Column("forward_seq", sa.Text(), nullable=False),
        sa.Column("forward_tm", sa.String(), nullable=False),
        sa.Column("forward_gc", sa.String(), nullable=False),
        sa.Column("reverse_seq", sa.Text(), nullable=False),
        sa.Column("reverse_tm", sa.String(), nullable=False),
        sa.Column("reverse_gc", sa.String(), nullable=False),
        sa.Column("pair_score", sa.String(), nullable=False),
        sa.Column("warnings_json", sa.dialects.postgresql.JSONB(), nullable=False),
        sa.Column("alignment_blocks_json", sa.dialects.postgresql.JSONB(), nullable=True),
    )


def downgrade():
    op.drop_table("primer_results")
    op.drop_table("primer_runs")
