# File: backend/app/api/v1/primers/models.py
# Version: v0.3.0
"""
SQLAlchemy models for primer design runs and results.

This version:
- Uses the project's shared Base (backend.app.db.base.Base) so Alembic --autogenerate sees these tables.
- SWITCHES from PostgreSQL JSONB to cross-dialect SQLAlchemy JSON for SQLite compatibility in dev.
"""

from __future__ import annotations

import uuid
from sqlalchemy import Column, String, Integer, DateTime, ForeignKey, Text, JSON
from sqlalchemy.orm import relationship
from sqlalchemy.sql import func

from backend.app.db.base import Base  # ✅ shared Base


class PrimerRun(Base):
    __tablename__ = "primer_runs"

    id = Column(String, primary_key=True, default=lambda: str(uuid.uuid4()))
    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)
    sequence_digest = Column(String, nullable=False, index=True)
    target_start = Column(Integer, nullable=False)
    target_end = Column(Integer, nullable=False)
    parameters_json = Column(JSON, nullable=False)  # ✅ cross-dialect
    status = Column(String, nullable=False, default="completed")

    result = relationship("PrimerResult", back_populates="run", uselist=False, cascade="all, delete-orphan")


class PrimerResult(Base):
    __tablename__ = "primer_results"

    id = Column(String, primary_key=True, default=lambda: str(uuid.uuid4()))
    run_id = Column(String, ForeignKey("primer_runs.id", ondelete="CASCADE"), nullable=False)

    forward_seq = Column(Text, nullable=False)
    forward_tm = Column(String, nullable=False)
    forward_gc = Column(String, nullable=False)

    reverse_seq = Column(Text, nullable=False)
    reverse_tm = Column(String, nullable=False)
    reverse_gc = Column(String, nullable=False)

    pair_score = Column(String, nullable=False)
    warnings_json = Column(JSON, nullable=False)            # ✅ cross-dialect
    alignment_blocks_json = Column(JSON, nullable=True)     # ✅ cross-dialect

    run = relationship("PrimerRun", back_populates="result")
