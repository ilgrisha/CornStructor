# File: backend/app/api/v1/primers/models.py
# Version: v0.1.0
"""
SQLAlchemy models for primer design runs and results.
"""

from __future__ import annotations

import uuid
from sqlalchemy import Column, String, Integer, DateTime, ForeignKey, Text
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import relationship, declarative_base
from sqlalchemy.sql import func

# If you already have a shared Base, import it instead:
# from app.db.base import Base
Base = declarative_base()


class PrimerRun(Base):
    __tablename__ = "primer_runs"

    id = Column(String, primary_key=True, default=lambda: str(uuid.uuid4()))
    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)
    sequence_digest = Column(String, nullable=False, index=True)
    target_start = Column(Integer, nullable=False)
    target_end = Column(Integer, nullable=False)
    parameters_json = Column(JSONB, nullable=False)
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
    warnings_json = Column(JSONB, nullable=False)

    # Optional cached alignment blocks
    alignment_blocks_json = Column(JSONB, nullable=True)

    run = relationship("PrimerRun", back_populates="result")
