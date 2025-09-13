# File: backend/app/db/models.py
# Version: v0.3.0
"""
ORM models for CornStructor.

Tables:
- Run: a single pipeline execution (aka "job") with its report URL and status.
- Design: persisted construction-tree design (sequence, params, serialized tree,
          and optional GA progress) linked 1:1 to a Run via `runs.design_id`.
"""
from __future__ import annotations

from datetime import datetime
from enum import Enum as PyEnum
from typing import Optional

from sqlalchemy import (
    String, Integer, DateTime, Enum, Text, UniqueConstraint, ForeignKey
)
from sqlalchemy.orm import declarative_base, Mapped, mapped_column, relationship

Base = declarative_base()


class RunStatus(str, PyEnum):
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"


class Design(Base):
    """Persisted construction-tree design.

    Stores the input DNA `sequence`, creation parameters as JSON
    (both global and per-level under a single JSON blob), the serialized
    construction tree JSON (`tree_json`), and optional GA progress JSON
    (`ga_progress_json`) to support on-the-fly report generation.
    """
    __tablename__ = "designs"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    created_at: Mapped[datetime] = mapped_column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at: Mapped[datetime] = mapped_column(DateTime, default=datetime.utcnow, nullable=False)

    # Core payload
    sequence: Mapped[str] = mapped_column(Text, nullable=False)
    sequence_len: Mapped[int] = mapped_column(Integer, nullable=False)
    params_json: Mapped[Optional[str]] = mapped_column(Text, nullable=True)   # JSON string (globals + levels)
    tree_json: Mapped[Optional[str]] = mapped_column(Text, nullable=True)     # JSON string (export_tree_to_json)
    ga_progress_json: Mapped[Optional[str]] = mapped_column(Text, nullable=True)  # JSON string (export_ga_progress_json)

    # Optional friendly label
    name: Mapped[Optional[str]] = mapped_column(String(120), nullable=True)

    # Reverse one-to-one to Run (Run.design_id -> designs.id)
    run: Mapped["Run"] = relationship(back_populates="design", uselist=False)


class Run(Base):
    """A single CornStructor run/pipeline execution."""
    __tablename__ = "runs"

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    job_id: Mapped[str] = mapped_column(String(32), nullable=False, unique=True, index=True)
    created_at: Mapped[datetime] = mapped_column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at: Mapped[datetime] = mapped_column(DateTime, default=datetime.utcnow, nullable=False)

    status: Mapped[RunStatus] = mapped_column(Enum(RunStatus), default=RunStatus.RUNNING, nullable=False)
    report_url: Mapped[Optional[str]] = mapped_column(String(512), nullable=True)
    exit_code: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)

    # Optional stats & metadata captured at start or completion
    sequence_len: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    params_json: Mapped[Optional[str]] = mapped_column(Text, nullable=True)
    note: Mapped[Optional[str]] = mapped_column(String(255), nullable=True)

    # One-to-one link to a Design
    design_id: Mapped[Optional[int]] = mapped_column(ForeignKey("designs.id", ondelete="SET NULL"), nullable=True)
    design: Mapped[Optional[Design]] = relationship(back_populates="run", lazy="joined")

    __table_args__ = (UniqueConstraint("job_id", name="uq_runs_job_id"),)

    def touch(self) -> None:
        """Update `updated_at` timestamp."""
        self.updated_at = datetime.utcnow()
