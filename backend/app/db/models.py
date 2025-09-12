# File: backend/app/db/models.py
# Version: v0.1.0
"""
ORM models for CornStructor.

Currently includes:
- Run: a single pipeline execution (aka "job") with its report URL and status.

We model `job_id` as a unique external identifier (short hex) used by the API
and filesystem folder names under `/reports/{job_id}`.
"""
from __future__ import annotations

from datetime import datetime
from sqlalchemy.orm import declarative_base, Mapped, mapped_column
from sqlalchemy import String, Integer, DateTime, Enum, Text, UniqueConstraint

Base = declarative_base()

from enum import Enum as PyEnum

class RunStatus(str, PyEnum):
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"

class Run(Base):
    """A single CornStructor run/pipeline execution."""
    __tablename__ = "runs"
    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    job_id: Mapped[str] = mapped_column(String(32), nullable=False, unique=True, index=True)
    created_at: Mapped[datetime] = mapped_column(DateTime, nullable=False, default=datetime.utcnow)
    updated_at: Mapped[datetime] = mapped_column(DateTime, nullable=False, default=datetime.utcnow)
    status: Mapped[RunStatus] = mapped_column(Enum(RunStatus), nullable=False, default=RunStatus.RUNNING)
    report_url: Mapped[str | None] = mapped_column(String(512), nullable=True)
    exit_code: Mapped[int | None] = mapped_column(Integer, nullable=True)
    sequence_len: Mapped[int | None] = mapped_column(Integer, nullable=True)
    params_json: Mapped[str | None] = mapped_column(Text, nullable=True)
    note: Mapped[str | None] = mapped_column(String(255), nullable=True)

    __table_args__ = (
        UniqueConstraint("job_id", name="uq_runs_job_id"),
    )

    def touch(self) -> None:
        """Update `updated_at` timestamp."""
        from datetime import datetime as _dt
        self.updated_at = _dt.utcnow()
