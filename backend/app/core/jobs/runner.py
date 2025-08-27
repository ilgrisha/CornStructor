# File: backend/app/core/jobs/runner.py
# Version: v0.4.1
"""
Async job runner: invokes the CLI pipeline and streams logs via an asyncio.Queue.

- POST /api/design/start -> creates a job and spawns the CLI
- GET  /api/design/{job_id}/logs -> SSE that drains the job's queue

v0.4.1
------
- Persist run lifecycle into the central SQL database (RUNNING -> COMPLETED/FAILED)
- Keep emitting RESULT links, including /reports/{job_id}/index.html when generated.
- Preserve writing FASTA to a file and invoking CLI with --fasta <file>.
"""
from __future__ import annotations

import asyncio
import shlex
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional

from backend.app.core.config import settings
from backend.app.core.visualization.run_index import write_run_index

# New: persistence
from backend.app.db.session import SessionLocal
from backend.app.services.run_store import record_run_start, record_run_completion

@dataclass
class Job:
    job_id: str
    outdir: Path
    queue: asyncio.Queue[str]
    task: Optional[asyncio.Task] = None

class JobManager:
    def __init__(self) -> None:
        self._jobs: Dict[str, Job] = {}

    def create(self) -> Job:
        job_id = uuid.uuid4().hex[:12]
        outdir = Path(settings.OUTPUT_DIR) / job_id
        outdir.mkdir(parents=True, exist_ok=True)
        job = Job(job_id=job_id, outdir=outdir, queue=asyncio.Queue())
        self._jobs[job_id] = job

        # Persist RUNNING row
        try:
            db = SessionLocal()
            record_run_start(db, job_id=job_id, sequence_len=None)
        finally:
            db.close()
        return job

    def get(self, job_id: str) -> Optional[Job]:
        return self._jobs.get(job_id)

    async def run_cli_job(self, job: Job, *, fasta_text: str, sequence_id: str = "seq1") -> None:
        """Spawn the CLI pipeline process, stream logs, and finalize artifacts."""
        # Write FASTA file
        fasta_path = job.outdir / f"{sequence_id}.fasta"
        fasta_clean = fasta_text.strip().upper()
        fasta_path.write_text(f">{sequence_id}\n{fasta_clean}\n", encoding="utf-8")

        # update sequence length if we can estimate
        seq_len = sum(1 for ch in fasta_clean if ch in "ACGTN")
        try:
            db = SessionLocal()
            record_run_start(db, job_id=job.job_id, sequence_len=seq_len)
        finally:
            db.close()

        cmd = [
            "python",
            "-m",
            "backend.app.cli.assemble_cli",
            "--fasta",
            str(fasta_path),
            "--levels",
            str(settings.LEVELS_PATH),
            "--globals",
            str(settings.GLOBALS_PATH),
            "--outdir",
            str(job.outdir),
            "--log-level",
            "INFO",
        ]

        await job.queue.put(f"RUN: {shlex.join(cmd)}")

        proc = await asyncio.create_subprocess_exec(
            *cmd, stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.PIPE
        )

        async def pump(stream, prefix: str):
            assert stream is not None
            while True:
                line = await stream.readline()
                if not line:
                    break
                await job.queue.put(f"{prefix}: {line.decode().rstrip()}")

        await asyncio.gather(pump(proc.stdout, "STDOUT"), pump(proc.stderr, "STDERR"))
        rc = await proc.wait()

        # Discover artifacts and write index
        for d in job.outdir.iterdir():
            if d.is_dir() and d.name.endswith("_cluster_reports"):
                await job.queue.put(f"RESULT: {settings.REPORTS_PUBLIC_BASE}/{job.job_id}/{d.name}/")

        index_path = write_run_index(job.outdir, job.job_id, reports_public_base=settings.REPORTS_PUBLIC_BASE)
        if index_path and index_path.exists():
            await job.queue.put(f"RESULT: {settings.REPORTS_PUBLIC_BASE}/{job.job_id}/index.html")

        # Persist completion
        try:
            db = SessionLocal()
            report_url = f"{settings.REPORTS_PUBLIC_BASE}/{job.job_id}/index.html" if index_path and index_path.exists() else None
            record_run_completion(db, job_id=job.job_id, report_url=report_url, exit_code=rc)
        finally:
            db.close()

        await job.queue.put(f"EXIT: {rc}")
        await job.queue.put("__EOF__")

    async def start(self, fasta_text: str) -> str:
        job = self.create()
        job.task = asyncio.create_task(self.run_cli_job(job, fasta_text=fasta_text))
        return job.job_id

job_manager = JobManager()
