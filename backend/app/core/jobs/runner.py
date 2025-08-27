# File: backend/app/core/jobs/runner.py
# Version: v0.2.0
"""
Async job runner: invokes the CLI pipeline and streams logs via an asyncio.Queue.

- POST /api/design/start -> creates a job and spawns the CLI
- GET  /api/design/{job_id}/logs -> SSE that drains the job's queue
"""
from __future__ import annotations

import asyncio
import shlex
import uuid
from pathlib import Path
from typing import Dict, Optional

from backend.app.core.config import settings


class Job:
    def __init__(self, job_id: str, outdir: Path) -> None:
        self.job_id = job_id
        self.outdir = outdir
        self.queue: asyncio.Queue[str] = asyncio.Queue()
        self.task: Optional[asyncio.Task] = None


class JobManager:
    def __init__(self) -> None:
        self._jobs: Dict[str, Job] = {}

    def create(self) -> Job:
        job_id = uuid.uuid4().hex[:12]
        outdir = settings.OUTPUT_DIR / job_id
        outdir.mkdir(parents=True, exist_ok=True)
        job = Job(job_id, outdir)
        self._jobs[job_id] = job
        return job

    def get(self, job_id: str) -> Optional[Job]:
        return self._jobs.get(job_id)

    async def run_cli_job(self, job: Job, fasta_text: str, sequence_id: str = "seq1") -> None:
        fasta_path = job.outdir / f"{sequence_id}.fasta"
        fasta_path.write_text(f">{sequence_id}\n{fasta_text.strip().upper()}\n", encoding="utf-8")

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
            "--log",
            "INFO",
        ]

        await job.queue.put(f"RUN: {shlex.join(cmd)}")

        proc = await asyncio.create_subprocess_exec(
            *cmd, stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.PIPE
        )

        async def _pump(stream, prefix: str):
            assert stream is not None
            while True:
                line = await stream.readline()
                if not line:
                    break
                await job.queue.put(f"{prefix}{line.decode().rstrip()}")

        out_task = asyncio.create_task(_pump(proc.stdout, ""))
        err_task = asyncio.create_task(_pump(proc.stderr, "ERR: "))

        rc = await proc.wait()
        await asyncio.gather(out_task, err_task)

        # Expose artifacts
        for p in sorted(job.outdir.glob("*.html")):
            await job.queue.put(f"RESULT: {settings.REPORTS_PUBLIC_BASE}/{job.job_id}/{p.name}")
        for d in job.outdir.iterdir():
            if d.is_dir() and d.name.endswith("_cluster_reports"):
                await job.queue.put(f"RESULT: {settings.REPORTS_PUBLIC_BASE}/{job.job_id}/{d.name}/")

        await job.queue.put(f"EXIT: {rc}")
        await job.queue.put("__EOF__")

    async def start(self, fasta_text: str) -> str:
        job = self.create()
        job.task = asyncio.create_task(self.run_cli_job(job, fasta_text=fasta_text))
        return job.job_id


job_manager = JobManager()
