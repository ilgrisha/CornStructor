# File: backend/app/core/pipeline/py_runner.py
# Version: v0.1.2
"""
In-process CornStructor pipeline runner.

v0.1.2:
- After exporting GA progress, save its JSON payload into the Design (DB) so that
  live reports can be generated on-the-fly without reading from the reports folder.
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable, Optional
import json
import logging

from backend.app.core.assembly.hierarchical_assembler import HierarchicalAssembler
from backend.app.config.config_level import load_levels_config, LevelConfig
from backend.app.config.config_global import load_global_config, GlobalConfig

from backend.app.core.export.json_exporter import export_tree_to_json
from backend.app.core.export.analysis_exporter import analyze_tree_to_json
from backend.app.core.export.ga_progress_exporter import export_ga_progress_json
from backend.app.core.export.fasta_exporter import export_fragments_to_fasta

from backend.app.core.visualization.analysis_html_report import export_analysis_html
from backend.app.core.visualization.ga_progress_html_report import export_ga_progress_html
from backend.app.core.visualization.cluster_html_report import export_all_levels
from backend.app.core.visualization.tree_html_exporter import export_tree_to_html

from backend.app.core.visualization.run_index_simple import write_run_index_simple

# NEW: persist GA progress to DB
from backend.app.db.session import SessionLocal
from backend.app.services.design_store import save_ga_progress


class _SSELogHandler(logging.Handler):
    def __init__(self, emit_fn: Callable[[str], None]) -> None:
        super().__init__()
        self._emit = emit_fn

    def emit(self, record: logging.LogRecord) -> None:
        try:
            msg = self.format(record)
        except Exception:
            msg = record.getMessage()
        self._emit(msg)


def execute_pipeline(
    *,
    outdir: Path,
    fasta_text: str,
    job_id: str,
    levels_path: Path,
    globals_path: Path,
    reports_public_base: str = "/reports",
    emit_log: Optional[Callable[[str], None]] = None,
) -> int:
    logger = logging.getLogger(f"cornstructor.job.{job_id}")
    logger.setLevel(logging.INFO)
    logger.handlers = []
    if emit_log is not None:
        h = _SSELogHandler(emit_log)
        h.setFormatter(logging.Formatter("%(message)s"))
        logger.addHandler(h)

    outdir.mkdir(parents=True, exist_ok=True)
    fasta_path = outdir / "input.fasta"

    lines = [ln.strip() for ln in fasta_text.splitlines() if ln.strip()]
    if lines and lines[0].startswith(">"):
        seq = "".join(lines[1:])
    else:
        seq = "".join(lines)
    seq = seq.upper()
    fasta_path.write_text(f">seq\n{seq}\n", encoding="utf-8")
    emit_log and emit_log(f"Saved FASTA to {fasta_path}")

    levels_cfg: dict[int, LevelConfig] = load_levels_config(levels_path)
    global_cfg: GlobalConfig = load_global_config(str(globals_path))
    emit_log and emit_log("Loaded levels/global configs")

    assembler = HierarchicalAssembler(levels_cfg, global_cfg, logger=logger)
    root = assembler.build(full_seq=seq, root_id=job_id)
    emit_log and emit_log("Built construction tree")

    # Exports
    export_tree_to_json(root, outdir / "tree.json", root_id=job_id)
    emit_log and emit_log("Exported tree.json")

    clusters_dir = outdir / "clusters"
    export_all_levels(root, clusters_dir, levels_cfg=levels_cfg, global_cfg=global_cfg)
    emit_log and emit_log("Exported per-level cluster HTML")

    ga_json = export_ga_progress_json(root, outdir / "ga_progress.json", root_id=job_id)
    export_ga_progress_html(ga_json, outdir / "ga_progress.html")
    emit_log and emit_log("Exported GA progress reports")

    # Persist GA progress JSON to DB for live rendering
    try:
        db = SessionLocal()
        try:
            save_ga_progress(db, job_id=job_id, ga_progress_json=json.dumps(ga_json))
            emit_log and emit_log("Saved GA progress to database")
        finally:
            db.close()
    except Exception as ex:
        emit_log and emit_log(f"WARN: failed to save GA progress to DB: {ex}")

    analysis_json = analyze_tree_to_json(root, outdir / "analysis.json", root_id=job_id)
    export_analysis_html(analysis_json, outdir / "analysis.html")
    emit_log and emit_log("Exported analysis reports")

    export_tree_to_html(root, outdir / "tree.html", root_id=job_id)
    export_fragments_to_fasta(root, outdir / "fragments.fasta", root_id=job_id)
    emit_log and emit_log("Exported tree HTML and fragments FASTA")

    # Robust index
    idx = write_run_index_simple(outdir, job_id, reports_public_base=reports_public_base)
    emit_log and emit_log(f"Generated index.html → {idx}")

    return 0
