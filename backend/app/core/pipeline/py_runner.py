# File: backend/app/core/pipeline/py_runner.py
# Version: v0.3.0
"""
In-process CornStructor pipeline runner.

v0.3.0
------
- Persist **config snapshots** (`globals`, `levels`) into `Design.params_json`.
- Persist **GA progress** JSON into `Design.ga_progress_json`.
- Attach serialized **tree.json** to the linked `Design` record.
- Mark the `Run` **COMPLETED/FAILED** and set `report_url` to the dynamic
  renderer at `/reports/{job_id}/index.html`.
- Emit stable log lines the UI expects:
    * "Saved config snapshots to database"
    * "Stored design tree in database"
    * "RESULT: /reports/<job_id>/index.html"
    * "EXIT: <code>"

This runner builds the construction tree with the given configuration files,
exports artifacts to a transient `outdir` (used for local inspection / dev),
and persists all data necessary to **re-generate reports on the fly** later
(from the database) after container restarts.
"""
from __future__ import annotations

import json
import logging
import traceback
from pathlib import Path
from typing import Callable, Optional

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

# DB interactions
from backend.app.db.session import SessionLocal
from backend.app.services.design_store import (
    attach_tree_to_design,
    mark_run_completed,
    save_config_snapshot,
    save_ga_progress,
)


class _SSELogHandler(logging.Handler):
    """Forward log records to a simple callable used by the log SSE stream."""
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
    """
    Run the CornStructor pipeline for a single job.

    Parameters
    ----------
    outdir : Path
        Transient output directory used while running (dev convenience).
    fasta_text : str
        Input sequence text (FASTA or raw). If FASTA, only the first record is used.
    job_id : str
        Stable job identifier (also used as root node id).
    levels_path : Path
        Path to the per-level configuration JSON used for this run.
    globals_path : Path
        Path to the global configuration JSON used for this run.
    reports_public_base : str
        Base URL where dynamic reports are served (default: "/reports").
    emit_log : Optional[Callable[[str], None]]
        Callback used to stream logs to the UI (SSE). If None, no streaming.

    Returns
    -------
    int
        0 on success; non-zero on failure. Also updates the Run row accordingly.
    """
    logger = logging.getLogger(f"cornstructor.job.{job_id}")
    logger.setLevel(logging.INFO)
    logger.handlers = []
    if emit_log is not None:
        h = _SSELogHandler(emit_log)
        h.setFormatter(logging.Formatter("%(message)s"))
        logger.addHandler(h)

    def _log(msg: str) -> None:
        if emit_log:
            emit_log(msg)
        else:
            logger.info(msg)

    rc = 1  # assume failure until the end
    tree_json_path: Optional[Path] = None

    try:
        # --- Prepare IO / input sequence ---
        outdir.mkdir(parents=True, exist_ok=True)
        fasta_path = outdir / "input.fasta"

        lines = [ln.strip() for ln in fasta_text.splitlines() if ln.strip()]
        if lines and lines[0].startswith(">"):
            seq = "".join(lines[1:])
        else:
            seq = "".join(lines)
        seq = seq.upper()
        fasta_path.write_text(f">seq\n{seq}\n", encoding="utf-8")
        _log(f"Saved FASTA to {fasta_path}")

        # --- Load configuration ---
        levels_cfg: dict[int, LevelConfig] = load_levels_config(levels_path)
        global_cfg: GlobalConfig = load_global_config(str(globals_path))
        _log("Loaded levels/global configs")

        # Persist exact config snapshots to DB so we can re-render on-the-fly later
        try:
            gl_text = Path(globals_path).read_text(encoding="utf-8")
            lv_text = Path(levels_path).read_text(encoding="utf-8")
            db = SessionLocal()
            try:
                save_config_snapshot(db, job_id=job_id, globals_json=gl_text, levels_json=lv_text)
                _log("Saved config snapshots to database")
            finally:
                db.close()
        except Exception as ex:
            _log(f"WARN: failed to save config snapshots: {ex}")

        # --- Build Tree ---
        assembler = HierarchicalAssembler(levels_cfg, global_cfg, logger=logger)
        root = assembler.build(full_seq=seq, root_id=job_id)
        _log("Built construction tree")

        # --- Export core artifacts to the transient outdir (for dev/local) ---
        tree_json_path = outdir / "tree.json"
        export_tree_to_json(root, tree_json_path, root_id=job_id)
        _log("Exported tree.json")

        clusters_dir = outdir / "clusters"
        export_all_levels(root, clusters_dir, levels_cfg=levels_cfg, global_cfg=global_cfg)
        _log("Exported per-level cluster HTML")

        ga_json = export_ga_progress_json(root, outdir / "ga_progress.json", root_id=job_id)
        export_ga_progress_html(ga_json, outdir / "ga_progress.html")
        _log("Exported GA progress reports")

        # Persist GA progress JSON to DB (used by dynamic report renderer)
        try:
            db = SessionLocal()
            try:
                save_ga_progress(db, job_id=job_id, ga_progress_json=json.dumps(ga_json))
                _log("Saved GA progress to database")
            finally:
                db.close()
        except Exception as ex:
            _log(f"WARN: failed to save GA progress to DB: {ex}")

        analysis_json = analyze_tree_to_json(root, outdir / "analysis.json", root_id=job_id)
        export_analysis_html(analysis_json, outdir / "analysis.html")
        _log("Exported analysis reports")

        export_tree_to_html(root, outdir / "tree.html", root_id=job_id)
        export_fragments_to_fasta(root, outdir / "fragments.fasta", root_id=job_id)
        _log("Exported tree HTML and fragments FASTA")

        # --- Persist tree JSON into the Design row for on-the-fly rendering ---
        try:
            if tree_json_path and tree_json_path.exists():
                payload = tree_json_path.read_text(encoding="utf-8")
                db = SessionLocal()
                try:
                    if attach_tree_to_design(db, job_id=job_id, tree_json=payload):
                        _log("Stored design tree in database")
                finally:
                    db.close()
        except Exception as ex:
            _log(f"WARN: failed to store tree.json to DB: {ex}")

        # --- Build a simple index (transient) and announce result URL ---
        idx = write_run_index_simple(outdir, job_id, reports_public_base=reports_public_base)
        _log(f"Generated index.html → {idx}")

        rc = 0  # success

    except Exception as ex:
        # Emit stack trace for diagnosis
        tb = traceback.format_exc()
        _log("ERROR: pipeline failed")
        _log(tb)
        rc = 1

    finally:
        # Update Run row with final status and dynamic report URL (if success)
        try:
            db = SessionLocal()
            try:
                report_url = f"{reports_public_base.rstrip('/')}/{job_id}/index.html" if rc == 0 else None
                mark_run_completed(db, job_id=job_id, report_url=report_url, exit_code=rc)
            finally:
                db.close()
        except Exception as ex:
            _log(f"WARN: failed to mark run completed: {ex}")

        # Emit final, stable lines the UI relies on
        if rc == 0:
            _log(f"RESULT: {reports_public_base.rstrip('/')}/{job_id}/index.html")
        else:
            _log("RESULT: (none)")
        _log(f"EXIT: {rc}")

    return rc
