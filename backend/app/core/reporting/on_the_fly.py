# File: backend/app/core/reporting/on_the_fly.py
# Version: v1.0.0
"""
On-the-fly report rendering helpers.

Given a Design row from the DB, produce all report artifacts into a temporary
folder using the same exporters as the batch pipeline, but without persisting
to the shared /reports. This supports dynamic serving after container restarts.

Strategy (Path B)
-----------------
- If Design.params_json contains snapshots of `globals` and `levels`, we write
  them to temp files and load config from them; otherwise we fall back to the
  current settings.LEVELS_PATH and GLOBALS_PATH.
- We reconstruct the tree deterministically by re-running HierarchicalAssembler
  on the stored sequence and config.
- If `design.tree_json` exists, we also write it out as tree.json (so index links
  work) and we prefer it for parity with the original run.
- If `design.ga_progress_json` exists, we use it to render ga_progress.html and
  ga_progress.json; otherwise we export GA progress from the reconstructed tree.
- We export analysis.html/json, cluster HTML, tree.html and fragments.fasta.
- Finally, we build an index.html that links to everything present.

Caching
-------
- A tiny LRU cache (MAX_CACHE) keeps a handful of recently rendered temp dirs
  to avoid recomputation when a user clicks around within the same run.
"""
from __future__ import annotations

import json
import tempfile
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional

from backend.app.core.config import settings
from backend.app.config.config_level import load_levels_config, LevelConfig
from backend.app.config.config_global import load_global_config, GlobalConfig
from backend.app.core.assembly.hierarchical_assembler import HierarchicalAssembler

from backend.app.core.export.json_exporter import export_tree_to_json
from backend.app.core.export.analysis_exporter import analyze_tree_to_json
from backend.app.core.export.ga_progress_exporter import export_ga_progress_json
from backend.app.core.export.fasta_exporter import export_fragments_to_fasta

from backend.app.core.visualization.cluster_html_report import export_all_levels
from backend.app.core.visualization.ga_progress_html_report import export_ga_progress_html
from backend.app.core.visualization.analysis_html_report import export_analysis_html
from backend.app.core.visualization.tree_html_exporter import export_tree_to_html
from backend.app.core.visualization.run_index_simple import write_run_index_simple

from backend.app.db.models import Design

MAX_CACHE = 8
RENDERED_INDEX = "index.html"


@dataclass
class RenderEntry:
    dir: Path


_cache: "OrderedDict[str, RenderEntry]" = OrderedDict()


def _lru_get(job_id: str) -> Optional[RenderEntry]:
    entry = _cache.get(job_id)
    if entry:
        # move to end (recently used)
        _cache.move_to_end(job_id)
    return entry


def _lru_put(job_id: str, entry: RenderEntry) -> None:
    _cache[job_id] = entry
    _cache.move_to_end(job_id)
    while len(_cache) > MAX_CACHE:
        old_job, old_entry = _cache.popitem(last=False)
        # best-effort cleanup
        try:
            if old_entry.dir.exists():
                for p in old_entry.dir.rglob("*"):
                    try:
                        p.unlink()
                    except Exception:
                        pass
                old_entry.dir.rmdir()
        except Exception:
            pass


def _write_config_snapshots(params_json: Optional[str], tmpdir: Path) -> tuple[Path, Path]:
    """
    From Design.params_json, extract 'globals' and 'levels' snapshots if present
    and write them to tmpdir. Returns (levels_path, globals_path).
    If snapshots are missing, return the current settings paths.
    """
    try:
        data = json.loads(params_json) if params_json else {}
    except Exception:
        data = {}
    gl = data.get("globals")
    lv = data.get("levels")

    if gl is not None and lv is not None:
        gl_path = tmpdir / "globals.snapshot.json"
        lv_path = tmpdir / "levels.snapshot.json"
        gl_path.write_text(json.dumps(gl), encoding="utf-8")
        lv_path.write_text(json.dumps(lv), encoding="utf-8")
        return lv_path, gl_path

    # Fallback to current config files on disk
    return settings.LEVELS_PATH, settings.GLOBALS_PATH


def _reconstruct_and_export(design: Design, outdir: Path, job_id: str) -> None:
    """Rebuild the tree and export all artifacts to outdir."""
    # 1) Save sequence as FASTA input for reproducibility (optional)
    (outdir / "input.fasta").write_text(f">seq\n{design.sequence}\n", encoding="utf-8")

    # 2) Resolve config: prefer snapshots inside params_json
    lv_path, gl_path = _write_config_snapshots(design.params_json, outdir)

    # 3) Load config and build tree
    levels_cfg: Dict[int, LevelConfig] = load_levels_config(lv_path)
    global_cfg: GlobalConfig = load_global_config(str(gl_path))
    assembler = HierarchicalAssembler(levels_cfg, global_cfg)
    root = assembler.build(full_seq=design.sequence, root_id=job_id)

    # 4) Export core artifacts
    export_tree_to_json(root, outdir / "tree.json", root_id=job_id)
    clusters_dir = outdir / "clusters"
    export_all_levels(root, clusters_dir, levels_cfg=levels_cfg, global_cfg=global_cfg)

    # analysis (json + html)
    analysis_json = analyze_tree_to_json(root, outdir / "analysis.json", root_id=job_id)
    export_analysis_html(analysis_json, outdir / "analysis.html")

    # tree html + fragments
    export_tree_to_html(root, outdir / "tree.html", root_id=job_id)
    export_fragments_to_fasta(root, outdir / "fragments.fasta", root_id=job_id)

    # GA progress: prefer stored, else compute fresh
    if design.ga_progress_json:
        (outdir / "ga_progress.json").write_text(design.ga_progress_json, encoding="utf-8")
        try:
            data = json.loads(design.ga_progress_json)
        except Exception:
            data = {}
        export_ga_progress_html(data, outdir / "ga_progress.html")
    else:
        ga_json = export_ga_progress_json(root, outdir / "ga_progress.json", root_id=job_id)
        export_ga_progress_html(ga_json, outdir / "ga_progress.html")

    # 5) Robust index
    write_run_index_simple(outdir, job_id, reports_public_base="/reports")


def ensure_rendered(*, job_id: str, design: Design) -> Path:
    """
    Ensure that a fully rendered artifact directory exists for `job_id`.
    Returns the output directory path.
    """
    existing = _lru_get(job_id)
    if existing and existing.dir.exists() and (existing.dir / RENDERED_INDEX).is_file():
        return existing.dir

    tmp = Path(tempfile.mkdtemp(prefix=f"cornstructor_{job_id}_"))
    _reconstruct_and_export(design, tmp, job_id)
    entry = RenderEntry(dir=tmp)
    _lru_put(job_id, entry)
    return tmp
