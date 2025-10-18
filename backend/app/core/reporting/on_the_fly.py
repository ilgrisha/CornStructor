# File: backend/app/core/reporting/on_the_fly.py
# Version: v1.5.0
"""
On-the-fly report rendering helpers (simplified).

v1.5.0
- Write JSON artifacts in friendly format (indent=2, sorted keys):
    globals.json, levels.json, tree.json, analysis.json, ga_progress.json
- Delegate CSV/FASTA/GenBank to core/export modules.
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

from backend.app.core.export.fasta_oriented_exporter import (
    export_fragments_fasta_oriented,
    export_oligos_fasta,
)
from backend.app.core.export.csv_exporter import export_csvs
from backend.app.core.export.genbank_exporter import export_genbank_from_tree

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
        _cache.move_to_end(job_id)
    return entry


def _lru_put(job_id: str, entry: RenderEntry) -> None:
    _cache[job_id] = entry
    _cache.move_to_end(job_id)
    while len(_cache) > MAX_CACHE:
        old_job, old_entry = _cache.popitem(last=False)
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


def _pretty_write(path: Path, data) -> None:
    """Write JSON with indent=2, sorted keys."""
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def _read_json(path: Path):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def _write_configs(params_json: Optional[str], tmpdir: Path) -> tuple[Path, Path]:
    """
    Persist ONLY levels.json and globals.json into tmpdir.
    Prefer snapshots embedded in Design.params_json; otherwise read from disk.
    Always pretty-print.
    """
    gl_path = tmpdir / "globals.json"
    lv_path = tmpdir / "levels.json"

    try:
        data = json.loads(params_json) if params_json else {}
    except Exception:
        data = {}

    gl = data.get("globals")
    lv = data.get("levels")

    if gl is None:
        try:
            gl = json.loads(Path(settings.GLOBALS_PATH).read_text(encoding="utf-8"))
        except Exception:
            gl = {}
    if lv is None:
        try:
            lv = json.loads(Path(settings.LEVELS_PATH).read_text(encoding="utf-8"))
        except Exception:
            lv = {}

    _pretty_write(gl_path, gl)
    _pretty_write(lv_path, lv)
    return lv_path, gl_path


def _reconstruct_and_export(design: Design, outdir: Path, job_id: str) -> None:
    # Input FASTA
    (outdir / "input.fasta").write_text(f">seq\n{design.sequence}\n", encoding="utf-8")

    # Configs
    lv_path, gl_path = _write_configs(design.params_json, outdir)

    levels_cfg: Dict[int, LevelConfig] = load_levels_config(lv_path)
    global_cfg: GlobalConfig = load_global_config(str(gl_path))
    assembler = HierarchicalAssembler(levels_cfg, global_cfg)
    root = assembler.build(full_seq=design.sequence, root_id=job_id)

    # JSON & HTML reports
    tree_json = export_tree_to_json(root, outdir / "tree.json", root_id=job_id) or _read_json(outdir / "tree.json")
    if tree_json is not None:
        _pretty_write(outdir / "tree.json", tree_json)

    clusters_dir = outdir / "clusters"
    export_all_levels(root, clusters_dir, levels_cfg=levels_cfg, global_cfg=global_cfg)

    analysis_json = analyze_tree_to_json(root, outdir / "analysis.json", root_id=job_id) or _read_json(outdir / "analysis.json")
    if analysis_json is not None:
        _pretty_write(outdir / "analysis.json", analysis_json)
    export_analysis_html(analysis_json or {}, outdir / "analysis.html")

    export_tree_to_html(root, outdir / "tree.html", root_id=job_id)

    # FASTA exports (strand-aware)
    export_fragments_fasta_oriented(root, outdir / "fragments.fasta", root_id=job_id)
    export_oligos_fasta(root, outdir / "oligos.fasta", root_id=job_id)

    # CSV exports (strand-aware + overlap Tm)
    export_csvs(
        root,
        outdir,
        root_id=job_id,
        tm_method=global_cfg.tm_method,
        tm_params=global_cfg.tm,
    )

    # GenBank export (oligos→oligos, fragments→synthons)
    export_genbank_from_tree(root, outdir / "assembly.gb", record_id=job_id, definition="CornStructor assembly.")

    # GA progress artifacts
    if design.ga_progress_json:
        _pretty_write(outdir / "ga_progress.json", json.loads(design.ga_progress_json))
        export_ga_progress_html(json.loads(design.ga_progress_json), outdir / "ga_progress.html")
    else:
        ga_json = export_ga_progress_json(root, outdir / "ga_progress.json", root_id=job_id) or _read_json(outdir / "ga_progress.json")
        if ga_json is not None:
            _pretty_write(outdir / "ga_progress.json", ga_json)
        export_ga_progress_html(ga_json or {}, outdir / "ga_progress.html")

    # Summary index
    write_run_index_simple(outdir, job_id, reports_public_base="/reports")


def ensure_rendered(*, job_id: str, design: Design) -> Path:
    existing = _lru_get(job_id)
    if existing and existing.dir.exists() and (existing.dir / RENDERED_INDEX).is_file():
        return existing.dir

    tmp = Path(tempfile.mkdtemp(prefix=f"cornstructor_{job_id}_"))
    _reconstruct_and_export(design, tmp, job_id)
    entry = RenderEntry(dir=tmp)
    _lru_put(job_id, entry)
    return tmp
