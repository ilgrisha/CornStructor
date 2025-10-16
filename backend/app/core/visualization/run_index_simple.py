# File: backend/app/core/visualization/run_index_simple.py
# Version: v0.1.0
"""
Simple, robust HTML index generator for a CornStructor run.

- Scans the job's output directory for known artifacts and lists whatever exists.
- No manifest required; no assumptions about a prior CLI.
- Always writes an index.html so links do not 404.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List


def _rel(p: Path, base: Path) -> str:
    try:
        return str(p.relative_to(base))
    except Exception:
        return p.name


def _collect_artifacts(outdir: Path) -> Dict[str, List[str]]:
    """
    Return mapping of artifact category -> list of relative paths found.
    """
    artifacts: Dict[str, List[str]] = {
        "HTML": [],
        "JSON": [],
        "FASTA": [],
        "Clusters": [],
        "Other": [],
    }

    # Singles we commonly produce
    singles = {
        "HTML": ["tree.html", "analysis.html", "ga_progress.html"],
        "JSON": ["tree.json", "analysis.json", "ga_progress.json"],
        "FASTA": ["fragments.fasta"],
    }
    for cat, names in singles.items():
        for n in names:
            p = outdir / n
            if p.is_file():
                artifacts[cat].append(_rel(p, outdir))

    # Clusters (per-level cluster HTML)
    clusters_dir = outdir / "clusters"
    if clusters_dir.is_dir():
        for p in sorted(clusters_dir.glob("*.html")):
            artifacts["Clusters"].append(_rel(p, outdir))

    # Anything else useful-looking (fallback)
    for p in sorted(outdir.glob("*")):
        if p.name in {"index.html", "input.fasta", "clusters"}:
            continue
        if p.is_file() and p.suffix.lower() not in {".py", ".log"}:
            rel = _rel(p, outdir)
            if not any(rel in lst for lst in artifacts.values()):
                artifacts["Other"].append(rel)

    # Remove empty categories
    return {k: v for k, v in artifacts.items() if v}

def write_run_index_simple(outdir: Path, job_id: str, reports_public_base: str = "/reports") -> Path:
    """
    Write a minimal, durable index.html that links to present artifacts.

    Returns the path to index.html.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    index_path = outdir / "index.html"

    artifacts = _collect_artifacts(outdir)
    base = reports_public_base.rstrip("/") + f"/{job_id}"

    def _section(title: str, items: List[str]) -> str:
        lis = "\n".join(f"<li><a href='{base}/{item}' target='_blank' rel='noopener'>{item}</a></li>"
                        for item in items)
        return f"<section><h2>{title}</h2><ul>{lis}</ul></section>"

    if artifacts:
        sections = "\n".join(_section(title, items) for title, items in artifacts.items())
        msg = ""
    else:
        sections = ""
        msg = "<p><em>No artifacts were detected for this run.</em></p>"

    html = f"""<!doctype html>
<meta charset="utf-8">
<title>CornStructor Report — {job_id}</title>
<style>
  body {{ font-family: system-ui, -apple-system, Segoe UI, Roboto, sans-serif; margin: 2rem; }}
  header {{ margin-bottom: 1rem; }}
  h1 {{ font-size: 1.6rem; margin: 0 0 .25rem 0; }}
  h2 {{ font-size: 1.1rem; margin: 1rem 0 .25rem 0; }}
  ul {{ margin: .25rem 0 1rem 1.2rem; }}
  .meta {{ color: #555; font-size: .9rem; }}
  .back a {{ text-decoration: none; }}
</style>
<header>
  <h1>Run {job_id}</h1>
  <div class="meta">Artifacts under <code>{base}</code></div>
</header>
{msg}
{sections}
<div class="back"><a href="/">&larr; Back to CornStructor</a></div>
"""
    index_path.write_text(html, encoding="utf-8")
    return index_path
