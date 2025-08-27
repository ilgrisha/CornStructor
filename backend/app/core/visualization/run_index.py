# File: backend/app/core/visualization/run_index.py
# Version: v0.4.0
"""
Run index HTML generator — Apple-style UI (no embedded history)

This module scans a job output directory for report artifacts and emits a single
modern ``index.html`` (SF font stack, glassy cards, subtle gradients, dark mode).

v0.4.0
------
- Removed the embedded "Previous reports" table and the local `runs.db`.
  History now lives in the main app database and is accessible from the app's
  main screen.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List
import html

@dataclass
class RunArtifact:
    """A single run artifact for the results page."""
    kind: str         # 'analysis' | 'ga' | 'clusters' | future kinds
    title: str        # User-facing label
    rel_path: str     # HREF relative to the job outdir (e.g., 'file.html' or 'dir/')

# ---------- Artifact Discovery ----------

def _find_artifacts(outdir: Path) -> List[RunArtifact]:
    """Probe a run outdir for known artifacts."""
    arts: List[RunArtifact] = []

    # Analysis HTML
    for p in sorted(outdir.glob("*_analysis.html")):
        arts.append(RunArtifact("analysis", f"Analysis · {p.stem}", p.name))

    # GA summary HTML
    for p in sorted(outdir.glob("*_ga_progress.html")):
        arts.append(RunArtifact("ga", f"GA Progress · {p.stem}", p.name))

    # Cluster reports directory
    for d in sorted(outdir.glob("*_cluster_reports")):
        if d.is_dir():
            arts.append(RunArtifact("clusters", f"Cluster reports · {d.name}", f"{d.name}/"))

    return arts

# ---------- HTML ----------

def _render_html(*, job_id: str, artifacts: List[RunArtifact], reports_public_base: str) -> str:
    """Render a minimal index page with links to artifacts."""
    # Build cards
    items_html = []
    for a in artifacts:
        icon = "🔎" if a.kind == "analysis" else ("🧬" if a.kind == "ga" else "🗂️")
        items_html.append(
            f"""
<a class='card' href='{html.escape(a.rel_path)}' target='_blank' rel='noopener'>
    <div class='icon'>{icon}</div>
    <div class='title'>{html.escape(a.title)}</div>
    <div class='chev'>→</div>
</a>"""
        )
    if not items_html:
        items_html.append("<div class='empty'>No artifacts found for this run.</div>")

    # Link back to app home
    back_link = f"<a class='back' href='/'>&larr; Back to CornStructor</a>"

    return f"""<!doctype html>
<html lang='en'>
<meta charset='utf-8'>
<meta name='viewport' content='width=device-width, initial-scale=1'>
<title>CornStructor · Run {html.escape(job_id)}</title>
<style>
  :root {{ --bg: #0b0f14; --fg: #e6edf3; --muted: #9fb0c0; --card:#10161d; --border:#263241; }}
  @media (prefers-color-scheme: light) {{ :root {{ --bg:#f8fafc; --fg:#0b1220; --muted:#4b5563; --card:#ffffff; --border:#e5e7eb; }} }}
  * {{ box-sizing:border-box; }}
  body {{ margin:0; padding:3rem 1rem; background:var(--bg); color:var(--fg); font:16px/1.5 -apple-system, BlinkMacSystemFont, 'SF Pro Text', Inter, Segoe UI, Roboto, sans-serif; }}
  .container {{ max-width: 860px; margin: 0 auto; }}
  .header {{ display:flex; align-items:center; justify-content:space-between; margin-bottom:1rem; }}
  h1 {{ margin:0; font-size:1.4rem; letter-spacing:.2px; }}
  .grid {{ display:grid; grid-template-columns: repeat(auto-fill, minmax(260px, 1fr)); gap: 12px; }}
  .card {{ display:flex; align-items:center; gap:12px; padding:14px 16px; background:var(--card); border:1px solid var(--border); border-radius:14px; text-decoration:none; color:inherit; transition: transform .06s ease; }}
  .card:hover {{ transform: translateY(-1px); }}
  .icon {{ font-size:20px; }}
  .title {{ font-weight:600; }}
  .chev {{ margin-left:auto; opacity:.6; }}
  .empty {{ opacity:.7; padding: 1rem; border:1px dashed var(--border); border-radius:12px; }}
  .back {{ display:inline-block; margin-top: 1.25rem; color: inherit; opacity:.8; text-decoration:none; }}
</style>
<div class='container'>
  <div class='header'>
    <h1>Run <code>{html.escape(job_id)}</code></h1>
  </div>
  <div class='grid'>
    {''.join(items_html)}
  </div>
  {back_link}
</div>
</html>"""

# ---------- Public API ----------

def write_run_index(outdir: Path, job_id: str, *, reports_public_base: str = "/reports") -> Path | None:
    """Generate a minimal `index.html` inside *outdir*.

    Returns the path to the generated ``index.html`` or ``None`` if *outdir* does not exist.
    """
    if not outdir.exists():
        return None

    # Discover artifacts for this run
    artifacts = _find_artifacts(outdir)

    # Render and write HTML
    html_text = _render_html(job_id=job_id, artifacts=artifacts, reports_public_base=reports_public_base)
    index_path = outdir / "index.html"
    index_path.write_text(html_text, encoding="utf-8")
    return index_path
