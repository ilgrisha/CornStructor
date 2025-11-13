# File: backend/app/core/visualization/report_index_html.py
# Version: v0.2.0

"""
Generate a single HTML "index" page with links to all outputs for a given run:
- Tree (HTML)
- Analysis (HTML)
- GA Progress (HTML)
- Per-level cluster reports (directory of HTMLs)
- Download link for fragments FASTA
- Download link for fragments CSV (sequence_id,sequence,sequence_length)
- (Optional) raw JSON assets, if present

Usage:
    export_run_index_html(root_id, outdir)

Writes:
    <outdir>/<root_id}_index.html
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Tuple
import datetime
import html


def _fmt_bytes(n: int) -> str:
    for unit in ["B", "KB", "MB", "GB"]:
        if n < 1024.0:
            return f"{n:.0f} {unit}"
        n /= 1024.0
    return f"{n:.0f} TB"


def _file_info(p: Path) -> Tuple[str, str]:
    """Return (size_str, mtime_str) for an existing file path."""
    try:
        st = p.stat()
        size = _fmt_bytes(st.st_size)
        mtime = datetime.datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d %H:%M")
        return size, mtime
    except Exception:
        return ("", "")


def export_run_index_html(root_id: str, outdir: Path) -> Path:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Expected top-level artifacts
    fasta = outdir / f"{root_id}_fragments.fasta"
    csvf  = outdir / f"{root_id}_fragments.csv"
    tree_html = outdir / f"{root_id}_tree.html"
    analysis_html = outdir / f"{root_id}_analysis.html"
    ga_html = outdir / f"{root_id}_ga_progress.html"

    # Optional JSONs
    tree_json = outdir / f"{root_id}_tree.json"
    analysis_json = outdir / f"{root_id}_analysis.json"
    ga_json = outdir / f"{root_id}_ga_progress.json"

    # Cluster reports directory (may contain multiple HTML files)
    cluster_dir = outdir / f"{root_id}_cluster_reports"

    quick_links: List[str] = []
    if tree_html.exists():
        sz, mt = _file_info(tree_html)
        quick_links.append(
            f'<li><a href="{tree_html.name}">Assembly Tree</a> '
            f'<span class="meta">({sz}, {mt})</span></li>'
        )
    if analysis_html.exists():
        sz, mt = _file_info(analysis_html)
        quick_links.append(
            f'<li><a href="{analysis_html.name}">Assembly Analysis</a> '
            f'<span class="meta">({sz}, {mt})</span></li>'
        )
    if ga_html.exists():
        sz, mt = _file_info(ga_html)
        quick_links.append(
            f'<li><a href="{ga_html.name}">GA Progress</a> '
            f'<span class="meta">({sz}, {mt})</span></li>'
        )

    # FASTA & CSV downloads
    fasta_link = ""
    if fasta.exists():
        sz, mt = _file_info(fasta)
        fasta_link = (
            f'<a class="btn" href="{fasta.name}" download>Download fragments FASTA</a> '
            f'<span class="meta">({sz}, {mt})</span>'
        )

    csv_link = ""
    if csvf.exists():
        sz, mt = _file_info(csvf)
        csv_link = (
            f'<a class="btn" href="{csvf.name}" download>Download sequences CSV</a> '
            f'<span class="meta">({sz}, {mt})</span>'
        )

    # Optional JSON assets section
    json_items: List[str] = []
    for p in (tree_json, analysis_json, ga_json):
        if p.exists():
            sz, mt = _file_info(p)
            json_items.append(
                f'<li><a href="{p.name}">{p.name}</a> '
                f'<span class="meta">({sz}, {mt})</span></li>'
            )

    # Cluster report HTMLs
    cluster_items: List[str] = []
    if cluster_dir.exists() and cluster_dir.is_dir():
        htmls = sorted(cluster_dir.glob("*.html")) or sorted(cluster_dir.rglob("*.html"))
        for hp in htmls:
            rel = hp.relative_to(outdir).as_posix()
            sz, mt = _file_info(hp)
            title = hp.stem
            cluster_items.append(
                f'<li><a href="{rel}">{html.escape(title)}</a> '
                f'<span class="meta">({sz}, {mt})</span></li>'
            )

    # Build HTML
    root_h = html.escape(root_id)
    clusters_block = (
        "<ul>\n" + "\n".join(cluster_items) + "\n</ul>"
        if cluster_items else "<p class='muted'>No per-level cluster reports found.</p>"
    )
    jsons_block = (
        "<ul>\n" + "\n".join(json_items) + "\n</ul>"
        if json_items else "<p class='muted'>No JSON assets found.</p>"
    )
    quick_block = (
        "<ul>\n" + "\n".join(quick_links) + "\n</ul>"
        if quick_links else "<p class='muted'>No HTML visualizations found yet.</p>"
    )
    fasta_block = fasta_link if fasta_link else "<p class='muted'>FASTA file not found.</p>"
    csv_block = csv_link if csv_link else "<p class='muted'>CSV file not found.</p>"

    html_doc = f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8" />
  <title>Assembly Reports: {root_h}</title>
  <style>
    body {{ font-family: -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif; margin: 22px; }}
    h1, h2 {{ margin: 0.2em 0 0.4em 0; }}
    section {{ border: 1px solid #ddd; padding: 12px 14px; margin: 14px 0; border-radius: 6px; }}
    ul {{ margin: 6px 0 0 18px; }}
    .meta {{ color: #666; font-size: 12px; margin-left: 6px; }}
    .muted {{ color: #777; }}
    .btn {{
      display: inline-block; padding: 8px 12px; border-radius: 6px; border: 1px solid #ccc;
      text-decoration: none; color: #0a5; font-weight: 600; background: #f5fff8;
    }}
    .btn:hover {{ background: #ecfff2; }}
  </style>
</head>
<body>
  <h1>Assembly Reports: {root_h}</h1>

  <section>
    <h2>Quick links</h2>
    {quick_block}
  </section>

  <section>
    <h2>Fragments FASTA</h2>
    {fasta_block}
  </section>

  <section>
    <h2>Sequences CSV</h2>
    {csv_block}
  </section>

  <section>
    <h2>Per-level Cluster Reports</h2>
    {clusters_block}
  </section>

  <section>
    <h2>Raw JSON Assets</h2>
    {jsons_block}
  </section>
</body>
</html>
"""
    out_path = outdir / f"{root_id}_index.html"
    out_path.write_text(html_doc, encoding="utf-8")
    return out_path
