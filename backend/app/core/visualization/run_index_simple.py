# File: backend/app/core/visualization/run_index_simple.py
# Version: v1.3.0
"""
Modern, grouped index page for a run, with a light "Apple-like" aesthetic.

Groups:
- FASTA      → input.fasta, fragments.fasta, oligos.fasta, assembly.gb
- Tables     → fragments.csv, oligos.csv
- Parameters → globals.json, levels.json
- Visuals    → analysis.html, tree.html, GA progress, clusters
- Bundle     → Download all (ZIP)
"""
from __future__ import annotations

from pathlib import Path


def _exists(p: Path) -> bool:
    try:
        return p.is_file()
    except Exception:
        return False


def _collect_clusters(clusters_dir: Path) -> list[str]:
    items: list[str] = []
    if clusters_dir.is_dir():
        for p in sorted(clusters_dir.glob("*.html")):
            items.append(p.name)
    return items


def write_run_index_simple(outdir: Path, job_id: str, *, reports_public_base: str = "/reports") -> None:
    base = f"{reports_public_base}/{job_id}"
    clusters = _collect_clusters(outdir / "clusters")

    def link(name: str, fname: str) -> str:
        p = outdir / fname
        if _exists(p):
            return f'<a class="item" href="{base}/{fname}" download>{name}</a>'
        return f'<span class="item disabled">{name}</span>'

    def link_nodl(name: str, fname: str) -> str:
        p = outdir / fname
        if _exists(p):
            return f'<a class="item" href="{base}/{fname}">{name}</a>'
        return f'<span class="item disabled">{name}</span>'

    # Bundle link is always available because the API will generate it on the fly.
    def link_bundle() -> str:
        return f'<a class="btn primary" href="{base}/bundle.zip">Download all (ZIP)</a>'

    clusters_html = ""
    if clusters:
        cluster_links = "\n".join(
            f'<a class="pill" href="{base}/clusters/{fn}">{fn}</a>' for fn in clusters
        )
        clusters_html = f"""
        <div class="card">
          <div class="card-title">Clusters</div>
          <div class="cluster-list">{cluster_links}</div>
        </div>
        """

    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>CornStructor report — {job_id}</title>
<style>
  :root {{
    --bg: #f5f5f7;
    --card: #ffffff;
    --text: #1d1d1f;
    --sub: #6e6e73;
    --link: #0071e3;
    --border: #e5e5ea;
    --pill: #f2f2f7;
    --btnbg: #e6f5ff;
    --btnbg-hover: #d9efff;
    --btnborder: #b6dfff;
  }}
  * {{ box-sizing: border-box; }}
  body {{
    margin: 0; padding: 24px;
    background: var(--bg);
    color: var(--text);
    font: 14px/1.45 -apple-system, BlinkMacSystemFont, "SF Pro Text", "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
  }}
  .wrap {{ max-width: 1100px; margin: 0 auto; }}
  header {{ margin: 8px 0 20px 0; display:flex; align-items:center; justify-content:space-between; gap:12px; flex-wrap:wrap; }}
  h1 {{
    font-weight: 700; letter-spacing: -.02em; margin: 0 0 4px 0;
    font-size: 22px;
  }}
  .sub {{ color: var(--sub); font-size: 13px; }}
  .grid {{
    display: grid; gap: 14px;
    grid-template-columns: repeat(auto-fill, minmax(280px, 1fr));
    align-items: start;
  }}
  .card {{
    background: var(--card);
    border: 1px solid var(--border);
    border-radius: 14px;
    padding: 14px;
    box-shadow: 0 8px 24px rgba(0,0,0,.05);
  }}
  .card-title {{
    font-weight: 650; margin-bottom: 10px;
  }}
  .item-list {{ display: flex; flex-direction: column; gap: 8px; }}
  a.item {{ color: var(--link); text-decoration: none; font-weight: 500; }}
  a.item:hover {{ text-decoration: underline; }}
  .item.disabled {{ color: #bbb; text-decoration: none; cursor: not-allowed; }}
  .pill {{
    display: inline-block; padding: 6px 10px; background: var(--pill);
    border: 1px solid var(--border); border-radius: 999px; text-decoration: none;
    color: var(--text); margin: 4px 6px 0 0; font-size: 12px;
  }}
  .pill:hover {{ border-color: #c7c7cc; }}
  footer {{ margin-top: 22px; color: var(--sub); font-size: 12px; }}

  .btn {{
    display: inline-block;
    padding: 8px 14px;
    border-radius: 999px;
    border: 1px solid var(--btnborder);
    background: var(--btnbg);
    color: #0b5cab;
    text-decoration: none;
    font-weight: 600;
    transition: background .15s ease, transform .05s ease;
  }}
  .btn:hover {{ background: var(--btnbg-hover); }}
  .btn:active {{ transform: translateY(1px); }}
  .btn.primary {{
    background: #d1fae5; border-color: #a7f3d0; color: #065f46;
  }}
  .btn.primary:hover {{ background: #a7f3d0; }}
</style>
</head>
<body>
  <div class="wrap">
    <header>
      <div>
        <h1>Report for run <code>{job_id}</code></h1>
        <div class="sub">Quick access to exports and visualizations.</div>
      </div>
      <div>
        {link_bundle()}
      </div>
    </header>

    <section class="grid">

      <div class="card">
        <div class="card-title">FASTA</div>
        <div class="item-list">
          {link("Input FASTA", "input.fasta")}
          {link("All fragments (FASTA)", "fragments.fasta")}
          {link("Oligos only (FASTA)", "oligos.fasta")}
          {link("GenBank (assembly.gb)", "assembly.gb")}
        </div>
      </div>

      <div class="card">
        <div class="card-title">Tables</div>
        <div class="item-list">
          {link("All fragments (CSV)", "fragments.csv")}
          {link("Oligos only (CSV)", "oligos.csv")}
        </div>
      </div>

      <div class="card">
        <div class="card-title">Parameters</div>
        <div class="item-list">
          {link("Globals (JSON)", "globals.json")}
          {link("Levels (JSON)", "levels.json")}
        </div>
      </div>

      <div class="card">
        <div class="card-title">Visualizations</div>
        <div class="item-list">
          {link_nodl("Analysis overview", "analysis.html")}
          {link_nodl("Assembly tree", "tree.html")}
          {link_nodl("GA progress (HTML)", "ga_progress.html")}
          {link("GA progress (JSON)", "ga_progress.json")}
          {link("Tree (JSON)", "tree.json")}
        </div>
      </div>

      {clusters_html}

    </section>

    <footer>
      CornStructor — generated locally. All files are transient per-run artifacts.
    </footer>
  </div>
</body>
</html>
"""
    (outdir / "index.html").write_text(html, encoding="utf-8")
