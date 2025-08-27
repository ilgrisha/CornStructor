# File: backend/app/core/visualization/run_index.py
# Version: v0.3.0
"""
Run index HTML generator — Apple-style UI + Local History

This module scans a job output directory for report artifacts and emits a single
modern, Apple-inspired ``index.html`` (SF font stack, glassy cards, subtle gradients,
automatic dark mode). It now also persists a *local history* of runs in a small
SQLite database and exposes an interactive "Previous reports" table.

What's new in v0.3.0
--------------------
- Creates/updates a lightweight SQLite database (``runs.db``) next to your run folders.
- Persists the current run {job_id, report_url, created_at}.
- Renders an Apple-style page with a button **"Open previous reports"** that toggles
  a searchable, sortable table of recent runs (client-side).
- Adds small UX niceties (copy link, quick filter, keyboard shortcuts).

Artifacts detected:
- ``*_analysis.html`` (overall analysis report)
- ``*_ga_progress.html`` (GA progress summary)
- ``*_cluster_reports/`` (per-level cluster visualizations)

Public path assumptions:
- Your static reports are served at ``{reports_public_base}/{job_id}/``.
- All artifact links on the page are *relative* to that directory.

API
---
write_run_index(outdir: Path, job_id: str, reports_public_base: str = "/reports") -> Optional[Path]
    Generate/overwrite ``index.html`` inside *outdir* and return its path.

Notes
-----
- The history database path defaults to ``outdir.parent / 'runs.db'``.
- If SQLite is unavailable or not writable, the page still renders; history just won't show.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple
import html
import json
import sqlite3
from datetime import datetime


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

    # Cluster report directories
    for d in sorted(outdir.glob("*_cluster_reports")):
        if d.is_dir():
            arts.append(RunArtifact("clusters", f"Cluster Reports · {d.name}", d.name + "/"))

    return arts


# ---------- Local History (SQLite) ----------

def _db_path_for(outdir: Path) -> Path:
    """Default SQLite DB path (sits next to job folders)."""
    return outdir.parent / "runs.db"


def _db_init(conn: sqlite3.Connection) -> None:
    """Ensure the 'runs' table exists."""
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS runs (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            job_id TEXT UNIQUE NOT NULL,
            report_url TEXT NOT NULL,
            created_at TEXT NOT NULL
        )
        """
    )
    conn.commit()


def _db_upsert_run(conn: sqlite3.Connection, job_id: str, report_url: str) -> None:
    """Insert or update a run row."""
    now = datetime.utcnow().isoformat(timespec="seconds") + "Z"
    conn.execute(
        """
        INSERT INTO runs (job_id, report_url, created_at)
        VALUES (?, ?, ?)
        ON CONFLICT(job_id) DO UPDATE SET
            report_url=excluded.report_url,
            created_at=excluded.created_at
        """,
        (job_id, report_url, now),
    )
    conn.commit()


def _db_fetch_recent(conn: sqlite3.Connection, limit: int = 100) -> List[Tuple[str, str, str]]:
    """Fetch recent runs as (job_id, report_url, created_at) newest first."""
    cur = conn.execute(
        "SELECT job_id, report_url, created_at FROM runs ORDER BY datetime(created_at) DESC LIMIT ?",
        (limit,),
    )
    return list(cur.fetchall())


def _load_history_json(db_path: Path, current_job_id: str) -> str:
    """
    Read recent history and return a JSON string with rows:
      [{job_id, report_url, created_at, is_current}]
    If the DB is not available, returns '[]'.
    """
    try:
        conn = sqlite3.connect(db_path)
        try:
            _db_init(conn)
            rows = _db_fetch_recent(conn)
        finally:
            conn.close()
    except Exception:
        return "[]"

    data = [
        {
            "job_id": job_id,
            "report_url": report_url,
            "created_at": created_at,
            "is_current": job_id == current_job_id,
        }
        for (job_id, report_url, created_at) in rows
    ]
    return json.dumps(data)


# ---------- Page Rendering ----------

def _render_css() -> str:
    """Apple-inspired, responsive CSS with light/dark support."""
    return r"""
:root {
  --radius: 16px;
  --radius-lg: 22px;
  --gap: 16px;
  --shadow-sm: 0 4px 12px rgba(0,0,0,.08);
  --shadow-md: 0 8px 24px rgba(0,0,0,.12);
  --ring: 0 0 0 1px rgba(0,0,0,.06) inset;

  /* Light (default) */
  --bg: #f5f5f7;
  --surface: rgba(255,255,255,.72);
  --card: rgba(255,255,255,.82);
  --ink: #111111;
  --ink-2: #3a3a3c;
  --muted: #6e6e73;
  --link: #0071e3;      /* Apple blue */
  --accent: #34c759;    /* Apple green */
  --danger: #ff3b30;
  --divider: rgba(0,0,0,.08);
  --pill: #eef2ff;
  --code: #1f2937;
  --grad-1: #ffffff;
  --grad-2: #f6f7fb;
}

@media (prefers-color-scheme: dark) {
  :root {
    --bg: #0b0b0d;
    --surface: rgba(22,22,24,.6);
    --card: rgba(28,28,30,.7);
    --ink: #f5f5f7;
    --ink-2: #e5e5ea;
    --muted: #a1a1a6;
    --link: #0a84ff;
    --accent: #30d158;
    --danger: #ff453a;
    --divider: rgba(255,255,255,.08);
    --pill: #1c1c1e;
    --code: #e5e7eb;
    --grad-1: #0f1115;
    --grad-2: #0b0c10;
  }
}

* { box-sizing: border-box; }
html, body {
  margin: 0;
  padding: 0;
  font: 15px/1.45 -apple-system, BlinkMacSystemFont, "SF Pro Text", "SF Pro Display",
        "Helvetica Neue", Helvetica, Arial, ui-sans-serif, system-ui, "Segoe UI", sans-serif;
  color: var(--ink);
  background: linear-gradient(180deg, var(--grad-1), var(--grad-2));
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

a { color: var(--link); text-decoration: none; }
a:hover { text-decoration: underline; }

.container {
  max-width: 1100px;
  margin: 48px auto;
  padding: 0 20px;
}

.header {
  display: flex;
  align-items: center;
  gap: 14px;
  margin-bottom: 18px;
}

.logo {
  width: 40px; height: 40px;
  display: grid; place-items: center;
  border-radius: 12px;
  background: linear-gradient(180deg, rgba(52,199,89,.20), rgba(10,132,255,.18));
  backdrop-filter: saturate(180%) blur(20px);
  box-shadow: var(--ring);
}

.title {
  font-size: 22px;
  font-weight: 600;
  letter-spacing: .2px;
}

.subtitle {
  color: var(--muted);
  margin-top: 4px;
  font-size: 13px;
}

.hero {
  background: var(--surface);
  border-radius: var(--radius-lg);
  padding: 18px;
  box-shadow: var(--shadow-md);
  border: 1px solid var(--divider);
}

.hero-inner {
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 12px;
  flex-wrap: wrap;
}

.badges {
  display: flex;
  align-items: center;
  gap: 8px;
  flex-wrap: wrap;
}

.badge {
  display: inline-flex;
  align-items: center;
  gap: 6px;
  font-size: 12px;
  padding: 6px 10px;
  border-radius: 999px;
  background: var(--card);
  border: 1px solid var(--divider);
  box-shadow: var(--ring);
  color: var(--ink-2);
}

.badge .dot {
  width: 8px; height: 8px; border-radius: 50%;
  background: var(--accent);
  box-shadow: 0 0 0 3px rgba(52,199,89,.12);
}

.actions a.button,
.actions button.button {
  display: inline-flex; align-items: center; gap: 8px;
  padding: 10px 14px;
  background: var(--link);
  color: white;
  border-radius: 10px;
  border: none;
  box-shadow: 0 8px 18px rgba(10,132,255,.32);
  transition: transform .06s ease;
  cursor: pointer;
}
.actions a.button:hover,
.actions button.button:hover { text-decoration: none; transform: translateY(-1px); }
.actions .secondary {
  background: transparent; color: var(--link);
  border: 1px solid var(--divider);
  box-shadow: var(--ring);
}

.grid {
  display: grid;
  grid-template-columns: 1fr;
  gap: var(--gap);
  margin-top: 16px;
}
@media (min-width: 980px) {
  .grid { grid-template-columns: 1.2fr .8fr; }
}

.card {
  background: var(--card);
  border-radius: var(--radius);
  border: 1px solid var(--divider);
  box-shadow: var(--shadow-sm);
  padding: 14px;
}

.card h2 {
  margin: 4px 0 10px;
  font-size: 15px;
  font-weight: 600;
}

.list {
  margin: 0; padding: 0; list-style: none;
}
.list li {
  display: flex; align-items: center; justify-content: space-between;
  padding: 10px 8px;
  border-radius: 12px;
  border: 1px solid transparent;
}
.list li + li { margin-top: 6px; }
.list li:hover {
  background: var(--surface);
  border-color: var(--divider);
}

.kind {
  font-size: 12px;
  padding: 4px 8px;
  border-radius: 999px;
  background: var(--pill);
  color: var(--ink-2);
  border: 1px solid var(--divider);
}

.tools {
  display: flex;
  gap: 8px;
}
.tool-btn {
  font-size: 12px;
  padding: 6px 10px;
  border-radius: 9px;
  border: 1px solid var(--divider);
  background: var(--surface);
  color: var(--ink-2);
  box-shadow: var(--ring);
  cursor: pointer;
}
.tool-btn:hover { filter: brightness(1.05); }

.table-wrap {
  margin-top: 10px;
  display: none; /* toggled by JS */
}
.table-head {
  display: flex; gap: 8px; align-items: center; justify-content: space-between; margin-bottom: 8px;
}
.table-head .search {
  flex: 1;
  display: flex; gap: 8px; align-items: center;
  background: var(--surface);
  border: 1px solid var(--divider);
  border-radius: 10px;
  padding: 6px 10px;
}
.table-head input {
  flex: 1; border: 0; outline: 0; background: transparent; color: var(--ink);
  font-size: 14px;
}
.table {
  width: 100%;
  border-collapse: collapse;
  background: var(--card);
  border-radius: 12px;
  overflow: hidden;
  border: 1px solid var(--divider);
}
.table th, .table td {
  font-size: 13px;
  text-align: left;
  padding: 10px 12px;
  border-bottom: 1px solid var(--divider);
}
.table th {
  user-select: none;
  cursor: pointer;
}
.table tr:hover td {
  background: var(--surface);
}
.tag {
  font-size: 11px; padding: 2px 6px; border-radius: 999px; border: 1px solid var(--divider);
  background: var(--pill); color: var(--ink-2);
}

.footer {
  margin-top: 18px;
  color: var(--muted);
  font-size: 12.5px;
  text-align: center;
}
.code {
  color: var(--code);
  background: var(--pill);
  border: 1px solid var(--divider);
  padding: 2px 6px;
  border-radius: 8px;
}
    """


def _render_icon_corn_dna() -> str:
    """Tiny inline SVG mark (abstract corn + helix) to keep page asset-free."""
    return """
<svg width="22" height="22" viewBox="0 0 24 24" fill="none" aria-hidden="true">
  <defs>
    <linearGradient id="g" x1="0" y1="0" x2="1" y2="1">
      <stop offset="0%" stop-color="#34C759"/>
      <stop offset="100%" stop-color="#0A84FF"/>
    </linearGradient>
  </defs>
  <path d="M7 17c4-4 6-6 10-10" stroke="url(#g)" stroke-width="2" stroke-linecap="round"/>
  <path d="M7 7c2.2 2 3.8 3.5 6 6" stroke="url(#g)" stroke-width="2" stroke-linecap="round"/>
  <path d="M12 3c3 0 5 2 5 5 0 6-7 11-11 11-3 0-5-2-5-5 0-6 7-11 11-11z"
        stroke="url(#g)" stroke-opacity=".35" fill="none" stroke-width="1.2"/>
</svg>
"""


def _render_js(history_json: str, current_job_id: str) -> str:
    """Client-side interactions (toggle table, filter, sort, copy link)."""
    # history_json is already a JSON string; we inject as-is.
    return f"""
<script>
const PREV_RUNS = {history_json};
const CURRENT_JOB = {json.dumps(current_job_id)};

function qs(sel, el=document) {{ return el.querySelector(sel); }}
function qsa(sel, el=document) {{ return Array.from(el.querySelectorAll(sel)); }}

function toggleHistory() {{
  const wrap = qs('#history-wrap');
  wrap.style.display = (wrap.style.display === 'none' || !wrap.style.display) ? 'block' : 'none';
}}

function copyLink(href) {{
  navigator.clipboard.writeText(href).then(() => {{
    const el = qs('#copy-toast');
    el.textContent = 'Link copied';
    el.style.display = 'inline-block';
    setTimeout(() => el.style.display = 'none', 1200);
  }});
}}

function renderTable(rows) {{
  const tbody = qs('#history-tbody');
  tbody.innerHTML = '';
  rows.forEach(r => {{
    const tr = document.createElement('tr');
    tr.innerHTML = `
      <td><code class="tag">${{r.is_current ? 'current' : ''}}</code></td>
      <td><a href="${{r.report_url}}" target="_blank" rel="noopener">${{r.job_id}}</a></td>
      <td>${{r.created_at}}</td>
      <td style="text-align:right;">
        <div class="tools">
          <button class="tool-btn" onclick="copyLink('${{r.report_url}}')">Copy link</button>
          <a class="tool-btn" href="${{r.report_url}}" target="_blank" rel="noopener">Open</a>
        </div>
      </td>`;
    tbody.appendChild(tr);
  }});
}}

let sortKey = 'created_at';
let sortAsc = false;

function sortBy(key) {{
  if (sortKey === key) sortAsc = !sortAsc;
  else {{ sortKey = key; sortAsc = (key !== 'created_at'); }}
  applyFilterAndSort();
}}

function applyFilterAndSort() {{
  const q = qs('#search-input').value.trim().toLowerCase();
  let rows = PREV_RUNS.slice();
  if (q) {{
    rows = rows.filter(r =>
      r.job_id.toLowerCase().includes(q) ||
      r.report_url.toLowerCase().includes(q) ||
      (r.created_at || '').toLowerCase().includes(q)
    );
  }}
  rows.sort((a,b) => {{
    let va = a[sortKey] || '', vb = b[sortKey] || '';
    if (sortKey === 'created_at') {{
      // lexical works with ISO 8601
    }}
    if (va < vb) return sortAsc ? -1 : 1;
    if (va > vb) return sortAsc ? 1 : -1;
    return 0;
  }});
  renderTable(rows);
}}

window.addEventListener('DOMContentLoaded', () => {{
  // init
  renderTable(PREV_RUNS);
  qs('#search-input').addEventListener('input', applyFilterAndSort);
  // keyboard shortcut: press "h" to toggle history
  window.addEventListener('keydown', (e) => {{
    if (['INPUT','TEXTAREA'].includes((e.target.tagName||''))) return;
    if (e.key.toLowerCase() === 'h') {{
      toggleHistory();
    }}
  }});
}});
</script>
"""


def _render_html(
    job_id: str,
    artifacts: List[RunArtifact],
    reports_public_base: str,
    history_json: str,
) -> str:
    """Build the Apple-style HTML."""
    title = f"CornStructor results — {job_id}"
    items = "\n".join(
        f"""<li>
  <div style="display:flex;align-items:center;gap:10px;min-width:0;">
    <a href="{html.escape(a.rel_path)}" target="_blank" rel="noopener" style="min-width:0;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;">
      {html.escape(a.title)}
    </a>
  </div>
  <span class="kind">{html.escape(a.kind)}</span>
</li>"""
        for a in artifacts
    ) or "<li><em>No report artifacts were found for this run.</em></li>"

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{html.escape(title)}</title>
  <meta name="color-scheme" content="light dark" />
  <style>{_render_css()}</style>
</head>
<body>
  <div class="container">
    <div class="header">
      <div class="logo">{_render_icon_corn_dna()}</div>
      <div>
        <div class="title">CornStructor Results</div>
        <div class="subtitle">Run <span class="code">{html.escape(job_id)}</span> · served at <span class="code">{html.escape(reports_public_base)}/{html.escape(job_id)}/</span></div>
      </div>
    </div>

    <div class="hero">
      <div class="hero-inner">
        <div class="badges">
          <span class="badge"><span class="dot"></span> Completed</span>
          <span class="badge">Artifacts: {len(artifacts)}</span>
          <span class="badge">History: <span id="hist-count">{html.escape(str(len(json.loads(history_json))))}</span></span>
          <span id="copy-toast" class="badge" style="display:none;background:var(--accent);color:#fff;border-color:transparent;">Link copied</span>
        </div>
        <div class="actions">
          <a class="button" href="./" target="_self" rel="noopener">Refresh</a>
          <a class="button secondary" href="../" target="_self" rel="noopener">Run Folder</a>
          <button class="button" onclick="toggleHistory()">Open previous reports</button>
        </div>
      </div>
    </div>

    <div class="grid">
      <div class="card">
        <h2>Reports & Visualizations</h2>
        <ul class="list">
          {items}
        </ul>
        <div style="margin-top:10px;" class="tools">
          <button class="tool-btn" onclick="copyLink(window.location.href)">Copy this page link</button>
          <a class="tool-btn" href="{html.escape(reports_public_base)}/{html.escape(job_id)}/index.html" target="_blank" rel="noopener">Open in new tab</a>
        </div>
      </div>

      <div class="card">
        <h2>Previous reports</h2>
        <div class="table-wrap" id="history-wrap">
          <div class="table-head">
            <div class="search">
              <span>🔎</span>
              <input id="search-input" placeholder="Filter by job ID, date, or URL…" />
            </div>
            <div class="tools">
              <button class="tool-btn" onclick="applyFilterAndSort()">Reset sort</button>
            </div>
          </div>
          <table class="table" id="history-table">
            <thead>
              <tr>
                <th onclick="sortBy('is_current')">Status</th>
                <th onclick="sortBy('job_id')">Run ID</th>
                <th onclick="sortBy('created_at')">Created (UTC)</th>
                <th>Actions</th>
              </tr>
            </thead>
            <tbody id="history-tbody"></tbody>
          </table>
          <div class="subtitle" style="margin-top:8px;">Tip: press <span class="code">H</span> to toggle this panel.</div>
        </div>
        <div class="subtitle">Use the button above to open/close the history.</div>
      </div>
    </div>

    <div class="footer">
      CornStructor • Designed for clarity
    </div>
  </div>

  {_render_js(history_json, job_id)}
</body>
</html>"""


# ---------- Main API ----------

def write_run_index(outdir: Path, job_id: str, reports_public_base: str = "/reports") -> Optional[Path]:
    """
    Create/overwrite ``index.html`` in *outdir* with Apple-style UI + local history.

    Parameters
    ----------
    outdir:
        Filesystem directory for this run (e.g., ``/var/reports/{job_id}``).
    job_id:
        Unique identifier of the run, used in the page title and metadata.
    reports_public_base:
        Public mount under which *outdir* is served (default ``/reports``).

    Returns
    -------
    Optional[Path]
        The path to the generated ``index.html`` or ``None`` if *outdir* does not exist.
    """
    if not outdir.exists():
        return None

    # Discover artifacts for this run
    artifacts = _find_artifacts(outdir)

    # Build the canonical report URL for this run (index page itself)
    report_url = f"{reports_public_base.rstrip('/')}/{job_id}/index.html"

    # Persist into local history DB
    db_path = _db_path_for(outdir)
    try:
        conn = sqlite3.connect(db_path)
        try:
            _db_init(conn)
            _db_upsert_run(conn, job_id=job_id, report_url=report_url)
        finally:
            conn.close()
    except Exception:
        # If DB unavailable, we still render the page (history panel will be empty)
        pass

    # Load recent history to render inside the page (client-side table)
    history_json = _load_history_json(db_path, current_job_id=job_id)

    # Render and write HTML
    html_text = _render_html(
        job_id=job_id,
        artifacts=artifacts,
        reports_public_base=reports_public_base,
        history_json=history_json,
    )
    index_path = outdir / "index.html"
    index_path.write_text(html_text, encoding="utf-8")
    return index_path
