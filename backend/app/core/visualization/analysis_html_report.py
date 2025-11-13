# File: backend/app/core/visualization/analysis_html_report.py
# Version: v0.1.0
"""
Render an HTML report for the construction-tree analysis JSON produced by
backend.app.core.export.analysis_exporter.analyze_tree_to_json.

Creates a single, self-contained HTML with:
- Header summary (root id, sequence length, global stats)
- Per-level summary table
- Per-overlap detailed table (sortable via simple JS)
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Any, List
import html
import json


def _fmt_float(x: float, nd: int = 3) -> str:
    try:
        return f"{float(x):.{nd}f}"
    except Exception:
        return str(x)


def export_analysis_html(analysis_json: Dict[str, Any] | Path, out_html: Path) -> None:
    """
    Accepts either an analysis dict or a path to the JSON file and writes an HTML report.
    """
    if isinstance(analysis_json, Path):
        data = json.loads(analysis_json.read_text(encoding="utf-8"))
    else:
        data = analysis_json

    root_id = html.escape(str(data.get("root_id", "unknown")))
    seq_len = int(data.get("sequence_length", 0))
    glb = data.get("global", {})
    levels = data.get("levels", {})
    overlaps: List[Dict[str, Any]] = list(data.get("overlaps", []))

    # Build per-level summary rows
    lvl_rows = []
    for lvl_str, summ in sorted(levels.items(), key=lambda kv: int(kv[0])):
        lvl_rows.append(
            """
            <tr>
              <td>{lvl}</td>
              <td>{cnt}</td>
              <td>{avg_nd}</td>
              <td>{min_nd}</td>
              <td>{avg_gc}</td>
              <td>{avg_tm}</td>
            </tr>
            """.format(
                lvl=html.escape(lvl_str),
                cnt=int(summ.get("count_overlaps", 0)),
                avg_nd=_fmt_float(summ.get("avg_norm_edit_distance", 0.0)),
                min_nd=_fmt_float(summ.get("min_norm_edit_distance", 0.0)),
                avg_gc=_fmt_float(summ.get("avg_gc_percent", 0.0), 2),
                avg_tm=_fmt_float(summ.get("avg_tm_celsius", 0.0), 2),
            )
        )

    # Build overlap detail rows
    ov_rows = []
    for ov in overlaps:
        ov_rows.append(
            """
            <tr>
              <td>{level}</td>
              <td>{parent}</td>
              <td>{left}</td>
              <td>{right}</td>
              <td>{start}</td>
              <td>{end}</td>
              <td>{length}</td>
              <td class="mono">{seq}</td>
              <td>{gc}</td>
              <td>{tm}</td>
              <td>{mind}</td>
            </tr>
            """.format(
                level=int(ov.get("level", -1)),
                parent=html.escape(str(ov.get("parent_id", ""))),
                left=html.escape(str(ov.get("left_id", ""))),
                right=html.escape(str(ov.get("right_id", ""))),
                start=int(ov.get("start", 0)),
                end=int(ov.get("end", 0)),
                length=int(ov.get("length", 0)),
                seq=html.escape(str(ov.get("seq", ""))),
                gc=_fmt_float(ov.get("gc_percent", 0.0), 2),
                tm=_fmt_float(ov.get("tm_celsius", 0.0), 2),
                mind=_fmt_float(ov.get("min_norm_edit_distance_to_others", 1.0), 3),
            )
        )

    html_doc = f"""
<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8" />
  <title>Assembly Analysis: {root_id}</title>
  <style>
    body {{ font-family: -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif; margin: 18px; }}
    h1, h2 {{ margin: 0.2em 0; }}
    table {{ border-collapse: collapse; width: 100%; margin: 12px 0; }}
    th, td {{ border: 1px solid #ccc; padding: 6px 8px; text-align: left; }}
    th {{ background: #f5f5f5; cursor: pointer; }}
    .mono {{ font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; font-size: 12px; }}
    .small {{ color: #555; }}
  </style>
  <script>
    // Simple column sort for tables
    function sortTable(id, col, numeric=false) {{
      const table = document.getElementById(id);
      const tbody = table.tBodies[0];
      const rows = Array.from(tbody.rows);
      const dir = table.getAttribute('data-sort-dir') === 'asc' ? 'desc' : 'asc';
      table.setAttribute('data-sort-dir', dir);
      rows.sort((a,b) => {{
        let va = a.cells[col].innerText;
        let vb = b.cells[col].innerText;
        if (numeric) {{ va = parseFloat(va) || 0; vb = parseFloat(vb) || 0; }}
        if (va < vb) return dir === 'asc' ? -1 : 1;
        if (va > vb) return dir === 'asc' ? 1 : -1;
        return 0;
      }});
      rows.forEach(r => tbody.appendChild(r));
    }}
  </script>
</head>
<body>
  <h1>Assembly Analysis: {root_id}</h1>
  <div class="small">Sequence length: {seq_len} bp</div>

  <h2>Global</h2>
  <table>
    <tr><th>Total overlaps</th><td>{int(glb.get('count_overlaps', 0))}</td></tr>
    <tr><th>Avg normalized edit distance</th><td>{_fmt_float(glb.get('avg_norm_edit_distance', 0.0))}</td></tr>
    <tr><th>Min normalized edit distance</th><td>{_fmt_float(glb.get('min_norm_edit_distance', 0.0))}</td></tr>
  </table>

  <h2>Per-level summary</h2>
  <table id="levels" data-sort-dir="asc">
    <thead>
      <tr>
        <th onclick="sortTable('levels',0,true)">Level</th>
        <th onclick="sortTable('levels',1,true)">Overlaps</th>
        <th onclick="sortTable('levels',2,true)">Avg NED</th>
        <th onclick="sortTable('levels',3,true)">Min NED</th>
        <th onclick="sortTable('levels',4,true)">Avg GC %</th>
        <th onclick="sortTable('levels',5,true)">Avg Tm °C</th>
      </tr>
    </thead>
    <tbody>
      {''.join(lvl_rows)}
    </tbody>
  </table>

  <h2>Overlaps</h2>
  <table id="overlaps" data-sort-dir="asc">
    <thead>
      <tr>
        <th onclick="sortTable('overlaps',0,true)">Level</th>
        <th onclick="sortTable('overlaps',1,false)">Parent</th>
        <th onclick="sortTable('overlaps',2,false)">Left</th>
        <th onclick="sortTable('overlaps',3,false)">Right</th>
        <th onclick="sortTable('overlaps',4,true)">Start</th>
        <th onclick="sortTable('overlaps',5,true)">End</th>
        <th onclick="sortTable('overlaps',6,true)">Length</th>
        <th>Sequence</th>
        <th onclick="sortTable('overlaps',8,true)">GC %</th>
        <th onclick="sortTable('overlaps',9,true)">Tm °C</th>
        <th onclick="sortTable('overlaps',10,true)">Min NED</th>
      </tr>
    </thead>
    <tbody>
      {''.join(ov_rows)}
    </tbody>
  </table>
</body>
</html>
"""
    out_html.write_text(html_doc, encoding="utf-8")
