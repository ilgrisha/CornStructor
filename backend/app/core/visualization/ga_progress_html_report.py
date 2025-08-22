# File: backend/app/core/visualization/ga_progress_html_report.py
# Version: v0.3.0
"""
Render an HTML report for GA progress JSON produced by
backend.app.core.export.ga_progress_exporter.export_ga_progress_json.

Changes in v0.3.0:
- Per-graph y-axis scaling uses best, mean, and mean±std with padding for better visualization.
- Added a sparse, faint grid (horizontal + light verticals at ticks).

Changes in v0.2.0:
- Sort graphs by level and fragment location.
- Render standard deviation as a shaded band (mean ± std) instead of a line.
- Add generation tick labels on the x-axis.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List
import json
import html


def _sort_key(cl: Dict[str, Any]) -> tuple:
    def coerce(v, default=10**12):
        try:
            return int(v)
        except Exception:
            return default
    return (
        coerce(cl.get("parent_level"), 10**9),
        coerce(cl.get("parent_start"), 10**10),
        coerce(cl.get("parent_end"), 10**11),
        str(cl.get("cluster_id", "")),
    )


def export_ga_progress_html(progress_json: Dict[str, Any] | Path, out_html: Path) -> None:
    if isinstance(progress_json, Path):
        data = json.loads(progress_json.read_text(encoding="utf-8"))
    else:
        data = progress_json

    root_id = html.escape(str(data.get("root_id", "unknown")))
    clusters: List[Dict[str, Any]] = list(data.get("clusters", []) or [])

    # Sort clusters by level/start/end
    clusters.sort(key=_sort_key)

    blocks = []
    for idx, cl in enumerate(clusters):
        cid = html.escape(str(cl.get("cluster_id", f"cluster_{idx}")))
        gens = cl.get("gens", []) or []
        best = cl.get("best", []) or []
        mean = cl.get("mean", []) or []
        std = cl.get("std", []) or []
        title = f"L{cl.get('parent_level','?')} [{cl.get('parent_start','?')},{cl.get('parent_end','?')}]"
        canvas_id = f"c{idx}"
        # JS arrays as JSON literals
        def to_js(arr):
            return json.dumps([float(x) for x in arr])
        js_g = to_js(gens if gens else list(range(1, max(len(best), len(mean), len(std)) + 1)))
        js_b = to_js(best)
        js_m = to_js(mean)
        js_s = to_js(std)
        block = """
<section>
  <h3>%s</h3>
  <div class="small">%s</div>
  <canvas id="%s" width="700" height="240"></canvas>
  <script>
    (function(){
      const gens = %s;
      const best = %s;
      const mean = %s;
      const stdd = %s;
      const c = document.getElementById('%s');
      const ctx = c.getContext('2d');
      const W = c.width, H = c.height;
      const P = 35; const PW = W - P*2; const PH = H - P*2;
      const n = gens.length || 1;
      const xmin = 1, xmax = Math.max(1, gens[gens.length-1] || n);
      // compute band bounds
      const upper = mean.map((m,i)=> (m||0) + (stdd[i]||0));
      const lower = mean.map((m,i)=> (m||0) - (stdd[i]||0));
      function allVals(){
        let out = [];
        for (let i=0;i<arguments.length;i++){ if (Array.isArray(arguments[i])) out = out.concat(arguments[i]); }
        return out;
      }
      // y-range from best, mean, and band bounds (per-graph) with padding
      const vals = allVals(best, mean, upper, lower);
      let rawYmin = Math.min.apply(null, vals);
      let rawYmax = Math.max.apply(null, vals);
      if (!isFinite(rawYmin) || !isFinite(rawYmax)) { rawYmin = 0; rawYmax = 1; }
      const rawRange = rawYmax - rawYmin;
      let ymin, ymax;
      if (rawRange <= 1e-9) {
        const half = Math.max(0.5, Math.abs(rawYmax || 1) * 0.05);
        ymin = rawYmin - half;
        ymax = rawYmax + half;
      } else {
        const pad = rawRange * 0.08; // 8%% padding
        ymin = rawYmin - pad;
        ymax = rawYmax + pad;
      }
      function xMap(g){ return P + (g - xmin) / (xmax - xmin || 1) * PW; }
      function yMap(v){ return H - P - (v - ymin) / (ymax - ymin || 1) * PH; }
      // axes
      ctx.strokeStyle = '#999'; ctx.lineWidth = 1; ctx.beginPath();
      ctx.moveTo(P, P); ctx.lineTo(P, H-P); ctx.lineTo(W-P, H-P); ctx.stroke();
      ctx.fillStyle = '#333'; ctx.font = '12px sans-serif';
      ctx.fillText('Gen', W-P-28, H-P+26);
      ctx.save(); ctx.translate(P-42, P+6); ctx.rotate(-Math.PI/2); ctx.fillText('Fitness', 0, 0); ctx.restore();
      
      // draw x-axis ticks with generation numbers (limit to ~10 labels)
      const maxTicks = 10;
      const step = Math.max(1, Math.ceil(n / maxTicks));
      
      // sparse, faint grid (behind data)
      ctx.strokeStyle = 'rgba(0,0,0,0.08)';
      ctx.lineWidth = 1;
      // vertical grid lines at tick positions
      for (let i=0; i<n; i+=step){
        const g = gens[i] || (i+1);
        const x = xMap(g);
        ctx.beginPath();
        ctx.moveTo(x, P);
        ctx.lineTo(x, H-P);
        ctx.stroke();
      }
      // horizontal grid lines (5 divisions)
      const gridHLines = 5;
      for (let k=0; k<=gridHLines; k++){
        const yv = ymin + (ymax - ymin) * (k / gridHLines);
        const y = yMap(yv);
        ctx.beginPath();
        ctx.moveTo(P, y);
        ctx.lineTo(W-P, y);
        ctx.stroke();
      }
      
      // tick marks and labels
      ctx.fillStyle = '#444';
      ctx.strokeStyle = '#bbb';
      for (let i=0; i<n; i+=step){
        const g = gens[i] || (i+1);
        const x = xMap(g);
        ctx.beginPath();
        ctx.moveTo(x, H-P);
        ctx.lineTo(x, H-P+4);
        ctx.stroke();
        ctx.fillText(String(g), x-6, H-P+16);
      }

      function drawLine(arr, color){
        if (!arr || arr.length === 0) return;
        ctx.strokeStyle = color; ctx.lineWidth = 2; ctx.beginPath();
        for (let i=0;i<gens.length;i++){ const x=xMap(gens[i]|| (i+1)); const y=yMap(arr[i]||arr[arr.length-1]||0); if(i===0)ctx.moveTo(x,y); else ctx.lineTo(x,y);} 
        ctx.stroke();
      }

      // band: mean ± std
      if (mean.length && stdd.length){
        ctx.beginPath();
        for (let i=0;i<gens.length;i++){
          const x = xMap(gens[i] || (i+1));
          const y = yMap(lower[i] || lower[lower.length-1] || 0);
          if (i===0) ctx.moveTo(x,y); else ctx.lineTo(x,y);
        }
        for (let i=gens.length-1;i>=0;i--){
          const x = xMap(gens[i] || (i+1));
          const y = yMap(upper[i] || upper[upper.length-1] || 0);
          ctx.lineTo(x,y);
        }
        ctx.closePath();
        ctx.fillStyle = 'rgba(231, 111, 81, 0.20)'; // translucent band
        ctx.fill();
        // outline the band lightly
        ctx.strokeStyle = 'rgba(231, 111, 81, 0.6)';
        ctx.lineWidth = 1;
        ctx.stroke();
      }

      // lines
      drawLine(best, '#2a9d8f');
      drawLine(mean, '#264653');

      // legend
      const labels = [
        ['Best','#2a9d8f','line'],
        ['Avg','#264653','line'],
        ['Avg ± Std','rgba(231, 111, 81, 0.20)','band']
      ];
      let lx = P+10, ly = P+10;
      labels.forEach(([name,color,type])=>{
        if (type==='band'){
          ctx.fillStyle = color; ctx.fillRect(lx, ly-6, 18, 12);
          ctx.strokeStyle = 'rgba(231, 111, 81, 0.6)'; ctx.strokeRect(lx, ly-6, 18, 12);
        } else {
          ctx.strokeStyle = color; ctx.lineWidth = 3; ctx.beginPath(); ctx.moveTo(lx, ly); ctx.lineTo(lx+18, ly); ctx.stroke();
        }
        ctx.fillStyle = '#222'; ctx.fillText(' '+name, lx+22, ly+4);
        lx += 120;
      });
    })();
  </script>
</section>
""" % (title, cid, canvas_id, js_g, js_b, js_m, js_s, canvas_id)
        blocks.append(block)

    html_doc = f"""
<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8" />
  <title>GA Progress: {root_id}</title>
  <style>
    body {{ font-family: -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif; margin: 18px; }}
    h1, h2 {{ margin: 0.2em 0; }}
    section {{ border: 1px solid #ddd; padding: 10px 12px; margin: 12px 0; border-radius: 6px; }}
    .small {{ color: #555; font-size: 12px; }}
  </style>
</head>
<body>
  <h1>GA Progress: {root_id}</h1>
  <div class="small">Clusters: {len(clusters)}</div>
  {''.join(blocks) if blocks else '<p>No GA progress data found.</p>'}
</body>
</html>
"""
    out_html.write_text(html_doc, encoding="utf-8")
