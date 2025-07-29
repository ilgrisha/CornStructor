# File: backend/app/core/visualization/oligo_html_visualizer.py
# Version: v0.6.2

"""
HTML visualization of full-length oligo clusters.
Shows, for each cluster:
- the exact original fragment (5′→3′)
- a list of oligo IDs, strands, & full sequences
- the top strand (sense oligos) & bottom (antisense rev-comp)
- gaps as dashes
- GA progress log
"""

from pathlib import Path
from typing import List, Tuple, Dict


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT","TGCA"))[::-1]


def color_nucleotide(nt: str) -> str:
    cmap = {"A":"#d73027","C":"#4575b4","G":"#1a9850","T":"#984ea3","-":"#888888"}
    return f"<span style='color:{cmap.get(nt,'black')};'>{nt}</span>"


def color_sequence(seq: str) -> str:
    return ''.join(color_nucleotide(n) for n in seq)


def render_oligo_list(
    sid: str,
    ci: int,
    cluster: List[Tuple[str,str,int,int]]
) -> str:
    items = []
    for oi,(seq,strand,s,e) in enumerate(cluster,1):
        oid = f"{sid}_c{ci}_o{oi}"
        lbl = "sense" if strand=="+" else "antisense"
        items.append(
            f"<li><strong>{oid}</strong> ({lbl}): "
            f"<span class='mono'>{color_sequence(seq)}</span></li>"
        )
    return "<ul class='oligo-list'>" + "".join(items) + "</ul>"


def visualize_cluster_html(
    sid: str,
    ci: int,
    cluster: List[Tuple[str,str,int,int]],
    full_seq: str
) -> str:
    # define fragment
    start = cluster[0][2]
    end   = cluster[-1][3]
    frag = full_seq[start:end]
    L = len(frag)

    # original fragment
    orig = f"5′→ {color_sequence(frag)} 3′"

    # init blank
    top = ['-']*L
    bot = ['-']*L

    # place each oligo
    for seq,strand,s,e in cluster:
        rel_s = s - start
        if strand=='+':
            for i,nt in enumerate(seq):
                top[rel_s+i] = nt
        else:
            rc = reverse_complement(seq)
            for i,nt in enumerate(rc):
                bot[rel_s+i] = nt

    top_line = "5′→ "+color_sequence("".join(top))+" 3′"
    bot_line = "3′← "+color_sequence("".join(bot))+"  5′"

    return (
        f"<div class='mono'>{orig}</div>"
        f"<div class='mono'>{top_line}</div>"
        f"<div class='mono'>{bot_line}</div>"
    )


def format_progress_log(log: List[float]) -> str:
    if not log: return ""
    items = "".join(f"<li>Gen{i+1}: {f:.1f}</li>" for i,f in enumerate(log))
    return "<h4>GA Progress</h4><ul>"+items+"</ul>"


def export_html_report(
    output_path: Path,
    clusters_by_sequence: Dict[str,List[List[Tuple[str,str,int,int]]]],
    ga_logs: Dict[str,List[float]] = None,
    full_sequences: Dict[str,str] = None
) -> None:
    html = [
        "<html><head><meta charset='UTF-8'><title>Oligo Design Visualization</title>",
        "<style>",
        ".mono{font-family:monospace;white-space:pre-wrap;}",
        ".oligo-list{margin:0 0 8px 16px;font-family:monospace;}",
        ".oligo-list li{margin-bottom:4px;}",
        "</style></head><body>",
        "<h1>Oligo Design Visualization</h1>"
    ]

    for sid, clusters in clusters_by_sequence.items():
        html.append(f"<h2>Sequence: {sid}</h2>")
        seq = full_sequences.get(sid,"")
        for ci,cluster in enumerate(clusters,1):
            html.append(f"<h3>Cluster {ci}</h3>")
            html.append(render_oligo_list(sid,ci,cluster))
            html.append("<div style='background:#f8f8f8;padding:8px;margin-bottom:12px;border:1px solid #ccc;'>")
            html.append(visualize_cluster_html(sid,ci,cluster,seq))
            html.append("</div>")
        if ga_logs and sid in ga_logs:
            html.append(format_progress_log(ga_logs[sid]))

    html.append("</body></html>")
    output_path.write_text("\n".join(html), encoding='utf-8')
