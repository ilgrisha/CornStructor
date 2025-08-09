# File: backend/app/core/visualization/oligo_html_visualizer.py
# Version: v0.6.8

"""
HTML visualization of full‐length oligo clusters.

Features:
  - Horizontal scrolling of sequences and alignments.
  - Per-cluster fragment visualization: original sequence, top/bottom strands.
  - Colored nucleotides.
  - Oligo list with IDs and orientation.
  - GA progress log (optional).
  - Responsive layout with maximum visible nucleotide width before scrolling.

Fixes:
  - Prevents sequence wrapping by using white-space: pre and overflow-x:auto.
  - Separates each aligned strand onto its own line.
  - Adds a distinguishable visual header.
"""

from pathlib import Path
from typing import List, Tuple, Dict

MAX_VISIBLE_NT = 150  # Controls visible nucleotide width before scrolling (approximate)


def reverse_complement(seq: str) -> str:
    complement = str.maketrans("ACGT", "TGCA")
    return seq.translate(complement)[::-1]


def color_nucleotide(nt: str) -> str:
    cmap = {
        "A": "#d73027",  # red
        "C": "#4575b4",  # blue
        "G": "#1a9850",  # green
        "T": "#984ea3",  # purple
        "-": "#888888"   # grey for gaps
    }
    return f"<span style='color:{cmap.get(nt, 'black')}'>{nt}</span>"


def color_sequence(seq: str) -> str:
    return "".join(color_nucleotide(n) for n in seq)


def render_oligo_list(
    sid: str,
    ci: int,
    cluster: List[Tuple[str, str, int, int, str, str]]
) -> str:
    items = []
    for oi, (seq, strand, start, end, *_overlaps) in enumerate(cluster, start=1):
        oid = f"{sid}_c{ci}_o{oi}"
        lbl = "sense" if strand == "+" else "antisense"
        items.append(
            "<li>"
            f"<strong>{oid}</strong> ({lbl}): "
            f"<span class='mono'>{color_sequence(seq)}</span>"
            "</li>"
        )
    return (
        "<div class='scroll-box'>"
        "<ul class='oligo-list mono'>"
        + "".join(items) +
        "</ul></div>"
    )


def visualize_cluster_html(
    sid: str,
    ci: int,
    cluster: List[Tuple[str, str, int, int, str, str]],
    full_seq: str
) -> str:
    fragment_start = cluster[0][2]
    fragment_end = cluster[-1][3]
    fragment_seq = full_seq[fragment_start:fragment_end]
    L = len(fragment_seq)

    # Clip view to MAX_VISIBLE_NT if needed
    orig_line = f"5′→ {color_sequence(fragment_seq)} 3′"

    top = ["-"] * L
    bot = ["-"] * L

    for seq, strand, start, end, *_ in cluster:
        rel = start - fragment_start
        if strand == "+":
            for i, nt in enumerate(seq):
                if 0 <= rel + i < L:
                    top[rel + i] = nt
        else:
            rc = reverse_complement(seq)
            for i, nt in enumerate(rc):
                if 0 <= rel + i < L:
                    bot[rel + i] = nt

    top_line = "5′→ " + color_sequence("".join(top)) + " 3′"
    bot_line = "3′← " + color_sequence("".join(bot)) + "  5′"

    return (
        "<div class='scroll-box'>"
            "<div class='mono'>"
                f"{orig_line}\n"
                f"{top_line}\n"
                f"{bot_line}"
            "</div>"
        "</div>"
        )


def format_progress_log(log: List[float]) -> str:
    if not log:
        return ""
    items = "".join(f"<li>Gen {i+1}: {f:.1f}</li>" for i, f in enumerate(log))
    return "<h4>GA Progress</h4><ul>" + items + "</ul>"


def export_html_report(
    output_path: Path,
    clusters_by_sequence: Dict[str, List[List[Tuple[str, str, int, int, str, str]]]],
    ga_logs: Dict[str, List[float]] = None,
    full_sequences: Dict[str, str] = None
) -> None:
    html = [
        "<html><head><meta charset='UTF-8'><title>Oligo Design Visualization</title>",
        "<style>",
        ".mono{font-family:monospace;white-space:pre;display:inline-block;}",
        ".oligo-list{margin:0 0 12px 16px;font-family:monospace;}",
        ".oligo-list li{margin-bottom:8px;}",
        f".scroll-box{{overflow-x:auto;white-space:nowrap;padding:8px;background:#f9f9f9;"
        f"border:1px solid #ccc;margin-bottom:16px;max-width:{MAX_VISIBLE_NT}ch;}}",
        "body{font-family:sans-serif;padding:20px;background:#fcfcfc;color:#333;}",
        "h1{font-size:28px;border-bottom:2px solid #ccc;padding-bottom:8px;}",
        "h2{margin-top:30px;color:#2a4d69;}",
        "h3{margin-top:20px;color:#3e5c76;}",
        "</style></head><body>",
        "<h1>🧬 Oligo Design Visualization</h1>"
    ]

    for sid, clusters in clusters_by_sequence.items():
        html.append(f"<h2>Sequence: {sid}</h2>")
        full_seq = full_sequences.get(sid, "")
        for ci, cluster in enumerate(clusters, start=1):
            html.append(f"<h3>Cluster {ci}</h3>")
            html.append(render_oligo_list(sid, ci, cluster))
            html.append(visualize_cluster_html(sid, ci, cluster, full_seq))
        if ga_logs and sid in ga_logs:
            html.append(format_progress_log(ga_logs[sid]))

    html.append("</body></html>")
    output_path.write_text("\n".join(html), encoding="utf-8")
