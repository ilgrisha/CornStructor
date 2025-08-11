# File: backend/app/core/visualization/cluster_html_report.py
# Version: v0.6.0

"""
Generates level-specific HTML visualizations for FragmentNode trees.

- Produces one HTML file per internal tree level (excluding leaves).
- Uses oligo_html_visualizer to render aligned strand diagrams per cluster.
- Includes GA progress logs once per level (if available).
- Creates an index.html linking to all level reports.
- Can be run as a standalone script with a JSON input path.
"""

from pathlib import Path
from typing import List, Tuple, Optional
import argparse

from backend.app.core.models.fragment_node import FragmentNode
from backend.app.core.models.construction_tree_loader import load_tree_from_json
from backend.app.core.visualization.oligo_html_visualizer import (
    render_oligo_list,
    visualize_cluster_html,
    format_progress_log,
    reverse_complement,   # unused directly here but kept if you log/debug
    MAX_VISIBLE_NT
)

# Type alias:
# (fragment_id, cluster, parent_seq, ga_log, abs_start)
ClusterBundle = Tuple[str, List[Tuple[str, str, int, int, str, str]], str, List[float], int]


def collect_clusters_by_level(root: FragmentNode) -> dict[int, List[ClusterBundle]]:
    """
    Traverse tree and group clusters by level.

    Returns:
        Dict[level, List of tuples:
            (
              fragment_id: str,
              cluster: List[(seq_full, strand, start_rel, end_rel, prev_ov, next_ov)],
              parent_seq: str,
              ga_log: List[float],
              abs_start: int          # absolute start of the parent node (0-based)
            )
        ]
    """
    levels: dict[int, List[ClusterBundle]] = {}

    def recurse(node: FragmentNode):
        children = getattr(node, "children", [])
        if not children:
            return

        parent_start = node.start  # absolute 0-based start of this (parent) fragment

        cluster = []
        for child in children:
            seq = child.seq
            strand = child.strand
            start_rel = child.start - parent_start
            end_rel   = child.end   - parent_start
            prev_ov = child.overlap_prev
            next_ov = child.overlap_next
            seq_full = seq if strand == "+" else reverse_complement(seq)
            cluster.append((seq_full, strand, start_rel, end_rel, prev_ov, next_ov))

        levels.setdefault(node.level + 1, []).append(
            (node.fragment_id, cluster, node.seq, node.ga_log, parent_start)
        )

        for c in children:
            recurse(c)

    recurse(root)
    return levels


def _cluster_bracket_abs(cluster: List[Tuple[str, str, int, int, str, str]], abs_start: int) -> str:
    """
    Compute [left – right] bracket in 1-based ABSOLUTE coordinates for a cluster block.
    """
    frag_left_rel  = min(c[2] for c in cluster)
    frag_right_rel = max(c[3] for c in cluster)
    left_disp  = abs_start + frag_left_rel + 1       # 1-based inclusive
    right_disp = abs_start + frag_right_rel          # 1-based inclusive
    return f"[{left_disp} – {right_disp}]"


def export_level_html(
    level: int,
    clusters: List[ClusterBundle],
    output_path: Path
) -> None:
    """
    Export all clusters at a given level to a single HTML file.

    Args:
        level: Level index (1-based from root’s children).
        clusters: List of (fragment_id, cluster, full_seq, ga_log, abs_start).
        output_path: Path to write the HTML file.
    """
    html = [
        "<html><head><meta charset='UTF-8'><title>Oligo Clusters at Level {}</title>".format(level),
        "<style>",
        ".mono{font-family:monospace;white-space:pre;display:block;}",
        ".oligo-list{margin:0 0 12px 16px;font-family:monospace;}",
        ".oligo-list li{margin-bottom:8px;}",
        f".scroll-box{{overflow-x:auto;white-space:nowrap;padding:8px;background:#f9f9f9;"
        f"border:1px solid #ccc;margin-bottom:16px;max-width:{MAX_VISIBLE_NT}ch;}}",
        "body{font-family:sans-serif;padding:20px;background:#fcfcfc;color:#333;}",
        "h1{font-size:28px;border-bottom:2px solid #ccc;padding-bottom:8px;}",
        "h2{margin-top:30px;color:#2a4d69;}",
        "h3{margin-top:20px;color:#3e5c76;}",
        "</style></head><body>",
        f"<h1>🧬 Oligo Clusters at Level {level}</h1>"
    ]

    # show GA once per level: pick the first non-empty log we find
    level_ga_log: Optional[List[float]] = None

    for i, (fragment_id, cluster, parent_seq, ga_log, abs_start) in enumerate(clusters, start=1):
        bracket = _cluster_bracket_abs(cluster, abs_start)
        html.append(f"<h2>Cluster {i} from fragment {fragment_id} {bracket}</h2>")
        # Oligo list with ABS badges (1-based)
        html.append(render_oligo_list(fragment_id, i, cluster, abs_start=abs_start))
        # Alignment + ruler (abs_start drives 1-based ruler display)
        html.append(visualize_cluster_html(fragment_id, i, cluster, parent_seq, abs_start=abs_start))

        if (not level_ga_log) and ga_log:
            level_ga_log = ga_log

    if level_ga_log:
        html.append("<hr>")
        html.append("<h2>GA Progress (representative for this level)</h2>")
        html.append(format_progress_log(level_ga_log))

    html.append("</body></html>")
    output_path.write_text("\n".join(html), encoding="utf-8")


def export_index_html(levels: dict[int, List[ClusterBundle]], out_path: Path) -> None:
    """
    Export an index.html linking to all per-level HTML reports.

    Args:
        levels: Mapping of level → cluster list.
        out_path: Output directory to write index.html.
    """
    html = [
        "<html><head><meta charset='UTF-8'><title>Cluster Index</title>",
        "<style>",
        "body{font-family:sans-serif;padding:20px;background:#fcfcfc;color:#333;}",
        "h1{font-size:28px;border-bottom:2px solid #ccc;padding-bottom:8px;}",
        "ul{padding-left:24px;}",
        "li{margin-bottom:8px;}",
        "a{color:#2a4d69;text-decoration:none;}",
        "a:hover{text-decoration:underline;}",
        "</style></head><body>",
        "<h1>📚 Assembly Cluster Reports</h1>",
        "<ul>"
    ]

    for level in sorted(levels.keys()):
        fname = f"level_{level}.html"
        count = len(levels[level])
        html.append(f"<li><a href='{fname}'>Level {level}</a> &mdash; {count} fragments</li>")

    html.append("</ul></body></html>")
    index_path = out_path / "index.html"
    index_path.write_text("\n".join(html), encoding="utf-8")


def export_all_levels(root: FragmentNode, out_dir: Path) -> None:
    """
    Export separate HTML files for each tree level (excluding leaves).

    Args:
        root: Root FragmentNode of the tree.
        out_dir: Output directory to write HTML files.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    level_clusters = collect_clusters_by_level(root)

    for level, clusters in sorted(level_clusters.items()):
        out_path = out_dir / f"level_{level}.html"
        export_level_html(level, clusters, out_path)

    export_index_html(level_clusters, out_dir)


def main():
    """
    Standalone CLI entry point.
    Usage:
        python cluster_html_report.py --json data/tree_fragments.json --out data/output
    """
    parser = argparse.ArgumentParser(description="Export per-level HTML cluster reports from a FragmentNode tree.")
    parser.add_argument("--json", required=True, type=Path, help="Path to input JSON file.")
    parser.add_argument("--out",  required=True, type=Path, help="Directory to save HTML files.")
    args = parser.parse_args()

    print(f"Loading tree from: {args.json}")
    root = load_tree_from_json(args.json)

    print(f"Exporting HTML cluster reports to: {args.out}")
    export_all_levels(root, args.out)

    print("Done.")


if __name__ == "__main__":
    main()
