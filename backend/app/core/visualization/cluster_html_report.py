# File: backend/app/core/visualization/cluster_html_report.py
# Version: v0.5.0

"""
Generates level-specific HTML visualizations for FragmentNode trees.

- Produces one HTML file per internal tree level (excluding leaves).
- Uses oligo_html_visualizer to render aligned strand diagrams per cluster.
- Includes GA progress logs if available.
- Creates an index.html linking to all level reports.
- Can be run as a standalone script with a JSON input path.
"""

from pathlib import Path
from typing import List, Tuple
import argparse

from backend.app.core.models.fragment_node import FragmentNode
from backend.app.core.models.construction_tree_loader import load_tree_from_json
from backend.app.core.visualization.oligo_html_visualizer import (
    render_oligo_list,
    visualize_cluster_html,
    format_progress_log,
    reverse_complement,
    MAX_VISIBLE_NT
)


def collect_clusters_by_level(root: FragmentNode) -> dict[int, List[Tuple[str, List[Tuple[str, str, int, int, str, str]], str, List[float]]]]:
    """
    Traverse tree and group clusters by level.

    Returns:
        Dict[level, List of tuples: (fragment_id, cluster, parent_seq, ga_log)]
    """
    levels: dict[int, List] = {}

    def recurse(node: FragmentNode):
        children = getattr(node, "children", [])
        if not children:
            return

        parent_start = node.start
        cluster = []
        for child in children:
            seq = child.seq
            strand = child.strand
            start = child.start - parent_start
            end = child.end - parent_start
            prev_ov = child.overlap_prev
            next_ov = child.overlap_next
            seq_full = seq if strand == "+" else reverse_complement(seq)
            cluster.append((seq_full, strand, start, end, prev_ov, next_ov))

        levels.setdefault(node.level + 1, []).append(
            (node.fragment_id, cluster, node.seq, node.ga_log)
        )

        for c in children:
            recurse(c)

    recurse(root)
    return levels


def export_level_html(
    level: int,
    clusters: List[Tuple[str, List[Tuple[str, str, int, int, str, str]], str, List[float]]],
    output_path: Path
) -> None:
    """
    Export all clusters at a given level to a single HTML file.

    Args:
        level: Level index (1-based).
        clusters: List of tuples per cluster: (fragment_id, cluster, full_seq, ga_log).
        output_path: Path to write the HTML file.
    """
    html = [
        "<html><head><meta charset='UTF-8'><title>Oligo Clusters at Level {}</title>".format(level),
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
        f"<h1>🧬 Oligo Clusters at Level {level}</h1>"
    ]

    for i, (fragment_id, cluster, parent_seq, ga_log) in enumerate(clusters, start=1):
        html.append(f"<h2>Cluster {i} from fragment {fragment_id}</h2>")
        html.append(render_oligo_list(fragment_id, i, cluster))
        html.append("<div class='scroll-box'>")
        html.append(visualize_cluster_html(fragment_id, i, cluster, parent_seq))
        html.append("</div>")
        if ga_log:
            html.append(format_progress_log(ga_log))

    html.append("</body></html>")
    output_path.write_text("\n".join(html), encoding="utf-8")


def export_index_html(levels: dict[int, List], out_path: Path) -> None:
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
    parser.add_argument("--out", required=True, type=Path, help="Directory to save HTML files.")
    args = parser.parse_args()

    print(f"Loading tree from: {args.json}")
    root = load_tree_from_json(args.json)

    print(f"Exporting HTML cluster reports to: {args.out}")
    export_all_levels(root, args.out)

    print("Done.")


if __name__ == "__main__":
    main()
