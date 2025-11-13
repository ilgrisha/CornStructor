# File: backend/app/core/visualization/cluster_html_report.py
# Version: v0.7.3
"""
Generates level-specific HTML visualizations for FragmentNode trees.

v0.7.3
- Removed obsolete reference to `enforce_span_across_junction`.
- Uses new config module locations: backend.app.config.config_level/global.
- Safer constraints panel rendering with consistent HTML escaping.
- Minor style polish and clearer section headers.

Summary
-------
- Produces one HTML file per internal tree level (excluding leaves).
- Uses `oligo_html_visualizer` to render aligned strand diagrams per cluster.
- Includes GA progress logs once per level (if available).
- Creates an `index.html` linking to all level reports.
- Optionally decorates each level page with constraints from `levels.json`
  and GA knobs / Tm method from `globals.json`.

Standalone usage:
    python -m backend.app.core.visualization.cluster_html_report \
        --json <tree.json> --out <dir> [--levels <levels.json>] [--globals <globals.json>]
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Tuple, Optional, Dict
import argparse
import html

from backend.app.core.models.fragment_node import FragmentNode
from backend.app.core.models.construction_tree_loader import load_tree_from_json
from backend.app.core.visualization.oligo_html_visualizer import (
    render_oligo_list,
    visualize_cluster_html,
    format_progress_log,
    reverse_complement,  # kept handy for debugging
    MAX_VISIBLE_NT,
)

# Optional config imports (used if CLI passes --levels / --globals)
from backend.app.config.config_level import LevelConfig, load_levels_config
from backend.app.config.config_global import GlobalConfig, load_global_config


# Type alias for per-cluster payload:
# (fragment_id, cluster, parent_seq, ga_log, abs_start)
ClusterBundle = Tuple[
    str,
    List[Tuple[str, str, int, int, str, str]],
    str,
    List[float],
    int,
]


def collect_clusters_by_level(root: FragmentNode) -> Dict[int, List[ClusterBundle]]:
    """
    Traverse the tree and group clusters by level.

    A "cluster" is the set of child fragments under a given parent node.
    We package each cluster along with its parent sequence, GA progress
    log (if recorded), and the absolute start coordinate of the parent
    to support 1-based absolute display.

    Returns:
        Dict[level, List[ClusterBundle]]
            Where each ClusterBundle is:
            (
              fragment_id: str,
              cluster: List[(seq_full, strand, start_rel, end_rel, prev_ov, next_ov)],
              parent_seq: str,
              ga_log: List[float],
              abs_start: int  # absolute start of the parent node (0-based)
            )
    """
    levels: Dict[int, List[ClusterBundle]] = {}

    def recurse(node: FragmentNode) -> None:
        children = getattr(node, "children", [])
        if not children:
            return

        parent_start = node.start  # absolute 0-based start of this (parent) fragment

        cluster: List[Tuple[str, str, int, int, str, str]] = []
        for child in children:
            seq = child.seq
            strand = child.strand
            start_rel = child.start - parent_start
            end_rel = child.end - parent_start
            prev_ov = child.overlap_prev
            next_ov = child.overlap_next
            seq_full = seq if strand == "+" else reverse_complement(seq)
            cluster.append((seq_full, strand, start_rel, end_rel, prev_ov, next_ov))

        # Child level is parent.level + 1
        levels.setdefault(node.level + 1, []).append(
            (node.fragment_id, cluster, node.seq, node.ga_log, parent_start)
        )

        for c in children:
            recurse(c)

    recurse(root)
    return levels


def _cluster_bracket_abs(
    cluster: List[Tuple[str, str, int, int, str, str]],
    abs_start: int,
) -> str:
    """
    Compute '[left – right]' bracket in 1-based ABSOLUTE coordinates
    for a cluster block.

    Args:
        cluster: List of tuples for children in the cluster.
        abs_start: Absolute 0-based start coordinate of the parent fragment.

    Returns:
        A display string like "[123 – 456]".
    """
    frag_left_rel = min(c[2] for c in cluster)
    frag_right_rel = max(c[3] for c in cluster)
    left_disp = abs_start + frag_left_rel + 1  # 1-based inclusive
    right_disp = abs_start + frag_right_rel  # 1-based inclusive
    return f"[{left_disp} – {right_disp}]"


def _render_constraints_panel(
    level_cfg: Optional[LevelConfig],
    global_cfg: Optional[GlobalConfig],
) -> str:
    """
    Render a small constraints panel if configs are supplied.

    Args:
        level_cfg: Constraints for the current level (may be None).
        global_cfg: Global GA/Tm parameters (may be None).

    Returns:
        An HTML string with a styled panel, or an empty string if nothing to show.
    """
    if level_cfg is None and global_cfg is None:
        return ""

    rows: List[str] = []

    def fmt(x: object) -> str:
        """Format value with HTML escaping; show '—' for None."""
        return "—" if x is None else html.escape(str(x))

    if level_cfg is not None:
        disallowed = (
            ", ".join(level_cfg.overlap_disallowed_motifs)
            if level_cfg.overlap_disallowed_motifs
            else "—"
        )
        allowed = (
            ", ".join(f"{k}:{v}" for k, v in level_cfg.overlap_allowed_motifs.items())
            if level_cfg.overlap_allowed_motifs
            else "—"
        )

        rows.extend(
            [
                f"<tr><th>Overlap length</th><td>{level_cfg.overlap_min_size}–{level_cfg.overlap_max_size} bp</td></tr>",
                f"<tr><th>Tm (°C)</th><td>{fmt(level_cfg.overlap_tm_min)} – {fmt(level_cfg.overlap_tm_max)}</td></tr>",
                f"<tr><th>GC fraction</th><td>{fmt(level_cfg.overlap_gc_min)} – {fmt(level_cfg.overlap_gc_max)}</td></tr>",
                f"<tr><th>Max same-base run</th><td>{fmt(level_cfg.overlap_run_min)} – {fmt(level_cfg.overlap_run_max)}</td></tr>",
                f"<tr><th>Max gap allowed</th><td>{html.escape(str(level_cfg.max_gap_allowed))}</td></tr>",
                f"<tr><th>Disallowed motifs</th><td>{html.escape(disallowed)}</td></tr>",
                f"<tr><th>Allowed motifs (w)</th><td>{html.escape(allowed)}</td></tr>",
            ]
        )

    if global_cfg is not None:
        ga = global_cfg.ga_params()  # attribute-access view
        rows.extend(
            [
                f"<tr><th>Tm method</th><td>{html.escape(global_cfg.tm_method)}</td></tr>",
                f"<tr><th>CPU workers fraction</th><td>{global_cfg.cpu_workers_fraction:.2f}</td></tr>",
                f"<tr><th>GA pop × gens</th><td>{ga.population_size} × {ga.num_generations}</td></tr>",
                f"<tr><th>Mutation / Crossover</th><td>{ga.mutation_rate_initial:.2f} / {ga.crossover_rate:.2f}</td></tr>",
                f"<tr><th>Elitism / Inject</th><td>{ga.elitism_count} / {ga.random_injection_rate:.2f}</td></tr>",
                f"<tr><th>Tournament size</th><td>{ga.tournament_size}</td></tr>",
            ]
        )

    return (
        "<div style='margin:10px 0 18px 0;padding:12px;background:#fff;border:1px solid #ddd;'>"
        "<h3 style='margin:0 0 8px 0;color:#3e5c76;'>Constraints</h3>"
        "<table class='constraints' style='border-collapse:collapse;font-size:14px'>"
        "<tbody>"
        + "\n".join(rows)
        + "</tbody></table></div>"
    )


def export_level_html(
    level: int,
    clusters: List[ClusterBundle],
    output_path: Path,
    level_cfg: Optional[LevelConfig] = None,
    global_cfg: Optional[GlobalConfig] = None,
) -> None:
    """
    Export all clusters at a given level to a single HTML file.

    Args:
        level: Level index (1-based from root’s children).
        clusters: List of (fragment_id, cluster, full_seq, ga_log, abs_start).
        output_path: Path to write the HTML file.
        level_cfg: Optional LevelConfig to display constraints used at this level.
        global_cfg: Optional GlobalConfig to display GA/assembly knobs.
    """
    html_lines = [
        f"<html><head><meta charset='UTF-8'><title>Oligo Clusters at Level {level}</title>",
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
        "table{border:1px solid #ddd;background:#fff;}",
        "th,td{font-weight:normal;padding:4px 8px;}",
        "th{color:#555;text-align:left;border-right:1px solid #eee;}",
        "</style></head><body>",
        f"<h1>🧬 Oligo Clusters at Level {level}</h1>",
    ]

    # Constraints panel (if configs provided)
    html_lines.append(_render_constraints_panel(level_cfg, global_cfg))

    # Show one representative GA progress log per level (first non-empty found)
    level_ga_log: Optional[List[float]] = None

    for i, (fragment_id, cluster, parent_seq, ga_log, abs_start) in enumerate(
        clusters, start=1
    ):
        bracket = _cluster_bracket_abs(cluster, abs_start)
        html_lines.append(
            f"<h2>Cluster {i} from fragment {html.escape(fragment_id)} {bracket}</h2>"
        )
        # Oligo list with ABS badges (1-based). `render_oligo_list` supports abs_start kw.
        html_lines.append(
            render_oligo_list(fragment_id, i, cluster, abs_start=abs_start)
        )
        # Alignment + ruler (abs_start drives 1-based ruler display)
        html_lines.append(
            visualize_cluster_html(
                fragment_id, i, cluster, parent_seq, abs_start=abs_start
            )
        )

        if level_ga_log is None and ga_log:
            level_ga_log = ga_log

    if level_ga_log:
        html_lines.append("<hr>")
        html_lines.append("<h2>GA Progress (representative for this level)</h2>")
        html_lines.append(format_progress_log(level_ga_log))

    html_lines.append("</body></html>")
    output_path.write_text("\n".join(html_lines), encoding="utf-8")


def export_index_html(levels: Dict[int, List[ClusterBundle]], out_path: Path) -> None:
    """
    Export an `index.html` linking to all per-level HTML reports.

    Args:
        levels: Mapping of level → cluster list.
        out_path: Output directory to write `index.html`.
    """
    html_lines = [
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
        "<ul>",
    ]

    for level in sorted(levels.keys()):
        fname = f"level_{level}.html"
        count = len(levels[level])
        html_lines.append(
            f"<li><a href='{fname}'>Level {level}</a> &mdash; {count} fragments</li>"
        )

    html_lines.append("</ul></body></html>")
    index_path = out_path / "index.html"
    index_path.write_text("\n".join(html_lines), encoding="utf-8")


def export_all_levels(
    root: FragmentNode,
    out_dir: Path,
    levels_cfg: Optional[Dict[int, LevelConfig]] = None,
    global_cfg: Optional[GlobalConfig] = None,
) -> None:
    """
    Export separate HTML files for each tree level (excluding leaves).

    Args:
        root: Root `FragmentNode` of the tree.
        out_dir: Output directory to write HTML files.
        levels_cfg: Optional dict[level → LevelConfig] to annotate constraints per level.
        global_cfg: Optional `GlobalConfig` to show GA/assembly parameters.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    level_clusters = collect_clusters_by_level(root)

    for level, clusters in sorted(level_clusters.items()):
        out_path = out_dir / f"level_{level}.html"
        lvl_cfg = levels_cfg.get(level) if levels_cfg else None
        export_level_html(level, clusters, out_path, lvl_cfg, global_cfg)

    export_index_html(level_clusters, out_dir)


def main() -> None:
    """
    Standalone CLI entry point.

    Usage:
        python -m backend.app.core.visualization.cluster_html_report \
            --json data/tree_fragments.json \
            --out data/cluster_reports \
            [--levels app/config/levels.json] \
            [--globals app/config/globals.json]
    """
    parser = argparse.ArgumentParser(
        description="Export per-level HTML cluster reports from a FragmentNode tree."
    )
    parser.add_argument(
        "--json",
        required=True,
        type=Path,
        help="Path to input JSON tree (export_tree_to_json output).",
    )
    parser.add_argument(
        "--out", required=True, type=Path, help="Directory to save HTML files."
    )
    parser.add_argument(
        "--levels",
        required=False,
        type=Path,
        help="Optional levels.json to annotate constraints.",
    )
    parser.add_argument(
        "--globals",
        required=False,
        type=Path,
        help="Optional globals.json to annotate GA knobs.",
    )
    args = parser.parse_args()

    print(f"Loading tree from: {args.json}")
    root = load_tree_from_json(args.json)

    levels_cfg: Optional[Dict[int, LevelConfig]] = None
    global_cfg: Optional[GlobalConfig] = None

    if args.levels:
        try:
            levels_cfg = load_levels_config(args.levels)
            print(f"Loaded levels config from: {args.levels}")
        except Exception as e:
            print(f"Warning: failed to load levels config {args.levels}: {e}")

    if args.globals:
        try:
            global_cfg = load_global_config(args.globals)
            print(f"Loaded globals config from: {args.globals}")
        except Exception as e:
            print(f"Warning: failed to load globals config {args.globals}: {e}")

    print(f"Exporting HTML cluster reports to: {args.out}")
    export_all_levels(root, args.out, levels_cfg, global_cfg)
    print("Done.")


if __name__ == "__main__":
    main()
