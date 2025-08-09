# File: backend/app/core/assembly/export_utils.py
# Version: v0.2.0

"""
Export utilities for HierarchicalAssembler.FragmentNode trees.

Functions:
  - export_fragments_to_fasta     : write all fragments to FASTA
  - export_tree_to_json           : nested JSON serialization
  - export_tree_to_html           : simple nested‐list HTML viewer
  - export_clusters_and_fragments_html : combined tree + per‐level cluster pages
"""

import json
from pathlib import Path
from typing import List, Dict, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from backend.app.core.assembly.hierarchical_assembler import FragmentNode
from backend.app.core.visualization.oligo_html_visualizer import (
    reverse_complement,
    render_oligo_list,
    visualize_cluster_html,
    format_progress_log,
)


def _traverse(node: FragmentNode) -> List[FragmentNode]:
    """Flatten the tree into a pre‐order list of nodes."""
    out = [node]
    for c in getattr(node, "children", []):
        out.extend(_traverse(c))
    return out


def export_fragments_to_fasta(root: FragmentNode,
                              fasta_path: Path,
                              root_id: str) -> None:
    """
    Write every fragment (intermediate + leaf) as a FASTA record.
    ID: <root_id>_lvl<level>_<start>_<end>
    """
    records: List[SeqRecord] = []
    for n in _traverse(root):
        rid = f"{root_id}_lvl{n.level}_{n.start}_{n.end}"
        records.append(SeqRecord(Seq(n.seq), id=rid, description=""))
    SeqIO.write(records, fasta_path, "fasta")


def export_tree_to_json(root: FragmentNode,
                        json_path: Path,
                        root_id: str) -> None:
    """
    Serialize the entire assembly tree (with overlaps) to nested JSON.
    """
    def node_dict(n: FragmentNode) -> Dict:
        return {
            "level":    n.level,
            "start":    n.start,
            "end":      n.end,
            "strand":   getattr(n, "strand", "+"),
            "sequence": n.seq,
            "overlap_prev": getattr(n, "overlap_prev", ""),
            "overlap_next": getattr(n, "overlap_next", ""),
            "ga_log":   getattr(n, "ga_log", []),
            "children": [node_dict(c) for c in getattr(n, "children", [])]
        }

    out = {root_id: node_dict(root)}
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)


def export_tree_to_html(root: FragmentNode,
                        html_path: Path,
                        root_id: str) -> None:
    """
    Render a simple nested‐list HTML of fragments and their overlaps.
    """
    def render_node(n: FragmentNode) -> str:
        header = f"Level {n.level}: {n.start}–{n.end} ({len(n.seq)} bp)"
        children_html = "".join(render_node(c) for c in getattr(n, "children", []))
        return (
            "<li>"
            f"<strong>{header}</strong>"
            + (f"<ul>{children_html}</ul>" if children_html else "")
            + "</li>"
        )

    html = [
        "<html><head><meta charset='UTF-8'><title>Assembly Tree</title></head><body>",
        f"<h1>Assembly Tree: {root_id}</h1>",
        "<ul>",
        render_node(root),
        "</ul>",
        "</body></html>"
    ]
    html_path.write_text("\n".join(html), encoding="utf-8")


def export_clusters_and_fragments_html(root: FragmentNode,
                                       html_path: Path,
                                       root_id: str) -> None:
    """
    Combined report:
      1) Top-level: nested‐list assembly tree.
      2) For each internal node: visualize its children as a 'cluster'
         (list + aligned strands + GA log).
    """
    html: List[str] = [
        "<html><head><meta charset='UTF-8'><title>Clusters & Fragments</title>",
        "<style>",
        ".mono{font-family:monospace;white-space:pre-wrap;}",
        ".oligo-list{margin:0 0 12px 16px;font-family:monospace;}",
        ".oligo-list li{margin-bottom:8px;}",
        "</style></head><body>",
        f"<h1>Assembly & Cluster Report: {root_id}</h1>",
        "<h2>Assembly Tree Overview</h2>"
    ]

    def render_tree(n: FragmentNode) -> str:
        hdr = f"{n.level}:{n.start}-{n.end} ({len(n.seq)} bp)"
        kids = "".join(render_tree(c) for c in getattr(n, "children", []))
        return f"<li>{hdr}" + (f"<ul>{kids}</ul>" if kids else "") + "</li>"

    html.append("<ul>")
    html.append(render_tree(root))
    html.append("</ul>")

    def traverse_and_render(n: FragmentNode):
        if not getattr(n, "children", []):
            return

        html.append(f"<h2>Level {n.level} Cluster at {n.start}–{n.end}</h2>")

        cluster: List[Tuple[str, str, int, int, str, str]] = []

        for child in n.children:
            prev_seq = getattr(child, "overlap_prev", "")
            next_seq = getattr(child, "overlap_next", "")
            strand   = getattr(child, "strand", "+")

            # Full fragment with overlaps
            full_seq = prev_seq + child.seq + next_seq
            full_seq = full_seq if strand == "+" else reverse_complement(full_seq)

            # Relative start/end for display
            rel_start = child.start - n.start
            rel_end   = child.end - n.start

            cluster.append((full_seq, strand, rel_start, rel_end, prev_seq, next_seq))

        html.append(render_oligo_list(root_id, n.level, cluster))
        html.append("<div style='background:#f8f8f8;padding:8px;margin-bottom:16px;border:1px solid #ccc;'>")
        html.append(visualize_cluster_html(root_id, n.level, cluster, n.seq))
        html.append("</div>")

        if getattr(n, "ga_log", []):
            html.append(format_progress_log(n.ga_log))

        for c in n.children:
            traverse_and_render(c)

    traverse_and_render(root)

    html.append("</body></html>")
    html_path.write_text("\n".join(html), encoding="utf-8")
