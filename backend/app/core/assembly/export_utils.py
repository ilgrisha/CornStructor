# File: backend/app/core/assembly/export_utils.py
# Version: v0.1.0

"""
Export utilities for HierarchicalAssembler.FragmentNode trees.

Functions:
  - export_fragments_to_fasta   : write all fragments to FASTA
  - export_tree_to_json         : nested JSON serialization
  - export_tree_to_html         : simple nested‐list HTML viewer
  - export_clusters_and_fragments_html : stub for combined cluster+fragment HTML
"""

import json
from pathlib import Path
from typing import List, Tuple, Dict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from backend.app.core.assembly.hierarchical_assembler import FragmentNode


def _traverse(node: FragmentNode) -> List[FragmentNode]:
    """Flatten the tree into a list (pre‐order)."""
    out = [node]
    for c in node.children:
        out.extend(_traverse(c))
    return out


def export_fragments_to_fasta(root: FragmentNode,
                              fasta_path: Path,
                              root_id: str) -> None:
    """
    Write every fragment (intermediate and leaf) as a FASTA record.
    ID format: <root_id>_lvl<level>_<start>_<end>
    """
    records: List[SeqRecord] = []
    for node in _traverse(root):
        rid = f"{root_id}_lvl{node.level}_{node.start}_{node.end}"
        records.append(SeqRecord(Seq(node.seq), id=rid, description=""))
    SeqIO.write(records, fasta_path, "fasta")


def export_tree_to_json(root: FragmentNode,
                        json_path: Path,
                        root_id: str) -> None:
    """
    Serialize the entire tree with overlaps into nested JSON.
    """
    def node_dict(n: FragmentNode) -> dict:
        return {
            "level":    n.level,
            "start":    n.start,
            "end":      n.end,
            "sequence": n.seq,
            "overlaps": [
                {"seq": ov[0], "start": ov[1], "end": ov[2]}
                for ov in n.overlaps
            ],
            "children": [node_dict(c) for c in n.children]
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
        header = (
            f"Level {n.level}: {n.start}–{n.end} "
            f"({len(n.seq)} bp)"
        )
        # overlaps list
        ov_items = "".join(
            f"<li>{ov[0]} [{ov[1]}–{ov[2]}]</li>"
            for ov in n.overlaps
        )
        # children
        child_html = "".join(render_node(c) for c in n.children)
        return (
            "<li>"
            f"<strong>{header}</strong>"
            + (f"<ul>{ov_items}</ul>" if ov_items else "")
            + (f"<ul>{child_html}</ul>" if child_html else "")
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
    Stub: combine your existing cluster‐visualizer with intermediate fragments.
    Currently just delegates to export_tree_to_html.
    """
    export_tree_to_html(root, html_path, root_id)
