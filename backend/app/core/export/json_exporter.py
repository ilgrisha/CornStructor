# File: backend/app/core/export/json_exporter.py
# Version: v0.3.0

"""
Export a FragmentNode tree to a clean JSON file.
"""

import json
from pathlib import Path
from typing import Dict, Any, Optional

from backend.app.core.models.fragment_node import FragmentNode
from backend.app.core.assembly.fitness_utils import compute_tm

_RC_MAP = str.maketrans("ACGTacgt", "TGCAtgca")


def _orient(seq: str, strand: str) -> str:
    if (strand or "+") == "-":
        return (seq or "").translate(_RC_MAP)[::-1]
    return seq or ""


def _tm_or_none(seq: str) -> Optional[float]:
    if not seq:
        return None
    try:
        return float(compute_tm(seq))
    except Exception:
        return None


def export_tree_to_json(root: FragmentNode,
                        json_path: Path,
                        root_id: str) -> None:
    """
    Export the full fragment tree to JSON for future reloading or analysis.
    """

    def serialize_node(node: FragmentNode) -> Dict[str, Any]:
        strand = node.strand or "+"
        seq_raw = node.seq or ""
        ovp_raw = node.overlap_prev or ""
        ovn_raw = node.overlap_next or ""
        ovp = _orient(ovp_raw, strand)
        ovn = _orient(ovn_raw, strand)
        ovp_tm = _tm_or_none(ovp)
        ovn_tm = _tm_or_none(ovn)

        return {
            "fragment_id":  node.fragment_id,
            "level":        node.level,
            "start":        node.start,
            "end":          node.end,
            "seq":          seq_raw,
            "strand":       strand,
            "overlap_prev": ovp_raw,
            "overlap_next": ovn_raw,
            "overlap_prev_tm": ovp_tm,
            "overlap_next_tm": ovn_tm,
            "is_oligo":     node.is_oligo,
            "ga_log":       node.ga_log,
            "children":     [serialize_node(child) for child in node.children]
        }

    # Top-level object includes root_id for reference
    tree_data = {
        "root_id": root_id,
        "tree": serialize_node(root)
    }

    json_path.write_text(json.dumps(tree_data, indent=2), encoding="utf-8")
