# File: backend/app/core/export/json_exporter.py
# Version: v0.3.0

"""
Export a FragmentNode tree to a clean JSON file.
"""

import json
from pathlib import Path
from typing import Dict, Any

from backend.app.core.models.fragment_node import FragmentNode


def export_tree_to_json(root: FragmentNode,
                        json_path: Path,
                        root_id: str) -> None:
    """
    Export the full fragment tree to JSON for future reloading or analysis.
    """

    def serialize_node(node: FragmentNode) -> Dict[str, Any]:
        return {
            "fragment_id":  node.fragment_id,
            "level":        node.level,
            "start":        node.start,
            "end":          node.end,
            "seq":          node.seq,
            "strand":       node.strand,
            "overlap_prev": node.overlap_prev,
            "overlap_next": node.overlap_next,
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
