# File: backend/app/core/models/construction_tree_loader.py
# Version: v0.2.0

"""
Load a saved FragmentNode tree from a JSON file for debugging or reuse.
"""

import json
from pathlib import Path
from typing import Any, Dict

from backend.app.core.models.fragment_node import FragmentNode


def load_fragment_node_from_dict(data: Dict[str, Any]) -> FragmentNode:
    """
    Recursively reconstruct a FragmentNode from a dictionary.
    """
    return FragmentNode(
        fragment_id=data.get("fragment_id", ""),
        level=data["level"],
        start=data["start"],
        end=data["end"],
        seq=data["seq"],
        strand=data["strand"],
        overlap_prev=data.get("overlap_prev", ""),
        overlap_next=data.get("overlap_next", ""),
        is_oligo=data.get("is_oligo", False),
        ga_log=data.get("ga_log", []),
        children=[
            load_fragment_node_from_dict(child) for child in data.get("children", [])
        ]
    )


def load_tree_from_json(json_path: Path) -> FragmentNode:
    """
    Load the FragmentNode tree from a JSON file previously saved with export_tree_to_json.

    Returns:
        root_node: the root FragmentNode of the reconstructed tree
    """
    with json_path.open("r", encoding="utf-8") as f:
        data = json.load(f)
        root_data = data.get("tree")
        return load_fragment_node_from_dict(root_data)
