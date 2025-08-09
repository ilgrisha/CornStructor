# File: backend/app/core/visualization/tree_html_exporter.py
# Version: v0.1.0

"""
HTML nested tree exporter for FragmentNode hierarchy.
"""

from pathlib import Path
from typing import List

from backend.app.core.models.fragment_node import FragmentNode


def export_tree_to_html(root: FragmentNode,
                        html_path: Path,
                        root_id: str) -> None:
    """
    Render a simple nested‐list HTML of fragments and their overlaps.
    """
    def render_node(n: FragmentNode) -> str:
        header = f"Level {n.level}: {n.start}–{n.end} ({len(n.seq)} bp)"
        ov_items = "".join(
            f"<li>{ov[0]} [{ov[1]}–{ov[2]}]</li>"
            for ov in getattr(n, "overlaps", [])
        )
        children_html = "".join(render_node(c) for c in getattr(n, "children", []))
        return (
            "<li>"
            f"<strong>{header}</strong>"
            + (f"<ul>{ov_items}</ul>" if ov_items else "")
            + (f"<ul>{children_html}</ul>" if children_html else "")
            + "</li>"
        )

    html: List[str] = [
        "<html><head><meta charset='UTF-8'><title>Assembly Tree</title></head><body>",
        f"<h1>Assembly Tree: {root_id}</h1>",
        "<ul>",
        render_node(root),
        "</ul>",
        "</body></html>"
    ]
    html_path.write_text("\n".join(html), encoding="utf-8")
