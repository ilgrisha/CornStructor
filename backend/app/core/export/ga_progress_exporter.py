# File: backend/app/core/export/ga_progress_exporter.py
# Version: v0.1.0
"""
Export per-generation GA optimization statistics for each cluster (parent fragment)
encountered during hierarchical assembly.

The GA (GAOverlapSelector) now records per-generation stats in `progress_detail`
with entries {gen, best, mean, std}. The hierarchical assembler attaches this as
`ga_detail` to each child along with a stable `ga_cluster_id` constructed from the
parent region.

This module traverses the final FragmentNode tree, deduplicates by ga_cluster_id,
aggregates per-cluster series, and writes a JSON file suitable for visualization.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Tuple, Set
import json

from backend.app.core.models.fragment_node import FragmentNode


def _collect_clusters(root: FragmentNode) -> List[Dict[str, Any]]:
    """
    Traverse the tree and collect unique clusters using node.ga_cluster_id.
    We store one entry per parent cluster. Since assembler attaches the same
    ga_detail to all siblings, we pick the first occurrence per cluster id.
    """
    clusters: List[Dict[str, Any]] = []
    seen: Set[str] = set()

    def consider(node: FragmentNode) -> None:
        cid = getattr(node, "ga_cluster_id", "") or ""
        detail = getattr(node, "ga_detail", None)
        if cid and detail and isinstance(detail, list) and cid not in seen:
            # Attempt to parse level/start/end from cluster id pattern:
            # "S-<root>_L-<level>_S-<start>_E-<end>"
            level = None
            start = None
            end = None
            try:
                parts = cid.split("_")
                for p in parts:
                    if p.startswith("L-"):
                        level = int(p[2:])
                    elif p.startswith("S-") and p[2:].isdigit():
                        # Beware: first S- is sequence id prefix; the later S-<start>
                        start = int(p[2:])
                    elif p.startswith("E-"):
                        end = int(p[2:])
            except Exception:
                pass

            gens = [int(round(float(d.get("gen", i+1)))) for i, d in enumerate(detail)]
            best = [float(d.get("best", 0.0)) for d in detail]
            mean = [float(d.get("mean", 0.0)) for d in detail]
            std = [float(d.get("std", 0.0)) for d in detail]

            clusters.append({
                "cluster_id": cid,
                "parent_level": level if level is not None else max(0, node.level - 1),
                "parent_start": start,
                "parent_end": end,
                "gens": gens,
                "best": best,
                "mean": mean,
                "std": std,
            })
            seen.add(cid)

    def walk(n: FragmentNode) -> None:
        consider(n)
        for c in getattr(n, "children", []) or []:
            walk(c)

    walk(root)
    return clusters


def export_ga_progress_json(root: FragmentNode, out_path: Path, root_id: str) -> Dict[str, Any]:
    """
    Build GA progress payload and write it as JSON. Returns the payload.
    """
    clusters = _collect_clusters(root)
    payload: Dict[str, Any] = {
        "root_id": root_id,
        "clusters": clusters,
    }
    out_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return payload
