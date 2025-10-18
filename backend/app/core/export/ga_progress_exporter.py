# File: backend/app/core/export/ga_progress_exporter.py
# Version: v0.2.0
"""
Export per-generation GA optimization statistics for each cluster (parent fragment)
encountered during hierarchical assembly.

v0.2.0
- Support BOTH storage schemes:
  * legacy: node.ga_detail (list of {gen,best,mean,std}), optional node.ga_cluster_id
  * current: node.ga_log (list of best-per-gen floats) stored on the cluster/parent node
- Auto-generate a stable cluster_id if missing:
  "S-<root_id>_L-<level>_S-<start>_E-<end>"
- Deduplicate clusters by cluster_id.

Notes
- When only ga_log is available, we synthesize mean==best and std==0 so that
  downstream charts can still plot all series.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Set
import json

from backend.app.core.models.fragment_node import FragmentNode


def _collect_clusters(root: FragmentNode, root_id: str) -> List[Dict[str, Any]]:
    """
    Traverse the tree and collect unique clusters.
    A node is considered a 'cluster' if it has either:
      - ga_detail: list[{gen,best,mean,std}], or
      - ga_log: list[float] (best per generation).
    """
    clusters: List[Dict[str, Any]] = []
    seen: Set[str] = set()

    def consider(node: FragmentNode) -> None:
        # Prefer detailed stats if present
        detail = getattr(node, "ga_detail", None)
        log_list = getattr(node, "ga_log", None)

        if detail and isinstance(detail, list) and len(detail) > 0:
            gens = [int(round(float(d.get("gen", i + 1)))) for i, d in enumerate(detail)]
            best = [float(d.get("best", 0.0)) for d in detail]
            mean = [float(d.get("mean", 0.0)) for d in detail]
            std  = [float(d.get("std", 0.0)) for d in detail]
        elif isinstance(log_list, (list, tuple)) and len(log_list) > 0:
            gens = list(range(1, len(log_list) + 1))
            best = [float(x) for x in log_list]
            # Synthesize mean/std to keep visualizations happy
            mean = best[:]
            std  = [0.0] * len(best)
        else:
            return  # not a GA cluster node

        # Build/resolve cluster id
        cid = getattr(node, "ga_cluster_id", "") or ""
        try:
            s = int(getattr(node, "start", -1))
            e = int(getattr(node, "end", -1))
        except Exception:
            s, e = -1, -1
        if not cid:
            cid = f"S-{root_id}_L-{int(getattr(node, 'level', 0))}_S-{s}_E-{e}"

        if cid in seen:
            return

        # Parent metadata
        parent_level = int(getattr(node, "level", 0))

        clusters.append({
            "cluster_id": cid,
            "parent_level": parent_level,
            "parent_start": s,
            "parent_end": e,
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
    Compatible with both ga_detail and ga_log storage forms.
    """
    clusters = _collect_clusters(root, root_id)
    payload: Dict[str, Any] = {
        "root_id": root_id,
        "clusters": clusters,
    }
    out_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return payload
