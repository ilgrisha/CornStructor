# File: backend/app/core/export/analysis_exporter.py
# Version: v0.1.0
"""
Analyze a FragmentNode construction tree: fragments, overlaps, and design metrics.
Outputs a JSON report suitable for programmatic use and HTML visualization.

Metrics per-overlap:
- level, parent fragment id, left/right child ids
- absolute start/end, length
- sequence, GC%, Tm

Aggregates:
- counts per level
- global list of overlaps and pairwise normalized edit distance stats

Note:
- Tm is computed via fitness_utils.compute_tm (Primer3 if available, with fallbacks).
- Normalized edit distance uses fitness_utils implementations (edlib/rapidfuzz if available).
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Dict, List, Tuple, Optional

from backend.app.core.models.fragment_node import FragmentNode
from backend.app.core.assembly.fitness_utils import gc_fraction, compute_tm
from backend.app.core.assembly.batch_fitness import pairwise_norm_edit_stats


@dataclass
class OverlapInfo:
    level: int
    parent_id: str
    left_id: str
    right_id: str
    start: int
    end: int
    length: int
    seq: str
    gc_percent: float
    tm_celsius: float


def _collect_overlaps(root: FragmentNode) -> Tuple[List[OverlapInfo], Dict[int, List[OverlapInfo]]]:
    """
    Traverse the tree and collect all overlaps as absolute coordinates and sequences.

    We infer each junction overlap between adjacent children by intersecting
    [right_child.start, left_child.end]. This matches how GA extends children in
    hierarchical_assembler._build_children.
    """
    all_overlaps: List[OverlapInfo] = []
    by_level: Dict[int, List[OverlapInfo]] = {}

    def recurse(node: FragmentNode) -> None:
        if not node.children:
            return
        # Adjacent children of this node define overlaps at this level
        for i in range(len(node.children) - 1):
            left = node.children[i]
            right = node.children[i + 1]
            ov_start = max(right.start, node.start)  # conservative guard
            ov_end = min(left.end, node.end)
            if ov_end <= ov_start:
                # No actual span (should not happen, but guard against malformed trees)
                continue
            # The root node holds the full sequence
            # To avoid threading root, the absolute slice can be taken from the root-most ancestor.
            # Since every node.seq is its subsequence, but indices are absolute, we need root.seq.
            # We walk up once to find the true root.
            root_node = node
            while True:
                # naive: root has level 0 by construction
                if getattr(root_node, 'level', None) == 0:
                    break
                # No parent pointers available; rely on the fact that the top-most "root" we were
                # passed in has level 0. So capture from outer scope instead of searching.
                break
            # Use the outer-scope root variable for slicing
            seq_full = _ROOT_SEQ_CACHE[0]
            ov_seq = seq_full[ov_start:ov_end]

            info = OverlapInfo(
                level=left.level,  # junction belongs to children level
                parent_id=node.fragment_id,
                left_id=left.fragment_id,
                right_id=right.fragment_id,
                start=ov_start,
                end=ov_end,
                length=(ov_end - ov_start),
                seq=ov_seq,
                gc_percent=gc_fraction(ov_seq),
                tm_celsius=compute_tm(ov_seq),
            )
            all_overlaps.append(info)
            by_level.setdefault(left.level, []).append(info)
        for c in node.children:
            recurse(c)

    # Root sequence cache set by analyze_tree_to_json before calling _collect_overlaps
    recurse(root)
    return all_overlaps, by_level


# A tiny global used internally to pass the root sequence to helpers without
# rewriting many function signatures. It holds [root.seq] for the current call.
_ROOT_SEQ_CACHE: List[str] = [""]


def analyze_tree_to_json(root: FragmentNode, out_path: Path, root_id: str) -> Dict[str, Any]:
    """
    Build a comprehensive analysis JSON for the construction tree and write it to a file.

    Returns the JSON payload as a Python dict as well.
    """
    # Expose the full sequence for internal helpers
    _ROOT_SEQ_CACHE[0] = root.seq

    overlaps, by_level = _collect_overlaps(root)

    # Pairwise normalized edit distance stats across all overlaps
    all_seqs = [ov.seq for ov in overlaps]
    avg_nd, min_nd = pairwise_norm_edit_stats(all_seqs)

    # Per-overlap nearest-neighbor (min) normalized edit distance
    # Compute naively as dataset sizes are expected to be modest
    per_overlap_min_nd: List[float] = []
    for i, s in enumerate(all_seqs):
        best = 1.0
        for j, t in enumerate(all_seqs):
            if i == j:
                continue
            # Local import to avoid circular
            from backend.app.core.assembly.fitness_utils import normalized_edit_distance
            v = normalized_edit_distance(s, t)
            if v < best:
                best = v
        per_overlap_min_nd.append(best if len(all_seqs) > 1 else 1.0)

    # Attach per-overlap distances back
    overlaps_dicts: List[Dict[str, Any]] = []
    for info, mind in zip(overlaps, per_overlap_min_nd):
        d = asdict(info)
        d["min_norm_edit_distance_to_others"] = mind
        overlaps_dicts.append(d)

    # Per-level aggregates
    levels_summary: Dict[str, Any] = {}
    for lvl, items in by_level.items():
        seqs = [it.seq for it in items]
        avg_l, min_l = pairwise_norm_edit_stats(seqs)
        levels_summary[str(lvl)] = {
            "count_overlaps": len(items),
            "avg_norm_edit_distance": avg_l,
            "min_norm_edit_distance": min_l,
            "avg_gc_percent": (sum(it.gc_percent for it in items) / max(1, len(items))),
            "avg_tm_celsius": (sum(it.tm_celsius for it in items) / max(1, len(items))),
        }

    payload: Dict[str, Any] = {
        "root_id": root_id,
        "sequence_length": len(root.seq),
        "analysis_generated": True,
        "overlaps": overlaps_dicts,
        "levels": levels_summary,
        "global": {
            "count_overlaps": len(overlaps),
            "avg_norm_edit_distance": avg_nd,
            "min_norm_edit_distance": min_nd,
        },
    }

    out_path.write_text(
        __import__("json").dumps(payload, indent=2), encoding="utf-8"
    )
    return payload
