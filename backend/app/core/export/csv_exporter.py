# File: backend/app/core/export/csv_exporter.py
# Version: v0.1.1

"""
Export leaf fragment (oligo) sequences to CSV.

v0.1.1
- Honor node.strand: if '-', export reverse complement; if '+', export as-is.

Fields:
- sequence_id
- sequence
- sequence_length

Behavior:
- Traverses the FragmentNode tree and collects *leaf* nodes. A node is
  considered a leaf if `is_oligo` is True OR it has no children.
- Writes a CSV with the three fields above in left-to-right order.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Tuple
import csv
from Bio.Seq import Seq

from backend.app.core.models.fragment_node import FragmentNode


def _oriented_seq(n: FragmentNode) -> str:
    """
    Return the sequence oriented for export:
    - '+' strand: raw sequence
    - '-' strand: reverse complement
    """
    raw = getattr(n, "seq", "") or ""
    if getattr(n, "strand", "+") == "-":
        return str(Seq(raw).reverse_complement())
    return raw


def _collect_leaf_sequences(node: FragmentNode) -> List[Tuple[str, str, int]]:
    out: List[Tuple[str, str, int]] = []
    children = getattr(node, "children", None) or []
    if getattr(node, "is_oligo", False) or len(children) == 0:
        seq_out = _oriented_seq(node)
        sid = getattr(node, "fragment_id", "") or ""
        out.append((sid, seq_out, len(seq_out)))
        return out
    for ch in children:
        out.extend(_collect_leaf_sequences(ch))
    return out


def export_fragments_to_csv(root: FragmentNode, out_path: Path, root_id: str) -> Path:
    """
    Write all leaf sequences to CSV at `out_path`.
    Returns the written path.
    """
    rows = _collect_leaf_sequences(root)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Write CSV
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sequence_id", "sequence", "sequence_length"])
        for sid, seq, L in rows:
            w.writerow([sid, seq, L])

    return out_path
