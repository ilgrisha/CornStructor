# File: backend/app/core/export/csv_exporter.py
# Version: v0.2.0
"""
CSV exporters for fragments and oligos with strand-aware sequences and overlap Tm.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

from backend.app.core.assembly.fitness_utils import compute_tm

try:
    from backend.app.core.models.fragment_node import FragmentNode  # type: ignore
except Exception:
    FragmentNode = object  # type: ignore


_RC_MAP = str.maketrans("ACGTacgt", "TGCAtgca")


def _rc(seq: str) -> str:
    return (seq or "").translate(_RC_MAP)[::-1]


def _oriented(seq: str, strand: str) -> str:
    return _rc(seq) if (strand or "+") == "-" else seq


def _walk(root: FragmentNode) -> Iterable[FragmentNode]:
    stack = [root]
    while stack:
        n = stack.pop()
        yield n
        if getattr(n, "children", None):
            for c in reversed(n.children):
                stack.append(c)


def _row(
    n: FragmentNode,
    root_id: str,
    *,
    tm_method: str,
    tm_params: dict,
) -> Dict[str, str]:
    strand = str(getattr(n, "strand", "+") or "+")
    seq_raw = str(getattr(n, "seq", "") or "")
    ovp_raw = str(getattr(n, "overlap_prev", "") or "")
    ovn_raw = str(getattr(n, "overlap_next", "") or "")

    seq = _oriented(seq_raw, strand)
    ovp = _oriented(ovp_raw, strand) if ovp_raw else ""
    ovn = _oriented(ovn_raw, strand) if ovn_raw else ""

    try:
        ovp_tm = f"{compute_tm(ovp, method=tm_method, **tm_params):.3f}" if ovp else ""
    except Exception:
        ovp_tm = ""
    try:
        ovn_tm = f"{compute_tm(ovn, method=tm_method, **tm_params):.3f}" if ovn else ""
    except Exception:
        ovn_tm = ""

    s = int(getattr(n, "start", 0) or 0)
    e = int(getattr(n, "end", 0) or 0)
    L = max(0, e - s)

    return {
        "root_id": str(root_id),
        "fragment_id": str(getattr(n, "fragment_id", "")),
        "level": str(getattr(n, "level", "")),
        "is_oligo": "1" if bool(getattr(n, "is_oligo", False)) else "0",
        "strand": strand,
        "start": str(s),
        "end": str(e),
        "length": str(L),
        "overlap_prev": ovp,
        "overlap_prev_tm": ovp_tm,
        "overlap_next": ovn,
        "overlap_next_tm": ovn_tm,
        "sequence": seq,
    }


def export_csvs(
    root: FragmentNode,
    outdir: Path,
    *,
    root_id: str,
    tm_method: str,
    tm_params: dict,
) -> Tuple[Path, Path]:
    """Write fragments.csv (all nodes) and oligos.csv (leaf nodes) with oriented sequences + Tm columns."""
    all_rows: List[Dict[str, str]] = []
    oligo_rows: List[Dict[str, str]] = []

    for n in _walk(root):
        row = _row(n, root_id, tm_method=tm_method, tm_params=tm_params)
        all_rows.append(row)
        if bool(getattr(n, "is_oligo", False)):
            oligo_rows.append(row)

    frag_csv = outdir / "fragments.csv"
    oligo_csv = outdir / "oligos.csv"

    if all_rows:
        headers = list(all_rows[0].keys())
        with frag_csv.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=headers)
            w.writeheader()
            w.writerows(all_rows)

    if oligo_rows:
        headers = list(oligo_rows[0].keys())
        with oligo_csv.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=headers)
            w.writeheader()
            w.writerows(oligo_rows)

    return frag_csv, oligo_csv
