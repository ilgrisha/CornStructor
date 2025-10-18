# File: backend/app/core/export/fasta_oriented_exporter.py
# Version: v0.1.0
"""
Strand-aware FASTA exporters.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, List

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


def export_fragments_fasta_oriented(root: FragmentNode, outpath: Path, *, root_id: str) -> Path:
    """FASTA of all nodes (tree walk), strand-aware."""
    lines: List[str] = []
    for n in _walk(root):
        fid = str(getattr(n, "fragment_id", ""))
        lvl = str(getattr(n, "level", ""))
        strand = str(getattr(n, "strand", "+") or "+")
        s = int(getattr(n, "start", 0) or 0)
        e = int(getattr(n, "end", 0) or 0)
        seq = _oriented(str(getattr(n, "seq", "") or ""), strand)
        header = (
            f">{fid} | root={root_id} | level={lvl} | strand={strand} | "
            f"start={s} | end={e} | len={max(0, e - s)}"
        )
        lines.append(header)
        for i in range(0, len(seq), 80):
            lines.append(seq[i : i + 80])
    outpath.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
    return outpath


def export_oligos_fasta(root: FragmentNode, outpath: Path, *, root_id: str) -> Path:
    """FASTA with only leaf/oligo sequences, strand-aware."""
    lines: List[str] = []
    for n in _walk(root):
        if not bool(getattr(n, "is_oligo", False)):
            continue
        fid = str(getattr(n, "fragment_id", ""))
        lvl = str(getattr(n, "level", ""))
        strand = str(getattr(n, "strand", "+") or "+")
        s = int(getattr(n, "start", 0) or 0)
        e = int(getattr(n, "end", 0) or 0)
        seq = _oriented(str(getattr(n, "seq", "") or ""), strand)
        header = (
            f">{fid} | root={root_id} | level={lvl} | strand={strand} | "
            f"start={s} | end={e} | len={max(0, e - s)}"
        )
        lines.append(header)
        for i in range(0, len(seq), 80):
            lines.append(seq[i : i + 80])
    outpath.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
    return outpath
