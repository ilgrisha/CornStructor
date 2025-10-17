# File: backend/app/core/export/fasta_exporter.py
# Version: v0.2.1

"""
FASTA export utilities for FragmentNode trees.

v0.2.1
- Honor node.strand: if '-', export reverse complement; if '+', export as-is.
"""

from pathlib import Path
from typing import List
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from backend.app.core.models.fragment_node import FragmentNode


def _traverse(node: FragmentNode) -> List[FragmentNode]:
    """Flatten the tree into a pre‐order list of nodes."""
    out = [node]
    for c in getattr(node, "children", []):
        out.extend(_traverse(c))
    return out


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


def export_fragments_to_fasta(root: FragmentNode,
                              fasta_path: Path,
                              root_id: str) -> None:
    """
    Write every fragment (intermediate + leaf) as a FASTA record.
    ID: <fragment_id>_<start>_<end>
    """
    records: List[SeqRecord] = []
    for n in _traverse(root):
        rid = f"{n.fragment_id}_{n.start}_{n.end}"
        seq_out = _oriented_seq(n)
        records.append(SeqRecord(Seq(seq_out), id=rid, description=""))
    SeqIO.write(records, fasta_path, "fasta")


def export_oligos_to_fasta(root: FragmentNode,
                           fasta_path: Path,
                           root_id: str) -> None:
    """
    Export only oligo (leaf) nodes as FASTA.
    ID: <fragment_id>_<start>_<end>
    """
    records: List[SeqRecord] = []
    for n in _traverse(root):
        if n.is_oligo:
            rid = f"{n.fragment_id}_{n.start}_{n.end}"
            seq_out = _oriented_seq(n)
            records.append(SeqRecord(Seq(seq_out), id=rid, description=""))
    SeqIO.write(records, fasta_path, "fasta")
