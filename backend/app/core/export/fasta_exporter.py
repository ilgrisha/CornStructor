# File: backend/app/core/export/fasta_exporter.py
# Version: v0.2.0

"""
FASTA export utilities for FragmentNode trees.
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
        records.append(SeqRecord(Seq(n.seq), id=rid, description=""))
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
            records.append(SeqRecord(Seq(n.seq), id=rid, description=""))
    SeqIO.write(records, fasta_path, "fasta")
