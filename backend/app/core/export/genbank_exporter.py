# File: backend/app/core/export/genbank_exporter.py
# Version: v0.4.1
"""
GenBank exporter for CornStructor trees.

Updates
-------
v0.4.1
- Fix: add required GenBank annotations (molecule_type='DNA', topology='linear')
  to avoid "missing molecule_type in annotations" error from Biopython.
- Keep root fragment included as a Synthon feature when 'is_oligo' is False.

v0.4.0
- Do not skip the root fragment: always include it as a Synthon feature.

v0.3.1
- Export ONLY the nodes present in the tree (no extra inferred features).
- Map:  is_oligo == True  → Oligo_*  (lead/lag based on strand)
        is_oligo == False → Synthon_* (lead/lag based on strand)
- Write each node ONCE (lead for '+' strand, lag for '-' strand). No duplicates.
- Feature locations mirror GenBank style with complement(start..end) for '-' strand.
- Use 1-based closed intervals against the full reference/record sequence.
- Include /Sequence qualifier oriented to the feature strand (RC for lag).

Notes
-----
- 'start'/'end' are 0-based, end-exclusive, relative to the full record.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

from backend.app.core.models.fragment_node import FragmentNode


def _rc(s: str) -> str:
    """Reverse-complement (A/C/G/T only)."""
    return str(Seq(s).reverse_complement())


@dataclass
class _Counters:
    oligo: int = 0
    synthon: int = 0


def _iter_nodes_preorder(root: FragmentNode) -> Iterable[FragmentNode]:
    """Yield all nodes including the root (preorder)."""
    stack = [root]
    while stack:
        n = stack.pop()
        yield n
        if getattr(n, "children", None):
            for c in reversed(n.children):
                stack.append(c)


def _mk_feature_for_node(n: FragmentNode, counters: _Counters) -> SeqFeature | None:
    """
    Create a single GenBank SeqFeature for a node, choosing *lead* for '+' or
    *lag* for '-' strand, and mapping type by is_oligo.
    """
    start0 = int(getattr(n, "start", 0) or 0)
    end0 = int(getattr(n, "end", 0) or 0)
    if end0 <= start0:
        return None

    strand = -1 if str(getattr(n, "strand", "+")) == "-" else 1
    is_oligo = bool(getattr(n, "is_oligo", False))
    frag_id = str(getattr(n, "fragment_id", "")) or ""

    if is_oligo:
        counters.oligo += 1
        nth = counters.oligo
        base_key = "oligo"
    else:
        counters.synthon += 1
        nth = counters.synthon
        base_key = "synthon"

    role = "lead" if strand == 1 else "lag"
    feat_key = f"{base_key}_{role}"

    seq = str(getattr(n, "seq", "") or "")
    oriented_seq = seq if strand == 1 else _rc(seq)

    # Biopython uses 0-based, end-exclusive; strand stored on FeatureLocation.
    loc = FeatureLocation(start=start0, end=end0, strand=strand)

    std_name = f"{base_key.capitalize()}_{nth}_{role}"
    label = f"{base_key.capitalize()}_{nth}_{role}"

    qualifiers = {
        "Sequence": [oriented_seq],
        "label": [label],
        "note": [f"Geneious type: {base_key.capitalize()}_{role}"],
        "standard_name": [std_name],
    }
    if frag_id:
        qualifiers["fragment_id"] = [frag_id]

    return SeqFeature(location=loc, type=feat_key, qualifiers=qualifiers)


def _build_record(
    root: FragmentNode,
    record_id: str,
    definition: str | None = None,
) -> SeqRecord:
    """
    Build a single-record GenBank SeqRecord from the tree.

    The record sequence is taken from the *root* node's seq (the full design).
    """
    seq_text = str(getattr(root, "seq", "") or "")
    if not seq_text:
        seq_len = max(0, int(getattr(root, "end", 0) or 0) - int(getattr(root, "start", 0) or 0))
        seq_text = "N" * seq_len

    rec = SeqRecord(
        Seq(seq_text),
        id=record_id or "CornStructor",
        name=(record_id or "CornStructor")[:16],
        description=definition or ".",
    )

    # REQUIRED annotations for GenBank writer
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"] = "linear"
    # Optional cosmetics (safe defaults)
    rec.annotations.setdefault("data_file_division", "UNC")
    rec.annotations.setdefault("date", None)

    counters = _Counters()
    for node in _iter_nodes_preorder(root):
        feat = _mk_feature_for_node(node, counters)
        if feat is not None:
            rec.features.append(feat)

    return rec


def export_genbank_from_tree(
    root: FragmentNode,
    out_path,
    *,
    record_id: str | None = None,
    definition: str | None = None,
) -> None:
    """
    Write a GenBank file containing the full reference sequence (root.seq) and
    a FEATURES entry for each node (oligos and synthons), once per node, chosen
    as lead or lag according to its 'strand'.
    """
    rid = record_id or str(getattr(root, "fragment_id", "") or "CornStructor")
    rec = _build_record(root, rid, definition=definition)

    if hasattr(out_path, "write"):
        SeqIO.write(rec, out_path, "genbank")
    else:
        SeqIO.write(rec, str(out_path), "genbank")
