# File: backend/app/cli/genbank_from_tree.py
# Version: v0.2.1
"""
CLI: Build a GenBank file from a CornStructor tree.json + a reference FASTA.

Usage
-----
python -m backend.app.cli.genbank_from_tree \
  --fasta /path/to/reference.fasta \
  --tree  /path/to/tree.json \
  --outdir ./out

Notes
-----
- Supports tree.json shapes:
    A) { "root_id": "...", "tree": { ...node... } }   (current)
    B) { "root_id": "...", "root": { ...node... } }   (older)
    C) { ...node... }                                 (bare node)
- Emits a single-record GenBank with FEATURES for all nodes (root included).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

from Bio import SeqIO

from backend.app.core.models.fragment_node import FragmentNode
from backend.app.core.export.genbank_exporter import export_genbank_from_tree


def _coerce_strand(s: Optional[str]) -> str:
    return "-" if (str(s or "+").strip() == "-") else "+"


def _node_from_json(obj: Dict[str, Any], full_seq: str) -> FragmentNode:
    """
    Build a FragmentNode from a JSON object. Missing 'seq' is derived
    from reference using [start:end] (0-based, end-exclusive).
    """
    start = int(obj.get("start", 0) or 0)
    end = int(obj.get("end", 0) or 0)
    seq = obj.get("seq")
    if not seq:
        s = max(0, min(start, len(full_seq)))
        e = max(s, min(end, len(full_seq)))
        seq = full_seq[s:e]

    n = FragmentNode(
        fragment_id=str(obj.get("fragment_id", "")),
        level=int(obj.get("level", 0) or 0),
        start=start,
        end=end,
        seq=str(seq or ""),
        strand=_coerce_strand(obj.get("strand")),
        overlap_prev=str(obj.get("overlap_prev", "") or ""),
        overlap_next=str(obj.get("overlap_next", "") or ""),
        is_oligo=bool(obj.get("is_oligo", False)),
    )

    children = obj.get("children") or []
    n.children = [_node_from_json(ch, full_seq) for ch in children] if isinstance(children, list) else []

    for k in ("ga_log", "ga_detail", "ga_cluster_id"):
        if k in obj:
            setattr(n, k, obj[k])

    return n


def _load_tree_json(tree_path: Path, full_seq: str) -> tuple[FragmentNode, str]:
    """
    Load CornStructor's tree.json supporting multiple shapes:

      { "root_id": "...", "tree": { ... } }
      { "root_id": "...", "root": { ... } }
      { ...node... }  # bare node fallback

    Returns (root_node, root_id).
    """
    data = json.loads(tree_path.read_text(encoding="utf-8"))

    if isinstance(data, dict) and "tree" in data:
        root_obj = data["tree"]
        root_id = str(data.get("root_id") or root_obj.get("fragment_id") or tree_path.stem)
    elif isinstance(data, dict) and "root" in data:
        root_obj = data["root"]
        root_id = str(data.get("root_id") or root_obj.get("fragment_id") or tree_path.stem)
    else:
        root_obj = data if isinstance(data, dict) else {}
        root_id = str(root_obj.get("fragment_id") or tree_path.stem)

    root_node = _node_from_json(root_obj, full_seq)
    if not getattr(root_node, "fragment_id", None):
        setattr(root_node, "fragment_id", root_id)

    return root_node, root_id


def _read_reference_fasta(fasta_path: Path) -> str:
    """Read the first record sequence from a FASTA file."""
    with fasta_path.open("r", encoding="utf-8") as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            return str(rec.seq)
    raise ValueError(f"No sequences found in FASTA: {fasta_path}")


def _default_output_name(root_id: str) -> str:
    base = root_id or "CornStructor"
    base = "".join(ch if (ch.isalnum() or ch in "-_") else "_" for ch in base)
    return f"{base}.gb"


def parse_args(argv: List[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="genbank-from-tree",
        description="Generate a GenBank file from CornStructor tree.json and a reference FASTA.",
    )
    p.add_argument("--fasta", required=True, type=Path, help="Path to reference FASTA")
    p.add_argument("--tree", required=True, type=Path, help="Path to CornStructor tree.json")
    p.add_argument("--outdir", required=True, type=Path, help="Directory to write the GenBank into")
    p.add_argument("--record-id", default=None, help="Override GenBank record id/name")
    p.add_argument("--definition", default=None, help="GenBank DEFINITION line text")
    return p.parse_args(argv)


def main(argv: List[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])

    if not args.fasta.exists():
        print(f"[ERROR] FASTA not found: {args.fasta}", file=sys.stderr)
        return 2
    if not args.tree.exists():
        print(f"[ERROR] tree.json not found: {args.tree}", file=sys.stderr)
        return 2
    args.outdir.mkdir(parents=True, exist_ok=True)

    try:
        full_seq = _read_reference_fasta(args.fasta)
    except Exception as e:
        print(f"[ERROR] Failed reading FASTA: {e}", file=sys.stderr)
        return 3

    try:
        root, parsed_root_id = _load_tree_json(args.tree, full_seq)
    except Exception as e:
        print(f"[ERROR] Failed parsing tree.json: {e}", file=sys.stderr)
        return 4

    root_id = args.record_id or parsed_root_id or getattr(root, "fragment_id", None) or "CornStructor"
    out_name = _default_output_name(root_id)
    out_path = args.outdir / out_name

    try:
        export_genbank_from_tree(
            root,
            out_path,
            record_id=str(root_id),
            definition=args.definition or "CornStructor assembly.",
        )
    except Exception as e:
        print(f"[ERROR] GenBank export failed: {e}", file=sys.stderr)
        return 5

    print(f"[OK] Wrote GenBank → {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
