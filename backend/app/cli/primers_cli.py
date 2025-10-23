# File: backend/app/cli/primers_cli.py
# Version: v0.2.0
"""
CLI entrypoint for primer design.

Usage:
    python -m app.cli.primers_cli \
        --fasta input.fasta \
        --start 150 \
        --end 420 \
        --outdir /path/to/out \
        [--params-json /path/to/params.json]

Behavior:
- If --params-json is omitted, loads current parameters from backend/app/config/primers_param.json
  (with fallback to primers_param_default.json).
- Writes primers.fasta to the output directory with Forward and Reverse sequences.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional
from backend.app.core.primer.designer import PrimerDesigner
from backend.app.core.primer.parameters import PrimerDesignParameters
from backend.app.config.config_primers import load_current_params

def read_fasta(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    lines = [ln.strip() for ln in text.splitlines() if ln and not ln.startswith(">")]
    return "".join(lines)

def write_primers_fasta(path: Path, fwd: str, rev: str) -> None:
    path.write_text(f">Forward_Primer\n{fwd}\n>Reverse_Primer\n{rev}\n", encoding="utf-8")

def load_params_from_file(path: Optional[Path]) -> PrimerDesignParameters:
    if path and path.exists():
        import json
        return PrimerDesignParameters.model_validate(json.loads(path.read_text(encoding="utf-8")))
    # default to stored parameters (or defaults)
    return load_current_params()

def main():
    p = argparse.ArgumentParser(description="PrimeLime-like primer design CLI for CornStructor.")
    p.add_argument("--fasta", required=True, help="Input FASTA file with template sequence.")
    p.add_argument("--start", type=int, required=True, help="1-based inclusive start of target window.")
    p.add_argument("--end", type=int, required=True, help="1-based inclusive end of target window.")
    p.add_argument("--outdir", required=True, help="Output directory.")
    p.add_argument("--params-json", required=False, help="Optional JSON file with PrimerDesignParameters.")
    args = p.parse_args()

    fasta_path = Path(args.fasta)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    seq = read_fasta(fasta_path)
    params = load_params_from_file(Path(args.params_json) if args.params_json else None)

    designer = PrimerDesigner()
    result = designer.design(seq, args.start, args.end, params)

    out_fa = outdir / "primers.fasta"
    write_primers_fasta(out_fa, result.forward_seq, result.reverse_seq)
    print(f"Wrote {out_fa}")

if __name__ == "__main__":
    main()
