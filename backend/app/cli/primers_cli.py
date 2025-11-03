# File: backend/app/cli/primers_cli.py
# Version: v1.4.0
"""
CLI for CornStructor primer design (strict parameters, exact anchoring supported).

- Parameters JSON must use the new schema (camelCase), including "templateExact".
- If --exact is given, it forces anchored design regardless of JSON.
- Writes primers.fasta and primers.json; with --debug also writes candidates.csv + README.txt.

Usage:
    python -m backend.app.cli.primers_cli \
        --fasta backend/data/input/region.fasta \
        --start 150 --end 420 \
        --outdir backend/data/out/primers \
        --params-json backend/app/config/primers_param.json \
        [--one-based] [--exact] [--debug]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Tuple

# ✅ FIXED: import the parameters model from parameters.py (not designer.py)
from backend.app.core.primer.parameters import PrimerDesignParameters
from backend.app.core.primer.designer import (
    design,
    design_exact,
    DesignDiagnostics,
)

# ---------- IO helpers ----------

def read_single_fasta(path: Path) -> Tuple[str, str]:
    """Return (name, sequence). Enforces exactly one FASTA record."""
    name = None
    seq_lines = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    raise ValueError(f"Multiple FASTA records found in {path}. Provide a single-sequence FASTA.")
                name = line[1:].strip() or "sequence"
            else:
                seq_lines.append(line)
    if name is None:
        raise ValueError(f"No FASTA header found in {path}.")
    return name, "".join(seq_lines).replace(" ", "").replace("\t", "").upper()


def write_fasta(out: Path, forward_seq: str, reverse_seq: str, tm_f: float, tm_r: float, gc_f: float, gc_r: float) -> None:
    out.write_text(
        (
            f">Forward_Primer tm={tm_f:.1f} gc={gc_f:.1f} len={len(forward_seq)}\n"
            f"{forward_seq}\n"
            f">Reverse_Primer tm={tm_r:.1f} gc={gc_r:.1f} len={len(reverse_seq)}\n"
            f"{reverse_seq}\n"
        ),
        encoding="utf-8",
    )


def dump_debug_csv(outdir: Path, diag: DesignDiagnostics) -> None:
    import csv
    cand_file = outdir / "candidates.csv"
    with cand_file.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["side", "pos", "length", "tm", "gc", "rejected", "reason", "seq"])
        for row in diag.forward_candidates + diag.reverse_candidates:
            w.writerow([row.side, row.pos, row.length, f"{row.tm:.2f}", f"{row.gc:.2f}", row.rejected, row.reason, row.seq])
    (outdir / "README.txt").write_text(
        "Diagnostics files:\n"
        "- candidates.csv: per-candidate filters and reasons\n"
        f"- pairs_scored: {diag.pair_scores_checked}\n"
        f"- message: {diag.message}\n",
        encoding="utf-8",
    )


# ---------- Main ----------

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--fasta", required=True, type=Path)
    p.add_argument("--start", required=True, type=int, help="Target region start (0-based, inclusive) unless --one-based")
    p.add_argument("--end", required=True, type=int, help="Target region end (exclusive) unless --one-based")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--params-json", required=True, type=Path, help="Path to JSON with PrimerDesignParameters (strict, camelCase)")
    p.add_argument("--one-based", action="store_true", help="Treat start/end as 1-based inclusive coordinates")
    p.add_argument("--debug", action="store_true", help="Write diagnostics CSV files")
    p.add_argument("--exact", action="store_true", help="Force anchored primers at start/end (overrides JSON)")
    args = p.parse_args()

    try:
        # FASTA
        name, seq = read_single_fasta(args.fasta)

        # Coordinates
        start = args.start
        end = args.end
        if args.one_based:
            if end < start:
                raise ValueError("End must be >= start for 1-based inclusive coordinates.")
            start = start - 1  # 1-based inclusive -> 0-based inclusive
            # end remains as-is and becomes 0-based exclusive

        if not (0 <= start < end <= len(seq)):
            raise ValueError(f"Coordinates out of bounds for {name}: start={start}, end={end}, len={len(seq)} (0-based)")

        # Parameters (strict, no legacy)
        params_data = json.loads(args.params_json.read_text(encoding="utf-8"))
        params = PrimerDesignParameters.model_validate(params_data)

        args.outdir.mkdir(parents=True, exist_ok=True)

        # Respect --exact override, otherwise follow JSON's templateExact
        force_exact = args.exact
        if force_exact and not params.templateExact:
            # Create a copy with templateExact True without mutating the source dict
            params = params.model_copy(update={"templateExact": True})

        if params.templateExact:
            pair, diag = design_exact(seq, start, end, params)
        else:
            # Will raise NotImplementedError by design to keep semantics explicit
            pair, diag = design(seq, start, end, params)

        # Outputs
        out_fa = args.outdir / "primers.fasta"
        write_fasta(
            out_fa,
            pair.forward.seq,
            pair.reverse.seq,
            pair.forward.tm,
            pair.reverse.tm,
            pair.forward.gc,
            pair.reverse.gc,
        )

        meta = {
            "sequence_name": name,
            "coords": {"start": start, "end": end, "length": len(seq)},
            "product_size": pair.product_size,
            "forward": {"pos": pair.forward.pos, "len": pair.forward.length, "tm": pair.forward.tm, "gc": pair.forward.gc},
            "reverse": {"pos": pair.reverse.pos, "len": pair.reverse.length, "tm": pair.reverse.tm, "gc": pair.reverse.gc},
            "score": pair.score,
            "params_source": str(args.params_json),
            "anchored": bool(params.templateExact),
        }
        (args.outdir / "primers.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

        if args.debug:
            dump_debug_csv(args.outdir, diag)

        mode = "anchored" if params.templateExact else "exploratory"
        print(f"[OK] Wrote {out_fa} (product {pair.product_size} bp, {mode}). Params: {args.params_json}")

    except Exception as ex:
        print(f"[ERROR] {ex}", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
