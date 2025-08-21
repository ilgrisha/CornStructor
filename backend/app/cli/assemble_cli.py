# File: backend/app/cli/assemble_cli.py
# Version: v0.1.2

"""
Command-line interface for hierarchical assembly.

v0.1.2:
- Updated imports: configuration loaders moved to backend.app.config.*
- No behavior change.

v0.1.1:
- Switch to export_all_levels (per-level report set) instead of the old single
  export_clusters_and_fragments_html function.
"""

import argparse
import logging
from pathlib import Path

from Bio import SeqIO

from backend.app.config.config_level import load_levels_config
from backend.app.config.config_global import load_global_config
from backend.app.core.assembly.hierarchical_assembler import (
    HierarchicalAssembler,
    AssemblerError,
)
from backend.app.core.export.fasta_exporter import export_fragments_to_fasta
from backend.app.core.export.json_exporter import export_tree_to_json
from backend.app.core.visualization.tree_html_exporter import export_tree_to_html
from backend.app.core.visualization.cluster_html_report import export_all_levels


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--fasta", required=True, type=Path)
    p.add_argument("--levels", required=True, type=Path)
    p.add_argument("--globals", required=False, type=Path, help="Path to globals.json")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    args = p.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level))
    log = logging.getLogger("assemble_cli")

    fasta_path = Path(args.fasta)
    levels_path = Path(args.levels)
    globals_path = Path(args.globals) if args.globals else None
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    levels = load_levels_config(levels_path)
    gconf = load_global_config(globals_path)

    # Process each FASTA record
    n_ok = 0
    n_fail = 0
    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        root_id = rec.id
        full_seq = str(rec.seq).upper()
        log.info("Processing %s (%d bp)", root_id, len(full_seq))

        try:
            asm = HierarchicalAssembler(levels, gconf, logger=log)
            root = asm.build(full_seq, root_id)

            export_fragments_to_fasta(root, outdir / f"{root_id}_fragments.fasta", root_id)
            export_tree_to_json(root, outdir / f"{root_id}_tree.json", root_id)
            export_tree_to_html(root, outdir / f"{root_id}_tree.html", root_id)

            # Per-level cluster reports in a directory
            report_dir = outdir / f"{root_id}_cluster_reports"
            export_all_levels(root, report_dir, levels_cfg=levels, global_cfg=gconf)

            log.info("✓ Completed %s", root_id)
            n_ok += 1
        except AssemblerError as e:
            log.error("✗ Failed %s: %s", root_id, e)
            n_fail += 1

    if n_fail > 0:
        raise SystemExit(1)
    raise SystemExit(0)


if __name__ == "__main__":
    main()
