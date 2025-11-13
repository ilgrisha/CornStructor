# File: backend/app/cli/assemble_cli.py
# Version: v0.2.4

"""
Command-line interface for hierarchical assembly.

v0.2.4:
- Export leaf sequences to CSV (<root>_fragments.csv) and include a download link
  on the index page.

v0.2.3:
- Backward compatible logging flag: accept both --log-level and legacy --log.
- Prints a small startup banner with resolved args so it's obvious the run began.
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
from backend.app.core.export.csv_exporter import export_fragments_to_csv
from backend.app.core.export.json_exporter import export_tree_to_json
from backend.app.core.visualization.tree_html_exporter import export_tree_to_html
from backend.app.core.visualization.cluster_html_report import export_all_levels
from backend.app.core.export.analysis_exporter import analyze_tree_to_json
from backend.app.core.visualization.analysis_html_report import export_analysis_html
from backend.app.core.export.ga_progress_exporter import export_ga_progress_json
from backend.app.core.visualization.ga_progress_html_report import export_ga_progress_html
from backend.app.core.visualization.report_index_html import export_run_index_html


def main() -> None:
    p = argparse.ArgumentParser(description="Hierarchical assembly CLI")
    p.add_argument("--fasta", required=True, type=Path, help="Input FASTA with one or more records")
    p.add_argument("--levels", required=True, type=Path, help="levels.json")
    p.add_argument("--globals", required=False, type=Path, help="globals.json")
    p.add_argument("--outdir", required=True, type=Path, help="Output directory")
    # Accept BOTH --log-level and legacy --log (they map to the same dest)
    p.add_argument("--log-level", dest="log_level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                   help="Logging level (default: INFO)")
    p.add_argument("--log", dest="log_level", choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                   help=argparse.SUPPRESS)
    args = p.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level))
    log = logging.getLogger("assemble_cli")

    # Startup banner
    log.info("=== const_cli assemble ===")
    log.info("FASTA=%s | LEVELS=%s | GLOBALS=%s | OUTDIR=%s | LOG=%s",
             str(args.fasta), str(args.levels), str(args.globals) if args.globals else "None",
             str(args.outdir), args.log_level)

    fasta_path = Path(args.fasta)
    levels_path = Path(args.levels)
    globals_path = Path(args.globals) if args.globals else None
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    levels = load_levels_config(levels_path)
    gconf = load_global_config(globals_path)

    n_ok = 0
    n_fail = 0
    index_paths = []

    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        root_id = rec.id
        full_seq = str(rec.seq).upper()
        log.info("Processing %s (%d bp)", root_id, len(full_seq))

        try:
            asm = HierarchicalAssembler(levels, gconf, logger=log)
            root = asm.build(full_seq, root_id)

            # Core exports (FASTA + CSV)
            export_fragments_to_fasta(root, outdir / f"{root_id}_fragments.fasta", root_id)
            export_fragments_to_csv(root,  outdir / f"{root_id}_fragments.csv",   root_id)

            # Tree & analysis
            export_tree_to_json(root, outdir / f"{root_id}_tree.json", root_id)
            export_tree_to_html(root, outdir / f"{root_id}_tree.html", root_id)

            analysis_json_path = outdir / f"{root_id}_analysis.json"
            analysis_html_path = outdir / f"{root_id}_analysis.html"
            analysis_payload = analyze_tree_to_json(root, analysis_json_path, root_id)
            export_analysis_html(analysis_payload, analysis_html_path)

            # GA progress: JSON + HTML
            ga_json_path = outdir / f"{root_id}_ga_progress.json"
            ga_html_path = outdir / f"{root_id}_ga_progress.html"
            ga_payload = export_ga_progress_json(root, ga_json_path, root_id)
            export_ga_progress_html(ga_payload, ga_html_path)

            # Per-level cluster reports
            report_dir = outdir / f"{root_id}_cluster_reports"
            export_all_levels(root, report_dir, levels_cfg=levels, global_cfg=gconf)

            # Index page with links to everything + FASTA/CSV download
            index_path = export_run_index_html(root_id, outdir)
            log.info("Index written: %s", index_path)
            index_paths.append(index_path)

            log.info("✓ Completed %s", root_id)
            n_ok += 1
        except AssemblerError as e:
            log.error("✗ Failed %s: %s", root_id, e)
            n_fail += 1

    # Final summary / links
    if index_paths:
        print("")  # spacer
        print("Report index:")
        for pth in index_paths:
            abs_path = Path(pth).resolve()
            file_url = f"file://{abs_path.as_posix()}"
            print(f"  • {abs_path}")
            print(f"    {file_url}")
        print("")

    if n_fail > 0:
        raise SystemExit(1)
    raise SystemExit(0)


if __name__ == "__main__":
    main()
