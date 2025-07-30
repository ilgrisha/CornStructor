# File: backend/app/core/visualization/oligo_html_visualizer.py
# Version: v0.6.5

"""
HTML visualization of full‐length oligo clusters.

For each cluster, this module renders:
  1. The exact original DNA fragment (5′→3′) colored by nucleotide.
  2. A list of oligo IDs, their strand orientation, and full oligo sequences.
  3. The “top” strand showing sense oligos aligned (5′→3′) with dashes for gaps.
  4. The “bottom” strand showing antisense oligos (reverse‐complemented) aligned (3′←5′) with gaps.
  5. A GA progress log indicating best fitness per generation.

All output uses a monospace font for alignment, and nucleotides are color‐coded.
"""

from pathlib import Path
from typing import List, Tuple, Dict


def reverse_complement(seq: str) -> str:
    """
    Compute the reverse complement of a DNA sequence.

    Args:
        seq: Input DNA sequence (A/C/G/T).

    Returns:
        The reverse complement (T/G/C/A, reversed).
    """
    # Translate bases then reverse the string
    complement = str.maketrans("ACGT", "TGCA")
    return seq.translate(complement)[::-1]


def color_nucleotide(nt: str) -> str:
    """
    Wrap a single nucleotide or dash in a colored <span>.

    Args:
        nt: Single character ('A','C','G','T','-').

    Returns:
        HTML <span> string with inline color styling.
    """
    # Define a distinct color for each base, gray for gaps
    cmap = {
        "A": "#d73027",   # red
        "C": "#4575b4",   # blue
        "G": "#1a9850",   # green
        "T": "#984ea3",   # purple
        "-": "#888888"    # grey for gaps
    }
    color = cmap.get(nt, "black")
    return f"<span style='color:{color};'>{nt}</span>"


def color_sequence(seq: str) -> str:
    """
    Apply color_nucleotide to each character in a sequence.

    Args:
        seq: String of nucleotides/dashes.

    Returns:
        Concatenated HTML string of colored spans.
    """
    return "".join(color_nucleotide(n) for n in seq)


def render_oligo_list(
    sid: str,
    ci: int,
    cluster: List[Tuple[str, str, int, int, str, str]]
) -> str:
    """
    Render a bullet‐list of all oligos in the cluster.

    Each cluster entry is a 6‐tuple:
      (full_seq, strand, start, end, overlap_prev, overlap_next)

    We only display:
      - Oligo ID (sid_c{ci}_o{oi})
      - Strand ("sense" or "antisense")
      - Full‐length oligo sequence colored and monospaced.

    Args:
        sid: Sequence identifier.
        ci:  Cluster index.
        cluster: List of tuples for each oligo.

    Returns:
        HTML <ul>...</ul> block.
    """
    items = []
    for oi, (seq, strand, start, end, *_overlaps) in enumerate(cluster, start=1):
        oid = f"{sid}_c{ci}_o{oi}"
        lbl = "sense" if strand == "+" else "antisense"
        # Wrap the colored sequence in a monospace span
        items.append(
            "<li>"
            f"<strong>{oid}</strong> ({lbl}): "
            f"<span class='mono'>{color_sequence(seq)}</span>"
            "</li>"
        )
    return "<ul class='oligo-list'>" + "".join(items) + "</ul>"


def visualize_cluster_html(
    sid: str,
    ci: int,
    cluster: List[Tuple[str, str, int, int, str, str]],
    full_seq: str
) -> str:
    """
    Visualize a single cluster as aligned top/bottom strands.

    - Extracts the exact genomic fragment spanning from the first oligo's start
      to the last oligo's end.
    - Initializes top/bottom arrays of '-' (gaps).
    - Fills in sense oligo sequences on the top, and reverse‐complemented
      antisense oligos on the bottom.
    - Colors and labels both strands with direction arrows.

    Args:
        sid: Sequence identifier (unused here, for context).
        ci:  Cluster index (unused here, for context).
        cluster: List of 6‐tuples per oligo.
        full_seq: The full original sequence for this sid.

    Returns:
        HTML block with three <div class='mono'> lines:
          1) original fragment
          2) top strand
          3) bottom strand
    """
    # Determine the fragment span in the original sequence
    fragment_start = cluster[0][2]
    fragment_end   = cluster[-1][3]
    fragment_seq   = full_seq[fragment_start:fragment_end]
    L = len(fragment_seq)

    # Original fragment line (5′→3′)
    orig_line = f"5′→ {color_sequence(fragment_seq)} 3′"

    # Prepare gap‐filled arrays
    top = ["-"] * L
    bot = ["-"] * L

    # Place each oligo into the appropriate strand
    for seq, strand, start, end, *_ in cluster:
        rel = start - fragment_start
        if strand == "+":
            # Fill top with sense sequence
            for i, nt in enumerate(seq):
                top[rel + i] = nt
        else:
            # Fill bottom with reverse‐complement of antisense
            rc = reverse_complement(seq)
            for i, nt in enumerate(rc):
                bot[rel + i] = nt

    # Color and join
    top_line = "5′→ " + color_sequence("".join(top)) + " 3′"
    bot_line = "3′← " + color_sequence("".join(bot)) + "  5′"

    return (
        f"<div class='mono'>{orig_line}</div>"
        f"<div class='mono'>{top_line}</div>"
        f"<div class='mono'>{bot_line}</div>"
    )


def format_progress_log(log: List[float]) -> str:
    """
    Render the GA fitness progression as an HTML list.

    Args:
        log: List of best-fitness values per generation.

    Returns:
        HTML <h4> + <ul>…</ul> block, or empty string if no log.
    """
    if not log:
        return ""
    items = "".join(f"<li>Gen {i+1}: {f:.1f}</li>" for i, f in enumerate(log))
    return "<h4>GA Progress</h4><ul>" + items + "</ul>"


def export_html_report(
    output_path: Path,
    clusters_by_sequence: Dict[str, List[List[Tuple[str, str, int, int, str, str]]]],
    ga_logs: Dict[str, List[float]] = None,
    full_sequences: Dict[str, str] = None
) -> None:
    """
    Generate a complete HTML report for all sequences and clusters.

    - Prints each sequence ID as a <h2>.
    - For each cluster: renders the oligo list and the aligned strands.
    - At the end of each sequence, appends the GA progress log if provided.

    Args:
        output_path: File path to write the HTML report.
        clusters_by_sequence: Mapping sid → list of clusters → list of oligo tuples.
        ga_logs: Optional mapping sid → fitness log list.
        full_sequences: Optional mapping sid → full DNA string.

    Side‐effects:
        Writes the assembled HTML to `output_path`.
    """
    # Begin HTML document and inject basic CSS for monospace and list styling
    html = [
        "<html><head><meta charset='UTF-8'><title>Oligo Design Visualization</title>",
        "<style>",
        ".mono{font-family:monospace;white-space:pre-wrap;}",
        ".oligo-list{margin:0 0 12px 16px;font-family:monospace;}",
        ".oligo-list li{margin-bottom:8px;}",
        "</style></head><body>",
        "<h1>Oligo Design Visualization</h1>"
    ]

    for sid, clusters in clusters_by_sequence.items():
        html.append(f"<h2>Sequence: {sid}</h2>")
        full_seq = full_sequences.get(sid, "")
        for ci, cluster in enumerate(clusters, start=1):
            html.append(f"<h3>Cluster {ci}</h3>")
            # Oligo summary list
            html.append(render_oligo_list(sid, ci, cluster))
            # Aligned-strands visualization
            html.append("<div style='background:#f8f8f8;padding:8px;margin-bottom:16px;border:1px solid #ccc;'>")
            html.append(visualize_cluster_html(sid, ci, cluster, full_seq))
            html.append("</div>")
        # Optionally include GA fitness progression
        if ga_logs and sid in ga_logs:
            html.append(format_progress_log(ga_logs[sid]))

    # Close document and write file
    html.append("</body></html>")
    output_path.write_text("\n".join(html), encoding="utf-8")
