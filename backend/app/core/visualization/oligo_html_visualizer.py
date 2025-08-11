# File: backend/app/core/visualization/oligo_html_visualizer.py
# Version: v0.7.6

"""
HTML visualization of full‐length oligo clusters.

Features:
  - Horizontal scrolling of sequences and alignments.
  - Per-cluster fragment visualization:
      * absolute-position ruler (labels + ticks), 1-based for display
      * original sequence
      * top/bottom aligned strands
  - Colored nucleotides.
  - Oligo list with IDs, orientation, and ABSOLUTE coordinate badges.
  - GA progress log (optional).
  - Responsive layout with maximum visible nucleotide width before scrolling.

Notes:
  - Lines are stacked (not side-by-side) but scroll together horizontally.
  - Pass `abs_start` to `visualize_cluster_html` (absolute start of the parent fragment).
  - If `abs_start` is None in render_oligo_list, badges show relative (“+”) coordinates.
  - The [left – right] bracket is NOT rendered inside the alignment; callers should place it in the cluster title.
"""

from pathlib import Path
from typing import List, Tuple, Dict, Optional

MAX_VISIBLE_NT = 150  # Controls visible nucleotide width before scrolling (approximate)
TICK_STEP      = 10   # tick every 10 bp
LABEL_STEP     = 50   # numeric label every 50 bp


def reverse_complement(seq: str) -> str:
    complement = str.maketrans("ACGT", "TGCA")
    return seq.translate(complement)[::-1]


def color_nucleotide(nt: str) -> str:
    cmap = {
        "A": "#d73027",  # red
        "C": "#4575b4",  # blue
        "G": "#1a9850",  # green
        "T": "#984ea3",  # purple
        "-": "#888888"   # grey for gaps
    }
    return f"<span style='color:{cmap.get(nt, 'black')}'>{nt}</span>"


def color_sequence(seq: str) -> str:
    return "".join(color_nucleotide(n) for n in seq)


def _build_position_ruler(abs_left_display: int, length: int) -> Tuple[str, str]:
    """
    Build two monospace lines (1-based display):
      - labels line: absolute bp numbers every LABEL_STEP
      - ticks  line: '|' every TICK_STEP
      - Ticks every 10 bp (|) AND at position 1.
      - Labels at 1 and every 50 bp (50, 100, 150, ...).
      - Additionally, the *first visible* label is forced to the next multiple of 10
        (ceil to 10), even if it's not a 50-multiple (e.g., 121 → label at 130).

    Args:
        abs_left_display: 1-based coordinate of the first visible base
        length:           number of visible characters

    Returns:
        (label_line, tick_line) as plain text; caller wraps in <div class='mono'>.
    """
    ticks  = [" "] * length
    labels = [" "] * length

    # ticks: at 1 and every 10 aligned to absolute coordinates
    for abs_pos in range(abs_left_display, abs_left_display + length):
        offset = abs_pos - abs_left_display
        if abs_pos == 1 or abs_pos % TICK_STEP == 0:
            ticks[offset] = "|"

    # First visible label: ceil to next multiple of 10 (unless abs_left_display == 1)
    if abs_left_display == 1:
        first_round10 = 1
    else:
        rem = abs_left_display % 10
        first_round10 = abs_left_display if rem == 0 else abs_left_display + (10 - rem)

    # Labels at 1, every 50, and first_round10 (once)
    first_label_placed = False
    for abs_pos in range(abs_left_display, abs_left_display + length):
        offset = abs_pos - abs_left_display
        is_label = False
        if abs_pos == 1:
            is_label = True
        elif abs_pos % LABEL_STEP == 0:
            is_label = True
        elif not first_label_placed and abs_pos == first_round10:
            is_label = True
            first_label_placed = True

        if is_label:
            num_str = str(abs_pos)
            for j, ch in enumerate(num_str):
                if offset + j < length:
                    labels[offset + j] = ch

    label_line = "    " + "".join(labels)
    tick_line  = "    " + "".join(ticks)
    return label_line, tick_line


def render_oligo_list(
    sid: str,
    ci: int,
    cluster: List[Tuple[str, str, int, int, str, str]],
    abs_start: Optional[int] = None
) -> str:
    """
    Render a bullet list of all oligos with ABSOLUTE coordinate badges (1-based).

    Args:
        sid:       sequence/fragment identifier for display
        ci:        cluster index
        cluster:   list of (seq_full, strand, start_rel, end_rel, prev_ov, next_ov)
        abs_start: absolute start (0-based) of the PARENT fragment. If None, badges
                   will show relative positions as +start..+end.

    Note:
        For absolute badges we display inclusive ranges:
          start_abs_display = abs_start + start_rel + 1
          end_abs_display   = abs_start + end_rel
    """
    badge_style = (
        "display:inline-block;margin-left:8px;padding:1px 6px;"
        "background:#eef;border:1px solid #99c;border-radius:10px;"
        "font-size:11px;color:#245;"
    )

    items = []
    for oi, (seq, strand, start_rel, end_rel, *_ov) in enumerate(cluster, start=1):
        oid = f"{sid}_c{ci}_o{oi}"
        lbl = "  sense  " if strand == "+" else "antisense"

        if abs_start is not None:
            start_abs_display = abs_start + start_rel + 1  # 1-based inclusive start
            end_abs_display   = abs_start + end_rel        # 1-based inclusive end
            badge = f"<span style='{badge_style}'>{start_abs_display}–{end_abs_display}</span>"
        else:
            # fallback: relative (user-visible +offsets)
            badge = f"<span style='{badge_style}'>+{start_rel+1}–+{end_rel}</span>"

        items.append(
            "<li>"
            f"<strong>{oid}</strong> ({lbl}) "
            f"{badge} : "
            f"<span class='mono'>{color_sequence(seq)}</span>"
            "</li>"
        )

    return (
        "<div class='scroll-box'>"
        "<ul class='oligo-list mono'>"
        + "".join(items) +
        "</ul></div>"
    )


def visualize_cluster_html(
    sid: str,
    ci: int,
    cluster: List[Tuple[str, str, int, int, str, str]],
    full_seq: str,
    abs_start: Optional[int] = None
) -> str:
    """
    Render a cluster with an absolute-position ruler (1-based display).

    Args:
        sid:        identifier used for oligo IDs (display only)
        ci:         cluster index (display only)
        cluster:    (seq_full, strand, start_rel, end_rel, prev_ov, next_ov),
                    start/end are RELATIVE to the parent fragment start.
        full_seq:   the parent fragment sequence for this cluster node (string)
        abs_start:  absolute genomic start (0-based) of the parent fragment;
                    used to display 1-based absolute positions.

    Note:
        The [left – right] bracket is intentionally NOT included here; callers should
        place it in the cluster title next to the fragment/cluster heading.
    """
    # Visible subrange (relative to parent)
    fragment_start_rel = min(c[2] for c in cluster)
    fragment_end_rel   = max(c[3] for c in cluster)
    fragment_seq       = full_seq[fragment_start_rel:fragment_end_rel]
    L = len(fragment_seq)

    # Absolute display coordinates (1-based)
    abs0 = abs_start if abs_start is not None else 0
    abs_left0      = abs0 + fragment_start_rel   # 0-based
    abs_left_disp  = abs_left0 + 1               # 1-based

    # Position ruler (labels & ticks)
    label_line, tick_line = _build_position_ruler(abs_left_disp, L)

    # Original and aligned lines (monospace)
    orig_line = f"5′→ {color_sequence(fragment_seq)} →3′"

    top = ["-"] * L
    bot = ["-"] * L
    for seq, strand, start_rel, _end_rel, *_ in cluster:
        rel = start_rel - fragment_start_rel
        if strand == "+":
            for i, nt in enumerate(seq):
                if 0 <= rel + i < L:
                    top[rel + i] = nt
        else:
            rc = reverse_complement(seq)
            for i, nt in enumerate(rc):
                if 0 <= rel + i < L:
                    bot[rel + i] = nt

    top_line = "5′→ " + color_sequence("".join(top)) + " →3′"
    bot_line = "3′← " + color_sequence("".join(bot)) + " ←5′"

    # Lines stacked; scroll together. Extra vertical space between blocks.
    return (
        "<div class='scroll-box'>"
        f"<div class='mono'>{label_line}</div>"
        f"<div class='mono'>{tick_line}</div>"
        f"<div class='mono'>{orig_line}</div>"
        "<div style='height:6px'></div>"
        f"<div class='mono'>{top_line}</div>"
        f"<div class='mono'>{bot_line}</div>"
        "<div style='height:6px'></div>"
        "</div>"
    )


def format_progress_log(log: List[float]) -> str:
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
    Generic one-page report (used less now that you have per-level pages).
    Keeps existing behavior; no abs_start context here, so badges are relative.
    """
    html = [
        "<html><head><meta charset='UTF-8'><title>Oligo Design Visualization</title>",
        "<style>",
        ".mono{font-family:monospace;white-space:pre;display:block;}",
        ".oligo-list{margin:0 0 12px 16px;font-family:monospace;}",
        ".oligo-list li{margin-bottom:8px;}",
        f".scroll-box{{overflow-x:auto;white-space:nowrap;padding:8px;background:#f9f9f9;"
        f"border:1px solid #ccc;margin-bottom:16px;max-width:{MAX_VISIBLE_NT}ch;}}",
        "body{font-family:sans-serif;padding:20px;background:#fcfcfc;color:#333;}",
        "h1{font-size:28px;border-bottom:2px solid #ccc;padding-bottom:8px;}",
        "h2{margin-top:30px;color:#2a4d69;}",
        "h3{margin-top:20px;color:#3e5c76;}",
        "</style></head><body>",
        "<h1>🧬 Oligo Design Visualization</h1>"
    ]

    for sid, clusters in clusters_by_sequence.items():
        html.append(f"<h2>Sequence: {sid}</h2>")
        full_seq = full_sequences.get(sid, "")
        for ci, cluster in enumerate(clusters, start=1):
            # Compute bracket from relative visible region (1-based display)
            frag_start_rel = min(c[2] for c in cluster)
            frag_end_rel   = max(c[3] for c in cluster)
            L = frag_end_rel - frag_start_rel
            left_disp  = frag_start_rel + 1
            right_disp = left_disp + L - 1
            bracket = f"[{left_disp} – {right_disp}]"

            html.append(f"<h3>Cluster {ci} {bracket}</h3>")
            # Relative badges only (no abs_start context here)
            html.append(render_oligo_list(sid, ci, cluster, abs_start=None))
            html.append(visualize_cluster_html(sid, ci, cluster, full_seq, abs_start=None))
        if ga_logs and sid in ga_logs:
            html.append(format_progress_log(ga_logs[sid]))

    html.append("</body></html>")
    output_path.write_text("\n".join(html), encoding="utf-8")
