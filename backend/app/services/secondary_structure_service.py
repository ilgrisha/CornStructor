# File: backend/app/services/secondary_structure_service.py
# Version: v0.1.0
"""
Secondary structure analysis service using ViennaRNA.

This module exposes a single function `analyze_stems` that:
1) Computes the minimum free energy (MFE) structure via ViennaRNA for a DNA sequence.
2) Marks nucleotides that participate in stems (paired in dot-bracket notation).
3) Extracts linear regions of consecutive paired nucleotides.
4) Merges adjacent stem regions if the unpaired gap between them is <= `merge_max_gap`.
5) Filters regions to keep only those whose final merged length >= `min_stem_len`.

All indices are 0-based, [start, end) half-open intervals to match the rest of CornStructor.
"""

from typing import List, Tuple

try:
    # ViennaRNA Python bindings
    import RNA  # type: ignore
except Exception as exc:  # pragma: no cover - import guard
    RNA = None  # Fallback handled in analyze_stems
    _IMPORT_ERROR = exc
else:
    _IMPORT_ERROR = None


def _rna_mfe_dotbracket(seq: str) -> str:
    """
    Compute MFE secondary structure (dot-bracket) using ViennaRNA for a DNA sequence.

    ViennaRNA is RNA-centric; to approximate DNA folding we apply a DNA-specific model
    by setting parameters on the fold compound:
      - md.temperature can be left default or set via env if needed
      - md.dangles, md.noLP, etc., use defaults suitable for MFE prediction

    NOTE: ViennaRNA supports soft settings for DNA via `md.betaScale`/parameter files.
    For simplicity and portability here, we use standard settings and treat the result
    as a heuristic stem detector for DNA.

    Args:
        seq: Uppercase DNA string (A/C/G/T). Other characters are passed-through to ViennaRNA.

    Returns:
        Dot-bracket structure string of same length as `seq`.
    """
    md = RNA.md()  # model details
    fc = RNA.fold_compound(seq, md)
    structure, mfe = fc.mfe()
    return structure


def _paired_flags_from_dotbracket(db: str) -> List[bool]:
    """
    Convert dot-bracket notation into per-nucleotide 'paired' flags.
    '(' and ')' are considered paired; '.' is unpaired.

    Args:
        db: dot-bracket structure string

    Returns:
        List[bool] of length len(db) where True denotes paired.
    """
    return [c in ("(", ")") for c in db]


def _runs_from_flags(flags: List[bool]) -> List[Tuple[int, int]]:
    """
    Extract consecutive True runs from a boolean list as [start, end) intervals.

    Args:
        flags: list of booleans

    Returns:
        List of (start, end) tuples (0-based, end-exclusive).
    """
    runs: List[Tuple[int, int]] = []
    n = len(flags)
    i = 0
    while i < n:
        if flags[i]:
            j = i + 1
            while j < n and flags[j]:
                j += 1
            runs.append((i, j))
            i = j
        else:
            i += 1
    return runs


def _merge_with_gap(runs: List[Tuple[int, int]], max_gap: int) -> List[Tuple[int, int]]:
    """
    Merge adjacent intervals if the unpaired gap between them is <= max_gap.

    Example:
        runs = [(10, 16), (18, 25)], max_gap=2  -> merge -> [(10, 25)]
        gap = next.start - prev.end

    Args:
        runs: sorted list of (start, end) intervals
        max_gap: maximum allowed unpaired length to merge across

    Returns:
        Merged list of (start, end).
    """
    if not runs:
        return []

    merged: List[Tuple[int, int]] = [runs[0]]
    for start, end in runs[1:]:
        prev_start, prev_end = merged[-1]
        gap = start - prev_end
        if gap <= max_gap:
            # merge
            merged[-1] = (prev_start, max(prev_end, end))
        else:
            merged.append((start, end))
    return merged


def analyze_stems(
    sequence: str,
    min_stem_len: int,
    merge_max_gap: int,
) -> List[Tuple[int, int]]:
    """
    Analyze DNA secondary structure stems and return linear intervals to visualize.

    Args:
        sequence: DNA sequence (uppercase 'A','C','G','T' preferred; others tolerated).
        min_stem_len: minimal number of *sequential nucleotides* (consecutive paired nts)
                      in a region to keep after merging.
        merge_max_gap: merge neighboring stem regions when the unpaired gap between them
                       is <= this value.

    Returns:
        List of (start, end) intervals (0-based, end-exclusive), sorted and non-overlapping.
    """
    if RNA is None:  # pragma: no cover - import guard path
        raise RuntimeError(
            "ViennaRNA (Python bindings) is not installed. Install the 'ViennaRNA' "
            "package and ensure system libraries are available."
        ) from _IMPORT_ERROR

    if not sequence:
        return []

    # Compute MFE structure
    db = _rna_mfe_dotbracket(sequence)

    # Convert to paired flags and extract runs
    flags = _paired_flags_from_dotbracket(db)
    runs = _runs_from_flags(flags)

    # Merge across short unpaired gaps
    merged = _merge_with_gap(runs, merge_max_gap)

    # Filter by minimal stem length
    filtered = [(s, e) for (s, e) in merged if (e - s) >= min_stem_len]

    return filtered
