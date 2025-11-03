# File: backend/app/core/primer/designer.py
# Version: v1.4.0
"""
Primer search & pairing (honors JSON parameters; exact anchoring supported).

What this file does
-------------------
- Consumes `PrimerDesignParameters` (Pydantic) from `parameters.py` with *strict* new keys.
- If `templateExact` is True (recommended), designs primers **anchored** to the window:
  - Forward primer starts exactly at `start`
  - Reverse primer ends exactly at `end`
  - Product size = `end - start` exactly
- Enforces:
  - length in [primerLengthMin, primerLengthMax]
  - Tm in [primerTmMin, primerTmMax]   (Wallace rule by default)
  - GC% in [primerGCMin, primerGCMax]
  - max homopolymer run ≤ primerHomopolymerMax
  - |Tm_f - Tm_r| ≤ primerTmDifferenceMax
- Produces rich diagnostics on failure.

Non-exact mode
--------------
- If `templateExact` is False, this module currently raises NotImplementedError,
  keeping semantics explicit. (You can plug your exploration/relaxation search
  here later.)

Coordinates
-----------
- `start` and `end` are 0-based, with `end` exclusive.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Iterable, List, Tuple, Optional, Dict

from backend.app.core.primer.parameters import PrimerDesignParameters as ParamModel


# --- Utilities -----------------------------------------------------------------------------------

def rc(seq: str) -> str:
    """Reverse-complement (A<->T, C<->G)."""
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def gc_content(seq: str) -> float:
    s = seq.upper()
    if not s:
        return 0.0
    gc = s.count("G") + s.count("C")
    return 100.0 * gc / len(s)


def has_homopolymer(seq: str, max_run: int) -> bool:
    """Return True if any base repeats more than max_run times consecutively."""
    if not seq:
        return False
    run_char = seq[0]
    run_len = 1
    for ch in seq[1:]:
        if ch == run_char:
            run_len += 1
            if run_len > max_run:
                return True
        else:
            run_char = ch
            run_len = 1
    return False


def wallace_tm(seq: str) -> float:
    """Simple Tm proxy (°C): 2*(A+T) + 4*(G+C)."""
    s = seq.upper()
    a = s.count("A")
    t = s.count("T")
    g = s.count("G")
    c = s.count("C")
    return 2.0 * (a + t) + 4.0 * (g + c)


# --- Diagnostics / DTOs --------------------------------------------------------------------------

@dataclass
class CandidateRow:
    side: str                # 'F' or 'R'
    pos: int                 # +strand index (F: start; R: slice start on +strand)
    length: int
    seq: str
    tm: float
    gc: float
    rejected: bool
    reason: str


@dataclass
class DesignDiagnostics:
    forward_candidates: List[CandidateRow] = field(default_factory=list)
    reverse_candidates: List[CandidateRow] = field(default_factory=list)
    pair_scores_checked: int = 0
    message: str = ""


@dataclass
class Primer:
    seq: str
    tm: float
    gc: float
    pos: int
    length: int


@dataclass
class PrimerPair:
    forward: Primer
    reverse: Primer
    product_size: int
    score: float


# --- Validation helpers --------------------------------------------------------------------------

def _candidate_ok(
    seq: str,
    p: ParamModel,
    tm_func: Callable[[str], float],
) -> Tuple[bool, float, float, str]:
    """Validate single-primer constraints against ParamModel."""
    L = len(seq)
    if not (p.primerLengthMin <= L <= p.primerLengthMax):
        return False, 0.0, 0.0, f"length({L}) not in [{p.primerLengthMin},{p.primerLengthMax}]"

    tm = tm_func(seq)
    if not (p.primerTmMin <= tm <= p.primerTmMax):
        return False, tm, 0.0, f"Tm {tm:.1f} not in [{p.primerTmMin:.1f},{p.primerTmMax:.1f}]"

    gc = gc_content(seq)
    if not (p.primerGCMin <= gc <= p.primerGCMax):
        return False, tm, gc, f"GC {gc:.1f}% not in [{p.primerGCMin:.1f},{p.primerGCMax:.1f}]"

    if has_homopolymer(seq, p.primerHomopolymerMax):
        return False, tm, gc, f"homopolymer > {p.primerHomopolymerMax}"

    # NOTE: primerThreePrimeEndLength / match thresholds, secondary-structure ΔG,
    # and off-target metrics are reserved for a future pass.
    return True, tm, gc, ""


def _score_pair(f: Primer, r: Primer, p: ParamModel, product_size: int) -> float:
    """Lower is better; balance closeness and ΔTm (size is fixed in exact mode)."""
    d1 = abs((f.tm + r.tm) / 2.0 - (p.primerTmMin + p.primerTmMax) / 2.0)
    d2 = abs(f.tm - r.tm)
    # Small soft penalty for being near GC bounds
    g1 = min(abs(f.gc - p.primerGCMin), abs(p.primerGCMax - f.gc))
    g2 = min(abs(r.gc - p.primerGCMin), abs(p.primerGCMax - r.gc))
    return d1 + 2.0 * d2 + 0.02 * (g1 + g2)


# --- Public API ----------------------------------------------------------------------------------

def design(
    sequence: str,
    start: int,
    end: int,
    params: ParamModel,
    tm_func: Optional[Callable[[str], float]] = None,
) -> Tuple[PrimerPair, DesignDiagnostics]:
    """
    Entry point: honors `templateExact`.
    """
    if params.templateExact:
        return design_exact(sequence, start, end, params, tm_func=tm_func)
    raise NotImplementedError("Non-exact mode is disabled. Set templateExact=true in parameters JSON.")


def design_exact(
    sequence: str,
    start: int,
    end: int,
    params: ParamModel,
    tm_func: Optional[Callable[[str], float]] = None,
) -> Tuple[PrimerPair, DesignDiagnostics]:
    """
    Anchor primers exactly at window edges:
      - Forward primer: sequence[start : start + Lf]
      - Reverse primer: rc(sequence[end - Lr : end])
    """
    if start < 0 or end <= start or end > len(sequence):
        raise ValueError(f"Invalid coordinates: start={start}, end={end}, len={len(sequence)} (0-based, end-exclusive)")

    tmf = tm_func or wallace_tm
    diag = DesignDiagnostics()
    s = sequence.upper()

    # Bound lengths by available sequence
    min_Lf = params.primerLengthMin
    max_Lf = min(params.primerLengthMax, len(s) - start)
    min_Lr = params.primerLengthMin
    max_Lr = min(params.primerLengthMax, end - 0)
    if max_Lf < min_Lf or max_Lr < min_Lr:
        raise ValueError("Not enough sequence to place anchored primers with the requested length bounds.")

    forward_ok: List[Primer] = []
    reverse_ok: List[Primer] = []

    # Enumerate forward candidates at `start`
    for Lf in range(min_Lf, max_Lf + 1):
        fseq = s[start : start + Lf]
        ok, tm, gc, reason = _candidate_ok(fseq, params, tmf)
        diag.forward_candidates.append(CandidateRow("F", start, Lf, fseq, tm, gc, not ok, reason))
        if ok:
            forward_ok.append(Primer(fseq, tm, gc, start, Lf))

    # Enumerate reverse candidates that end at `end`
    for Lr in range(min_Lr, max_Lr + 1):
        if end - Lr < 0:
            continue
        r_raw = s[end - Lr : end]
        rseq = rc(r_raw)
        ok, tm, gc, reason = _candidate_ok(rseq, params, tmf)
        diag.reverse_candidates.append(CandidateRow("R", end - Lr, Lr, rseq, tm, gc, not ok, reason))
        if ok:
            reverse_ok.append(Primer(rseq, tm, gc, end - Lr, Lr))

    product_size = end - start
    best: Optional[PrimerPair] = None
    checked = 0
    for f in forward_ok:
        for r in reverse_ok:
            if abs(f.tm - r.tm) > params.primerTmDifferenceMax:
                continue
            score = _score_pair(f, r, params, product_size)
            checked += 1
            if (best is None) or (score < best.score):
                best = PrimerPair(f, r, product_size, score)
    diag.pair_scores_checked = checked

    if best:
        diag.message = "OK"
        return best, diag

    # Explain failure
    def top_reasons(rows: List[CandidateRow], side: str) -> List[Tuple[str, int]]:
        counts: Dict[str, int] = {}
        for r in rows:
            if r.side == side and r.rejected and r.reason:
                counts[r.reason] = counts.get(r.reason, 0) + 1
        return sorted(counts.items(), key=lambda x: x[1], reverse=True)[:5]

    f_reasons = top_reasons(diag.forward_candidates, "F")
    r_reasons = top_reasons(diag.reverse_candidates, "R")
    diag.message = "No valid anchored primer pair"

    hints = [
        "No valid anchored primer pair at exact boundaries.",
        f"- Forward@{start}: tested {len(diag.forward_candidates)}, ok={sum(1 for r in diag.forward_candidates if not r.rejected)}",
        f"- Reverse@end={end}: tested {len(diag.reverse_candidates)}, ok={sum(1 for r in diag.reverse_candidates if not r.rejected)}",
        "- Top forward rejection reasons: " + (", ".join(f"{k} x{v}" for k, v in f_reasons) or "n/a"),
        "- Top reverse rejection reasons: " + (", ".join(f"{k} x{v}" for k, v in r_reasons) or "n/a"),
        "Try widening primerTmMin/primerTmMax, primerGCMin/primerGCMax, ",
        "increasing primerLengthMin..primerLengthMax, or relaxing primerTmDifferenceMax.",
    ]
    raise ValueError("\n".join(hints))
