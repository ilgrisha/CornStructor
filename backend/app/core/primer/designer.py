# File: backend/app/core/primer/designer.py
# Version: v0.2.0
"""
PrimerDesigner orchestrates candidate generation, filtering, off-target checks,
and composite scoring to select the best primer pair.

This version implements the full pipeline using local modules:
- generator: create forward/reverse candidates in the target window
- thermodynamics: Tm (NN) and GC%
- constraints: homopolymer, self-/cross-dimer, 3' complementarity
- offtarget: ungapped approximate local attachments
- scoring: composite score selection
"""

from __future__ import annotations

import hashlib
import math
import random
from dataclasses import dataclass
from typing import List, Tuple

from .parameters import PrimerDesignParameters
from .schemas import PrimerDesignResponse, PrimerSeqInfo
from .generator import generate_forward_candidates, generate_reverse_candidates
from .thermodynamics import tm_nearest_neighbor, gc_percent, revcomp, hairpin_dg_proxy
from .constraints import (
    longest_homopolymer,
    check_dimer_risk,
    three_prime_match_count,
)
from .offtarget import count_offtargets
from .scoring import ScoreContext, score_pair


@dataclass
class PrimerDesignResult:
    """Internal result used before persistence/DTO conversion."""
    forward_seq: str
    forward_tm: float
    forward_gc: float
    reverse_seq: str
    reverse_tm: float
    reverse_gc: float
    pair_score: float
    warnings: List[str]


class PrimerDesigner:
    """Main orchestration class for primer design."""
    def __init__(self) -> None:
        pass

    def design(
        self,
        sequence: str,
        start: int,
        end: int,
        params: PrimerDesignParameters,
    ) -> PrimerDesignResult:
        """
        Execute full design pipeline:
          1) Candidate generation (forward/reverse)
          2) Single-primer filters: length, Tm band, GC band, homopolymer, ΔG proxy
          3) Pairwise filters: ΔTm, self-/cross-dimer, 3' complementarity
          4) Off-target counting (full + 3' window)
          5) Composite scoring and selection
        """
        target = sequence[start - 1 : end].upper()
        if params.randomSeed is not None:
            random.seed(params.randomSeed)

        # 1) Generate candidates
        f_candidates = generate_forward_candidates(target, params.primerLengthMin, params.primerLengthMax)
        r_candidates = generate_reverse_candidates(target, params.primerLengthMin, params.primerLengthMax)

        def _single_ok(seq: str) -> Tuple[bool, List[str], float, float]:
            warns: List[str] = []
            tm = tm_nearest_neighbor(seq)
            gc = gc_percent(seq)
            if not (params.primerTmMin <= tm <= params.primerTmMax):
                return False, warns, tm, gc
            if not (params.primerGCMin <= gc <= params.primerGCMax):
                return False, warns, tm, gc
            if longest_homopolymer(seq) > params.primerHomopolymerMax:
                warns.append("Homopolymer exceeds limit.")
                return False, warns, tm, gc
            dg = hairpin_dg_proxy(seq)
            if dg < params.primerSecondaryStructureDeltaGMin:
                warns.append("Secondary structure ΔG proxy too strong.")
                return False, warns, tm, gc
            return True, warns, tm, gc

        # Filter singles
        f_pool: List[Tuple[str, float, float]] = []
        r_pool: List[Tuple[str, float, float]] = []
        for s in f_candidates:
            ok, _, tm, gc = _single_ok(s)
            if ok:
                f_pool.append((s, tm, gc))
        for s in r_candidates:
            ok, _, tm, gc = _single_ok(s)
            if ok:
                r_pool.append((s, tm, gc))

        if not f_pool or not r_pool:
            # Relaxed fallback: take extremes to return something usable
            # (still deterministic via sorting)
            f_pool = f_pool or [(f_candidates[0], tm_nearest_neighbor(f_candidates[0]), gc_percent(f_candidates[0]))]
            r_pool = r_pool or [(r_candidates[-1], tm_nearest_neighbor(r_candidates[-1]), gc_percent(r_candidates[-1]))]

        # 2) Pairwise evaluation
        best: Tuple[float, Tuple[str, float, float, str, float, float, List[str]]] | None = None

        ctx = ScoreContext(
            tm_min=params.primerTmMin,
            tm_max=params.primerTmMax,
            gc_min=params.primerGCMin,
            gc_max=params.primerGCMax,
            wTm=params.weights.wTm,
            wGC=params.weights.wGC,
            wDimer=params.weights.wDimer,
            w3p=params.weights.w3p,
            wOff=params.weights.wOff,
            pair_dt_max=params.primerTmDifferenceMax,
        )

        for f_seq, f_tm, f_gc in f_pool:
            # Precompute some self metrics
            f_self_consec, f_self_total = check_dimer_risk(f_seq, f_seq)
            for r_seq, r_tm, r_gc in r_pool:
                # Pair constraint: ΔTm
                if abs(f_tm - r_tm) > params.primerTmDifferenceMax + 5.0:
                    # small slack to avoid over-pruning at this stage; final scoring will penalize
                    continue

                # Cross-dimer metrics
                r_rc = revcomp(r_seq)
                c_consec, c_total = check_dimer_risk(f_seq, r_seq)

                # 3' end complementarity checks (self & cross)
                f_3p_self = three_prime_match_count(f_seq, f_seq, params.primerThreePrimeEndLength)
                f_3p_cross = three_prime_match_count(f_seq, r_rc, params.primerThreePrimeEndLength)

                # Off-target counts within the target region (both strands)
                off_full, off_3p = count_offtargets(
                    f_seq, target, params.primerTargetMatchMax, params.primerThreePrimeEndLength, params.primerThreePrimeTargetMatchMax
                    if hasattr(params, "primerThreePrimeTargetMatchMax") else params.primerThreePrimePrimerMatchMax
                )

                score = score_pair(
                    ctx,
                    f_tm, f_gc,
                    r_tm, r_gc,
                    f_self_consec, f_self_total,
                    c_consec, c_total,
                    f_3p_self, f_3p_cross,
                    off_full, off_3p,
                )

                pack = (f_seq, f_tm, f_gc, r_seq, r_tm, r_gc, [])  # warnings filled below as needed
                if best is None or score < best[0]:
                    best = (score, pack)

        assert best is not None
        score, (f_seq, f_tm, f_gc, r_seq, r_tm, r_gc, warns) = best

        # Collect warnings vs. params bands
        if f_gc < params.primerGCMin or f_gc > params.primerGCMax:
            warns.append("Forward primer GC out of bounds.")
        if r_gc < params.primerGCMin or r_gc > params.primerGCMax:
            warns.append("Reverse primer GC out of bounds.")
        if abs(f_tm - r_tm) > params.primerTmDifferenceMax:
            warns.append("ΔTm between primers exceeds recommended maximum.")

        return PrimerDesignResult(
            forward_seq=f_seq,
            forward_tm=f_tm,
            forward_gc=f_gc,
            reverse_seq=r_seq,
            reverse_tm=r_tm,
            reverse_gc=r_gc,
            pair_score=score,
            warnings=warns,
        )

    @staticmethod
    def digest_sequence(seq: str) -> str:
        return hashlib.sha256(seq.encode("utf-8")).hexdigest()
