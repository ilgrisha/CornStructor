# File: backend/app/api/v1/primers/router.py
# Version: v0.2.0
"""
Primer endpoints:
- POST /design
- POST /target-analysis
- POST /alignments
- GET /runs
- GET /runs/{run_id}
- GET /parameters         ← returns current primer design parameters
- PUT /parameters         ← validates & persists new parameters
"""

from __future__ import annotations

from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session
from typing import List

from backend.app.core.primer.designer import PrimerDesigner
from backend.app.core.primer.analysis import TargetAnalyzer
from backend.app.core.primer.alignments import AlignmentRenderer
from backend.app.core.primer.schemas import (
    PrimerDesignRequest,
    PrimerDesignResponse,
    TargetAnalysisRequest,
    TargetAnalysisResponse,
    AlignmentsRequest,
    AlignmentsResponse,
    PrimerRunRecord,
    PrimerSeqInfo,
)
from backend.app.core.primer.parameters import PrimerDesignParameters
from backend.app.config.config_primers import load_current_params, save_current_params, ensure_current_exists

from .models import PrimerRun, PrimerResult
from .deps import db_session

router = APIRouter(prefix="/api/v1/primers", tags=["primers"])


@router.get("/parameters", response_model=PrimerDesignParameters)
def get_parameters():
    """
    Return the current editable primer design parameters.
    If not initialized, create primers_param.json from defaults and return it.
    """
    _, params = ensure_current_exists()
    return params


@router.put("/parameters", response_model=PrimerDesignParameters)
def update_parameters(payload: PrimerDesignParameters):
    """
    Validate and persist new primer design parameters into primers_param.json.
    """
    save_current_params(payload)
    return payload


@router.post("/design", response_model=PrimerDesignResponse)
def design_primers(
    payload: PrimerDesignRequest,
    db: Session = Depends(db_session),
):
    """
    Design primers on the selected target window (Start/End 1-based inclusive).
    If `parameters` is omitted, the server uses the stored parameters
    from backend/app/config/primers_param.json (with default fallback).
    """
    # Extract
    seq = payload.sequence.replace("\n", "").replace("\r", "")
    start = payload.seqTargetStartPosition
    end = payload.seqTargetEndPosition
    params = payload.parameters or load_current_params()

    if start < 1 or end < start or end > len(seq):
        raise HTTPException(status_code=400, detail="Invalid target coordinates.")

    # Design
    designer = PrimerDesigner()
    result = designer.design(seq, start, end, params)

    # Persist
    digest = designer.digest_sequence(seq)
    run = PrimerRun(
        sequence_digest=digest,
        target_start=start,
        target_end=end,
        parameters_json=params.model_dump(),
        status="completed",
    )
    db.add(run)
    db.flush()

    res = PrimerResult(
        run_id=run.id,
        forward_seq=result.forward_seq,
        forward_tm=str(result.forward_tm),
        forward_gc=str(result.forward_gc),
        reverse_seq=result.reverse_seq,
        reverse_tm=str(result.reverse_tm),
        reverse_gc=str(result.reverse_gc),
        pair_score=str(result.pair_score),
        warnings_json=result.warnings,
        alignment_blocks_json=None,
    )
    db.add(res)
    db.commit()
    db.refresh(run)

    return PrimerDesignResponse(
        runId=run.id,
        forwardPrimer=PrimerSeqInfo(sequence=result.forward_seq, tm=result.forward_tm, gc=result.forward_gc),
        reversePrimer=PrimerSeqInfo(sequence=result.reverse_seq, tm=result.reverse_tm, gc=result.reverse_gc),
        pairScore=result.pair_score,
        warnings=result.warnings,
    )


@router.post("/target-analysis", response_model=TargetAnalysisResponse)
def target_analysis(
    payload: TargetAnalysisRequest,
):
    """Analyze GC windows, repeats, homopolymers, structure flags."""
    analyzer = TargetAnalyzer()
    data = analyzer.analyze(payload.sequence, payload.parameters)
    return TargetAnalysisResponse(**data)


@router.post("/alignments", response_model=AlignmentsResponse)
def alignments(payload: AlignmentsRequest):
    """Return ASCII blocks for forward/target, reverse/targetRC, self- and cross-dimer."""
    r = AlignmentRenderer()
    return AlignmentsResponse(
        forwardTarget=r.forward_vs_target(payload.forward, payload.target),
        reverseTarget=r.reverse_vs_target_rc(payload.reverse, payload.target),
        selfDimer=r.self_dimer(payload.forward),
        crossDimer=r.cross_dimer(payload.forward, payload.reverse),
    )


@router.get("/runs", response_model=List[PrimerRunRecord])
def list_runs(db: Session = Depends(db_session)):
    """List recent primer design runs (lightweight view)."""
    runs = db.query(PrimerRun).order_by(PrimerRun.created_at.desc()).limit(100).all()
    results = []
    for r in runs:
        results.append(
            PrimerRunRecord(
                id=r.id,
                createdAt=r.created_at.isoformat() if r.created_at else "",
                targetStart=r.target_start,
                targetEnd=r.target_end,
                sequenceDigest=r.sequence_digest,
                status=r.status,
                result=None,  # lightweight
            )
        )
    return results


@router.get("/runs/{run_id}", response_model=PrimerRunRecord)
def get_run(run_id: str, db: Session = Depends(db_session)):
    """Return a full run, including result if present."""
    run = db.get(PrimerRun, run_id)
    if not run:
        raise HTTPException(status_code=404, detail="Run not found.")
    out = PrimerRunRecord(
        id=run.id,
        createdAt=run.created_at.isoformat() if run.created_at else "",
        targetStart=run.target_start,
        targetEnd=run.target_end,
        sequenceDigest=run.sequence_digest,
        status=run.status,
        result=None,
    )
    if run.result:
        res = run.result
        out.result = PrimerDesignResponse(
            runId=run.id,
            forwardPrimer=PrimerSeqInfo(sequence=res.forward_seq, tm=float(res.forward_tm), gc=float(res.forward_gc)),
            reversePrimer=PrimerSeqInfo(sequence=res.reverse_seq, tm=float(res.reverse_tm), gc=float(res.reverse_gc)),
            pairScore=float(res.pair_score),
            warnings=list(res.warnings_json or []),
        )
    return out
