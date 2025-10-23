# File: backend/app/core/primer/schemas.py
# Version: v0.2.0
"""
DTOs for requests and responses used by Primer endpoints and services.

Update v0.2.0:
- `PrimerDesignRequest.parameters` is now optional.
  If omitted, backend will use the currently stored parameters
  (from backend/app/config/primers_param.json, with fallback to defaults).
"""

from __future__ import annotations

from pydantic import BaseModel, Field, conint
from typing import List, Optional, Dict, Any
from .parameters import PrimerDesignParameters, TargetAnalysisParameters


class PrimerDesignRequest(BaseModel):
    """Request to design primers within a target window of the given sequence."""
    sequence: str = Field(..., description="Full template sequence (raw or FASTA content pre-parsed to raw).")
    seqTargetStartPosition: conint(ge=1) = Field(..., description="1-based inclusive.")
    seqTargetEndPosition: conint(ge=1) = Field(..., description="1-based inclusive.")
    # Now optional: if not provided, server uses stored parameters
    parameters: Optional[PrimerDesignParameters] = None


class PrimerSeqInfo(BaseModel):
    """Basic primer properties."""
    sequence: str
    tm: float
    gc: float


class PrimerDesignResponse(BaseModel):
    """Result of primer design run."""
    runId: str
    forwardPrimer: PrimerSeqInfo
    reversePrimer: PrimerSeqInfo
    pairScore: float
    warnings: List[str] = Field(default_factory=list)


class TargetAnalysisRequest(BaseModel):
    sequence: str
    parameters: TargetAnalysisParameters = Field(default_factory=TargetAnalysisParameters)


class TargetAnalysisResponse(BaseModel):
    gcDistribution: List[float]
    homopolymers: List[Dict[str, Any]]
    repeats: List[Dict[str, Any]]
    structureWarnings: List[str]


class AlignmentsRequest(BaseModel):
    forward: str
    reverse: str
    target: str


class AlignmentsResponse(BaseModel):
    forwardTarget: str
    reverseTarget: str
    selfDimer: str
    crossDimer: str


class PrimerRunRecord(BaseModel):
    """For GET /runs and GET /runs/{run_id}."""
    id: str
    createdAt: str
    targetStart: int
    targetEnd: int
    sequenceDigest: str
    status: str
    result: Optional[PrimerDesignResponse] = None
