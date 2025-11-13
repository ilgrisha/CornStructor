# File: backend/app/api/v1/secondary_structure.py
# Version: v0.1.0
"""
API router for DNA secondary structure analysis using ViennaRNA.

POST /api/v1/analysis/stems
  - Body: SecondaryStructureRequest
  - Returns: SecondaryStructureResponse with regions of stems as [start, end) intervals

Security: No authentication is enforced here by default. If your project uses
dependency-based auth (e.g., JWT), wrap the endpoint with Depends as appropriate.
"""

from fastapi import APIRouter, HTTPException
from ...schemas.secondary_structure import (
    SecondaryStructureRequest,
    SecondaryStructureResponse,
    FeatureRegion,
)
from ...services.secondary_structure_service import analyze_stems

router = APIRouter(prefix="/analysis", tags=["analysis"])


@router.post("/stems", response_model=SecondaryStructureResponse)
async def analyze_stems_endpoint(payload: SecondaryStructureRequest) -> SecondaryStructureResponse:
    """
    Analyze stems via ViennaRNA and return merged regions for visualization.

    This endpoint follows CornStructor's convention of [start, end) 0-based indices and
    the "FeatureRegion" shape used elsewhere in the frontend overlays.
    """
    seq = payload.sequence.upper()
    try:
        intervals = analyze_stems(
            sequence=seq,
            min_stem_len=payload.min_stem_len,
            merge_max_gap=payload.merge_max_gap,
        )
    except RuntimeError as exc:
        # Typically raised if ViennaRNA is not installed or misconfigured.
        raise HTTPException(status_code=500, detail=str(exc)) from exc

    regions = [FeatureRegion(start=s, end=e) for (s, e) in intervals]
    return SecondaryStructureResponse(length=len(seq), regions=regions)
