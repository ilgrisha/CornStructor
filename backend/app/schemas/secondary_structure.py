# File: backend/app/schemas/secondary_structure.py
# Version: v0.1.0
"""
Pydantic schemas for secondary structure (stems) analysis.
"""

from typing import List, Literal
from pydantic import BaseModel, Field, constr


class SecondaryStructureRequest(BaseModel):
    """Request payload for secondary structure analysis."""
    sequence: constr(strip_whitespace=True, min_length=1) = Field(
        ...,
        description="DNA sequence to analyze (A/C/G/T preferred; case-insensitive).",
        examples=["ACGTTGCA..."],
    )
    min_stem_len: int = Field(
        4,
        ge=1,
        description="Minimal length (bp) for a stem region (after merging) to keep.",
    )
    merge_max_gap: int = Field(
        2,
        ge=0,
        description="Merge adjacent stem regions if the unpaired gap between them is <= this value.",
    )


class FeatureRegion(BaseModel):
    """Region to visualize on the sequence viewer."""
    kind: Literal["stems"] = "stems"
    start: int = Field(..., ge=0, description="0-based inclusive index.")
    end: int = Field(..., ge=0, description="0-based exclusive index; must be >= start.")


class SecondaryStructureResponse(BaseModel):
    """Response containing merged stem regions suitable for visualization."""
    length: int = Field(..., ge=0, description="Length of the input sequence.")
    regions: List[FeatureRegion] = Field(default_factory=list)
