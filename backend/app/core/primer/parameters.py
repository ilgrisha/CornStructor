# File: backend/app/core/primer/parameters.py
# Version: v1.0.0
"""
Pydantic model for primer design parameters (strict, no legacy keys).

This model mirrors the JSON schema you provided, including `templateExact`.
Some fields (e.g., 3' end match, secondary structure ΔG, off-target thresholds)
are not currently enforced by `designer.py` but are retained for future use.

Usage:
    from backend.app.core.primer.parameters import PrimerDesignParameters
"""

from __future__ import annotations

from typing import Optional
from pydantic import BaseModel, Field, conint, confloat


class Weights(BaseModel):
    """Scoring weights (reserved for future multi-factor scoring)."""
    wTm: float = 1.0
    wGC: float = 0.5
    wDimer: float = 1.0
    w3p: float = 1.25
    wOff: float = 1.0


class PrimerDesignParameters(BaseModel):
    # Lengths
    primerLengthMin: conint(ge=6) = Field(18, description="Minimum primer length")
    primerLengthMax: conint(gt=6) = Field(28, description="Maximum primer length")

    # Temperatures (range)
    primerTmMin: confloat(ge=0) = Field(55.0, description="Minimum acceptable primer Tm (°C)")
    primerTmMax: confloat(ge=0) = Field(62.0, description="Maximum acceptable primer Tm (°C)")

    # GC content (%)
    primerGCMin: confloat(ge=0, le=100) = Field(40.0, description="Minimum GC percentage")
    primerGCMax: confloat(ge=0, le=100) = Field(65.0, description="Maximum GC percentage")

    # Structure/sequence constraints
    primerHomopolymerMax: conint(ge=1) = Field(4, description="Max run of identical bases")
    primerThreePrimeEndLength: conint(ge=0) = Field(5, description="3' segment length for match checks")
    primerThreePrimePrimerMatchMax: conint(ge=0) = Field(2, description="Max matches allowed in 3' segment")
    primerTargetMatchMax: conint(ge=0) = Field(5, description="Max matches allowed to target outside site")
    primerTargetMatchNumberMax: conint(ge=0) = Field(2, description="Max number of such matches")
    primerSecondaryStructureDeltaGMin: float = Field(-9.0, description="Min ΔG for hairpins/dimers (not yet enforced)")

    # Pairing constraint
    primerTmDifferenceMax: confloat(ge=0) = Field(3.0, description="Max |Tm_f - Tm_r| (°C)")

    # Scoring weights (reserved)
    weights: Weights = Field(default_factory=Weights)

    # Reproducibility
    randomSeed: Optional[int] = Field(None, description="Optional RNG seed")

    # NEW: exact anchoring
    templateExact: bool = Field(True, description="If true, forward starts at 'start' and reverse ends at 'end'")

    # (Optional) bounds sanity
    def model_post_init(self, __context) -> None:  # pydantic v2 hook
        if self.primerLengthMax < self.primerLengthMin:
            raise ValueError("primerLengthMax must be >= primerLengthMin")
        if self.primerTmMax < self.primerTmMin:
            raise ValueError("primerTmMax must be >= primerTmMin")
        if self.primerGCMax < self.primerGCMin:
            raise ValueError("primerGCMax must be >= primerGCMin")
