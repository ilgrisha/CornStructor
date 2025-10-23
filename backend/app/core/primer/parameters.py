# File: backend/app/core/primer/parameters.py
# Version: v0.1.0
"""
Pydantic parameter models for primer design and target analysis.

These models are the single source of truth for default values and validation.
"""

from __future__ import annotations

from pydantic import BaseModel, Field, conint, confloat
from .constants import (
    DEFAULT_MIN_LEN,
    DEFAULT_MAX_LEN,
    DEFAULT_TM_MIN,
    DEFAULT_TM_MAX,
    DEFAULT_GC_MIN,
    DEFAULT_GC_MAX,
    DEFAULT_HOMOPOLYMER_MAX,
    DEFAULT_3P_END_LEN,
    DEFAULT_3P_PRIMER_MATCH_MAX,
    DEFAULT_TARGET_MATCH_MAX,
    DEFAULT_TARGET_MATCH_NUM_MAX,
    DEFAULT_SECSTRUCT_DG_MIN,
    DEFAULT_WEIGHTS,
)


class PrimerDesignWeights(BaseModel):
    """Weights used in composite scoring."""
    wTm: float = Field(DEFAULT_WEIGHTS["wTm"], ge=0)
    wGC: float = Field(DEFAULT_WEIGHTS["wGC"], ge=0)
    wDimer: float = Field(DEFAULT_WEIGHTS["wDimer"], ge=0)
    w3p: float = Field(DEFAULT_WEIGHTS["w3p"], ge=0)
    wOff: float = Field(DEFAULT_WEIGHTS["wOff"], ge=0)


class PrimerDesignParameters(BaseModel):
    """Tunable knobs for primer generation and filtering."""
    primerLengthMin: conint(ge=10, le=60) = DEFAULT_MIN_LEN
    primerLengthMax: conint(ge=10, le=60) = DEFAULT_MAX_LEN

    primerTmMin: confloat(ge=40.0, le=100.0) = DEFAULT_TM_MIN
    primerTmMax: confloat(ge=40.0, le=100.0) = DEFAULT_TM_MAX

    primerGCMin: confloat(ge=0.0, le=100.0) = DEFAULT_GC_MIN
    primerGCMax: confloat(ge=0.0, le=100.0) = DEFAULT_GC_MAX

    primerHomopolymerMax: conint(ge=2, le=12) = DEFAULT_HOMOPOLYMER_MAX

    primerThreePrimeEndLength: conint(ge=3, le=12) = DEFAULT_3P_END_LEN
    primerThreePrimePrimerMatchMax: conint(ge=0, le=12) = DEFAULT_3P_PRIMER_MATCH_MAX

    primerTargetMatchMax: conint(ge=0, le=40) = DEFAULT_TARGET_MATCH_MAX
    primerTargetMatchNumberMax: conint(ge=0, le=40) = DEFAULT_TARGET_MATCH_NUM_MAX

    primerSecondaryStructureDeltaGMin: float = DEFAULT_SECSTRUCT_DG_MIN

    # Pair constraints
    primerTmDifferenceMax: float = 3.0

    # Scoring weights
    weights: PrimerDesignWeights = Field(default_factory=PrimerDesignWeights)

    # Optional seed for deterministic tie-breaking
    randomSeed: int | None = None


class TargetAnalysisParameters(BaseModel):
    """Parameters for target region analysis and visualization."""
    targetGCMin: confloat(ge=0.0, le=100.0) = 35.0
    targetGCMax: confloat(ge=0.0, le=100.0) = 65.0
    targetSlidingWindowSize: conint(ge=5, le=200) = 20
    targetHomopolymerMax: conint(ge=2, le=12) = 4
    targetRepeatMax: conint(ge=2, le=20) = 5
