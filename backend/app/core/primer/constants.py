# File: backend/app/core/primer/constants.py
# Version: v0.1.0
"""
Constants and defaults for the Primer subsystem.

These values mirror the PrimeLime defaults so that logic, workflow, and parameters
can be shared across applications. Adjust as necessary for CornStructor.
"""

from __future__ import annotations

DEFAULT_MIN_LEN = 18
DEFAULT_MAX_LEN = 25

DEFAULT_TM_MIN = 55.0
DEFAULT_TM_MAX = 62.0

DEFAULT_GC_MIN = 40.0
DEFAULT_GC_MAX = 65.0

DEFAULT_HOMOPOLYMER_MAX = 4

DEFAULT_3P_END_LEN = 5
DEFAULT_3P_PRIMER_MATCH_MAX = 2

DEFAULT_TARGET_MATCH_MAX = 5
DEFAULT_TARGET_MATCH_NUM_MAX = 2

# Placeholder for ΔG/secondary structure defaults; implemented in later step
DEFAULT_SECSTRUCT_DG_MIN = -9.0  # kcal/mol (example placeholder)

# Default scoring weights (composite score components)
DEFAULT_WEIGHTS = {
    "wTm": 1.0,
    "wGC": 0.5,
    "wDimer": 1.0,
    "w3p": 1.25,
    "wOff": 1.0,
}
