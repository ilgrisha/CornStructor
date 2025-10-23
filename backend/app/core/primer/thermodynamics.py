# File: backend/app/core/primer/thermodynamics.py
# Version: v0.1.0
"""
Thermodynamics utilities for primer properties.

Implements:
- Nearest-neighbor Tm (BioPython)
- GC percentage
- Reverse complement
- Basic ΔG proxy for secondary structure screening (lightweight heuristic)

Notes:
- For production-grade ΔG hairpin/dimer energies, integrate a full model
  (e.g., ViennaRNA or detailed NN parameter set). Here we implement a
  pragmatic heuristic that's fast and effective as a first pass.
"""

from __future__ import annotations

from typing import Tuple

try:
    # BioPython >=1.78 MeltingTemp
    from Bio.SeqUtils import MeltingTemp as mt
except Exception:  # pragma: no cover - optional import guard
    mt = None  # type: ignore


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(table)[::-1]


def gc_percent(seq: str) -> float:
    if not seq:
        return 0.0
    s = seq.upper()
    gc = sum(1 for c in s if c in ("G", "C"))
    return 100.0 * gc / len(s)


def tm_nearest_neighbor(seq: str, Na: float = 50.0, divalent: float = 0.0, dntp: float = 0.0, dna_conc: float = 250.0) -> float:
    """
    Melting temperature using NN method (°C).

    Args:
        seq: primer sequence (A/C/G/T)
        Na: monovalent salt concentration (mM)
        divalent: divalent salt (mM)
        dntp: dNTP concentration (mM)
        dna_conc: primer strand concentration (nM) - used by mt.Tm_NN as annealing conc

    Returns:
        Tm in °C
    """
    if not seq:
        return 0.0
    if mt is None:  # fallback: Wallace rule
        return 2.0 * sum(1 for c in seq if c.upper() in ("A", "T")) + 4.0 * sum(1 for c in seq if c.upper() in ("G", "C"))
    # Convert to molar units expected by BioPython:
    # mt.Tm_NN assumes concentrations in mol/L; inputs here are in mM/nM.
    _Na = Na / 1000.0
    _Mg = divalent / 1000.0
    _dntp = dntp / 1000.0
    _dna = dna_conc * 1e-9
    # SantaLucia parameters
    return float(mt.Tm_NN(seq, Na=_Na, Mg=_Mg, dNTPs=_dntp, dnac1=_dna, dnac2=_dna))


def hairpin_dg_proxy(seq: str, window: int = 4) -> float:
    """
    Very lightweight ΔG proxy for hairpins/self-structure.

    Heuristic:
    - Scan small windows at 3' and internal regions for complementary runs.
    - Penalize longer consecutive complementary matches more.
    - Return a negative value in kcal/mol scale (approx), where more negative == more stable (worse).

    This is NOT a physical ΔG; it's a fast screen to exclude obvious strong structures.

    Returns:
        Proxy ΔG (negative is stronger hairpin).
    """
    if not seq or len(seq) < window + 3:
        return 0.0
    s = seq.upper()
    rc = revcomp(s)
    # Simple scan: best local consecutive match when aligning s against rc with small shifts
    best = 0
    for shift in range(2, min(20, len(s))):  # a few offsets to simulate loop sizes
        consec = 0
        local_best = 0
        for i in range(min(len(s), len(rc) - shift)):
            if s[i] == rc[i + shift]:
                consec += 1
                local_best = max(local_best, consec)
            else:
                consec = 0
        best = max(best, local_best)
    # Map best consecutive comp length to an approximate ΔG proxy
    # Roughly: -0.8 kcal/mol per consecutive base pair (toy mapping)
    return -0.8 * best
