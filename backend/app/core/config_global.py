# File: backend/app/core/config_global.py
# Version: v0.1.1

"""
Global configuration for genetic algorithm and assembly behavior.

This JSON is passed alongside the per-level config. Example:

{
  "ga": {
    "population_size": 240,
    "num_generations": 8,
    "mutation_rate_initial": 0.35,
    "crossover_rate": 0.6,
    "elitism_count": 2,
    "random_injection_rate": 0.05,
    "tournament_size": 3,
    "enable_hill_climb": true,
    "hill_climb_period": 5,
    "hill_climb_junc": 2,
    "hill_climb_cands": 10,
    "allowed_motif_scale": 1.0
  },
  "cpu_workers_fraction": 0.9,
  "enforce_span_across_junction": true,
  "max_total_adjustment_rounds": 2,
  "division_adjust_attempts_per_junction": 4,
  "division_adjust_step_bp": 5,
  "tm_method": "NN",
  "tm": {
    "Na": 0.05,
    "K": 0.0,
    "Mg": 0.0015,
    "dNTPs": 0.0,
    "dnac1": 2.5e-7,
    "dnac2": 2.5e-7,
    "saltcorr": 7
  }
}
"""

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Dict

from backend.app.core.assembly.ga_overlap_selector import GAParameters


@dataclass
class GlobalConfig:
    """
    Container for global (cross-level) knobs used by the assembler and GA.

    Notes:
      - `ga` uses default_factory to avoid mutable default issues with dataclasses.
      - Tm parameters mirror Biopython's MeltingTemp NN parameters.
    """
    # GA knobs
    ga: GAParameters = field(default_factory=GAParameters)

    # CPU usage fraction (0 < f ≤ 1)
    cpu_workers_fraction: float = 0.9

    # Require overlaps to span the junction within ±max_gap_allowed
    enforce_span_across_junction: bool = True

    # Division-point adjustment policy
    max_total_adjustment_rounds: int = 2
    division_adjust_attempts_per_junction: int = 4
    division_adjust_step_bp: int = 5

    # Biopython Tm method & parameters
    tm_method: str = "NN"  # or "Wallace"
    Na: float = 50e-3
    K: float = 0.0
    Mg: float = 0.0
    dNTPs: float = 0.0
    dnac1: float = 250e-9
    dnac2: float = 250e-9
    saltcorr: int = 7

    def tm_params_dict(self) -> Dict[str, float]:
        """Return a dict of Tm parameters for Biopython MeltingTemp calls."""
        return {
            "Na": self.Na, "K": self.K, "Mg": self.Mg, "dNTPs": self.dNTPs,
            "dnac1": self.dnac1, "dnac2": self.dnac2, "saltcorr": self.saltcorr
        }

    def ga_params(self) -> GAParameters:
        """Return the GAParameters object."""
        return self.ga


def load_global_config(path: Optional[Path]) -> GlobalConfig:
    """
    Load a GlobalConfig from JSON. Missing fields fall back to defaults.
    """
    if path is None:
        return GlobalConfig()
    data = json.loads(path.read_text(encoding="utf-8"))

    # GA params
    ga_data = data.get("ga", {})
    ga = GAParameters(
        population_size=ga_data.get("population_size", 200),
        num_generations=ga_data.get("num_generations", 3),
        mutation_rate_initial=ga_data.get("mutation_rate_initial", 0.35),
        crossover_rate=ga_data.get("crossover_rate", 0.6),
        elitism_count=ga_data.get("elitism_count", 1),
        random_injection_rate=ga_data.get("random_injection_rate", 0.05),
        tournament_size=ga_data.get("tournament_size", 3),
        enable_hill_climb=ga_data.get("enable_hill_climb", True),
        hill_climb_period=ga_data.get("hill_climb_period", 5),
        hill_climb_junc=ga_data.get("hill_climb_junc", 2),
        hill_climb_cands=ga_data.get("hill_climb_cands", 10),
        allowed_motif_scale=ga_data.get("allowed_motif_scale", 1.0),
    )

    gc = GlobalConfig(
        ga=ga,
        cpu_workers_fraction=float(data.get("cpu_workers_fraction", 0.9)),
        enforce_span_across_junction=bool(data.get("enforce_span_across_junction", True)),
        max_total_adjustment_rounds=int(data.get("max_total_adjustment_rounds", 2)),
        division_adjust_attempts_per_junction=int(data.get("division_adjust_attempts_per_junction", 4)),
        division_adjust_step_bp=int(data.get("division_adjust_step_bp", 5)),
        tm_method=str(data.get("tm_method", "NN")),
        Na=float(data.get("tm", {}).get("Na", 50e-3)),
        K=float(data.get("tm", {}).get("K", 0.0)),
        Mg=float(data.get("tm", {}).get("Mg", 0.0)),
        dNTPs=float(data.get("tm", {}).get("dNTPs", 0.0)),
        dnac1=float(data.get("tm", {}).get("dnac1", 250e-9)),
        dnac2=float(data.get("tm", {}).get("dnac2", 250e-9)),
        saltcorr=int(data.get("tm", {}).get("saltcorr", 7)),
    )
    return gc
