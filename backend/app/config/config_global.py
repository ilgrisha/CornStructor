# File: backend/app/config/config_global.py
# Version: v0.2.3
"""
Global configuration loader.

v0.2.3
- Change `ga_params()` to return an attribute-accessible view object (GAParamsView)
  instead of a dict, so callers like GAOverlapSelector can use dot-access
  (e.g., `.mutation_rate_initial`). Keeps `ga_params_dict` for dict-based callers.

v0.2.2
- Added `ga_params()` and `ga_params_dict` alias.

v0.2.1
- Added `tm_params_dict()` and `tm_params` property.

v0.2.0
- Moved from backend.app.core.config_global -> backend.app.config.config_global.
- **Removed** `enforce_span_across_junction` (it was a no-op in GA v1.2.0).
- Kept GA knobs and Tm chemistry parameters. Added optional Tris buffer.

This module reads `globals.json` and returns a strongly-typed GlobalConfig.
If `path` is None, we default to the `globals.json` file sitting next to this
module (backend/app/config/globals.json).
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Dict, Any


# ----------------------------- Data classes ----------------------------------


@dataclass(frozen=True)
class GAConfig:
    """Genetic Algorithm configuration knobs."""
    population_size: int
    num_generations: int
    mutation_rate_initial: float
    crossover_rate: float
    elitism_count: int
    random_injection_rate: float
    tournament_size: int = 3
    enable_hill_climb: bool = True
    hill_climb_period: int = 5
    hill_climb_junc: int = 2
    hill_climb_cands: int = 10
    allowed_motif_scale: float = 1.0


@dataclass(frozen=True)
class GAParamsView:
    """
    Read-only, attribute-accessible view of GA parameters **including**
    cpu_workers_fraction so downstream code can do `ga.cpu_workers_fraction`.
    """
    population_size: int
    num_generations: int
    mutation_rate_initial: float
    crossover_rate: float
    elitism_count: int
    random_injection_rate: float
    tournament_size: int
    enable_hill_climb: bool
    hill_climb_period: int
    hill_climb_junc: int
    hill_climb_cands: int
    allowed_motif_scale: float
    cpu_workers_fraction: float


@dataclass(frozen=True)
class GlobalConfig:
    """
    Global configuration covering GA knobs and thermodynamic parameters.
    """
    ga: GAConfig
    cpu_workers_fraction: float = 0.9

    # Removed: enforce_span_across_junction

    max_total_adjustment_rounds: int = 2
    division_adjust_attempts_per_junction: int = 4
    division_adjust_step_bp: int = 5

    # Melting temperature method & chemistry
    tm_method: str = "NN"
    Na: float = 50e-3
    K: float = 0.0
    Tris: float = 0.0
    Mg: float = 0.0
    dNTPs: float = 0.0
    dnac1: float = 250e-9
    dnac2: float = 250e-9
    saltcorr: int = 7

    # --- Compatibility helpers: Tm --------------------------------------------
    def tm_params_dict(self) -> Dict[str, Any]:
        """
        Return a dict with Tm-related parameters expected by GA/fitness code.

        Keys are intentionally named to match existing callers:
          - 'tm_method' : str
          - 'Na', 'K', 'Tris', 'Mg', 'dNTPs' : floats (M)
          - 'dnac1', 'dnac2' : strand concentrations (M)
          - 'saltcorr' : int (method switch for salt correction in Biopython)
        """
        return {
            "tm_method": self.tm_method,
            "Na": self.Na,
            "K": self.K,
            "Tris": self.Tris,
            "Mg": self.Mg,
            "dNTPs": self.dNTPs,
            "dnac1": self.dnac1,
            "dnac2": self.dnac2,
            "saltcorr": self.saltcorr,
        }

    @property
    def tm_params(self) -> Dict[str, Any]:
        """Back-compat alias returning the same dict as `tm_params_dict()`."""
        return self.tm_params_dict()

    # --- Compatibility helpers: GA --------------------------------------------
    def ga_params(self) -> GAParamsView:
        """
        Return an attribute-accessible view object with all GA knobs plus
        `cpu_workers_fraction`. This satisfies code paths that do
        `ga = global_cfg.ga_params()` and later access attributes like
        `ga.mutation_rate_initial` or `ga.cpu_workers_fraction`.
        """
        g = self.ga
        return GAParamsView(
            population_size=g.population_size,
            num_generations=g.num_generations,
            mutation_rate_initial=g.mutation_rate_initial,
            crossover_rate=g.crossover_rate,
            elitism_count=g.elitism_count,
            random_injection_rate=g.random_injection_rate,
            tournament_size=g.tournament_size,
            enable_hill_climb=g.enable_hill_climb,
            hill_climb_period=g.hill_climb_period,
            hill_climb_junc=g.hill_climb_junc,
            hill_climb_cands=g.hill_climb_cands,
            allowed_motif_scale=g.allowed_motif_scale,
            cpu_workers_fraction=self.cpu_workers_fraction,
        )

    @property
    def ga_params_dict(self) -> Dict[str, Any]:
        """
        Back-compat alias for callers that prefer a plain dict of GA params.
        Includes `cpu_workers_fraction`.
        """
        gp = self.ga_params()
        return {
            "population_size": gp.population_size,
            "num_generations": gp.num_generations,
            "mutation_rate_initial": gp.mutation_rate_initial,
            "crossover_rate": gp.crossover_rate,
            "elitism_count": gp.elitism_count,
            "random_injection_rate": gp.random_injection_rate,
            "tournament_size": gp.tournament_size,
            "enable_hill_climb": gp.enable_hill_climb,
            "hill_climb_period": gp.hill_climb_period,
            "hill_climb_junc": gp.hill_climb_junc,
            "hill_climb_cands": gp.hill_climb_cands,
            "allowed_motif_scale": gp.allowed_motif_scale,
            "cpu_workers_fraction": gp.cpu_workers_fraction,
        }


# ------------------------------- Loader ---------------------------------------


def load_global_config(path: Optional[Path | str] = None) -> GlobalConfig:
    """
    Load the global configuration from JSON.

    Args:
        path: Optional path to globals.json. If None, uses the file adjacent to
              this module: backend/app/config/globals.json.
    """
    if path is None:
        path = Path(__file__).with_name("globals.json")
    p = Path(path)
    data = json.loads(p.read_text())

    ga_data = data.get("ga", {}) or {}

    ga = GAConfig(
        population_size=int(ga_data.get("population_size", 200)),
        num_generations=int(ga_data.get("num_generations", 5)),
        mutation_rate_initial=float(ga_data.get("mutation_rate_initial", 0.35)),
        crossover_rate=float(ga_data.get("crossover_rate", 0.6)),
        elitism_count=int(ga_data.get("elitism_count", 1)),
        random_injection_rate=float(ga_data.get("random_injection_rate", 0.05)),
        tournament_size=int(ga_data.get("tournament_size", 3)),
        enable_hill_climb=bool(ga_data.get("enable_hill_climb", True)),
        hill_climb_period=int(ga_data.get("hill_climb_period", 5)),
        hill_climb_junc=int(ga_data.get("hill_climb_junc", 2)),
        hill_climb_cands=int(ga_data.get("hill_climb_cands", 10)),
        allowed_motif_scale=float(ga_data.get("allowed_motif_scale", 1.0)),
    )

    tm_block = data.get("tm", {}) or {}
    gc = GlobalConfig(
        ga=ga,
        cpu_workers_fraction=float(data.get("cpu_workers_fraction", 0.9)),
        max_total_adjustment_rounds=int(data.get("max_total_adjustment_rounds", 2)),
        division_adjust_attempts_per_junction=int(data.get("division_adjust_attempts_per_junction", 4)),
        division_adjust_step_bp=int(data.get("division_adjust_step_bp", 5)),
        tm_method=str(data.get("tm_method", "NN")),
        Na=float(tm_block.get("Na", 50e-3)),
        K=float(tm_block.get("K", 0.0)),
        Tris=float(tm_block.get("Tris", 0.0)),
        Mg=float(tm_block.get("Mg", 0.0)),
        dNTPs=float(tm_block.get("dNTPs", 0.0)),
        dnac1=float(tm_block.get("dnac1", 250e-9)),
        dnac2=float(tm_block.get("dnac2", 250e-9)),
        saltcorr=int(tm_block.get("saltcorr", 7)),
    )
    return gc
