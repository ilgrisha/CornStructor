# File: backend/app/config/config_global.py
# Version: v0.3.2

"""
Global configuration loader (attribute-access "view" objects).

v0.3.2
- Back-compat shims for assembler:
  * `tm_params_dict()` → returns the thermo params dict (mM/nM)
  * `tm_method_str()`  → returns the selected Tm method string
- Keep legacy `load_global_config(path)` import for the CLI.

v0.3.1
- Added legacy shim `load_global_config(path)` so older CLIs can import it.

v0.3.0
- Introduced GA perf knobs from JSON: `ga.workers` and `ga.batch_chunk_size`.
- Back-compat derivation of workers from `cpu_workers_fraction` when missing or 0.
"""

from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass
from typing import Any, Dict, Optional


# --------------------------
# GA params (attribute view)
# --------------------------

@dataclass(frozen=True)
class GAParamsView:
    # Core GA knobs
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

    # Perf knobs
    workers: int              # >=0; 0 means "auto" (selector decides)
    batch_chunk_size: int     # chunk size for batch helpers

    @staticmethod
    def _cpu_count() -> int:
        try:
            return os.cpu_count() or 1
        except Exception:
            return 1

    @classmethod
    def from_dict(
        cls,
        d: Dict[str, Any],
        *,
        cpu_workers_fraction: Optional[float] = None,
    ) -> "GAParamsView":
        """
        Build a GA params view from a dict. If 'workers' is 0 or missing,
        derive it from `cpu_workers_fraction` (if provided).
        """
        workers_raw = int(d.get("workers", 0))
        chunk = int(d.get("batch_chunk_size", 128))

        if workers_raw <= 0 and cpu_workers_fraction is not None:
            frac = max(0.05, min(1.0, float(cpu_workers_fraction)))
            derived = max(1, math.floor(cls._cpu_count() * frac))
            workers = derived
        else:
            workers = max(0, workers_raw)  # keep 0 for selector auto

        return cls(
            population_size=int(d["population_size"]),
            num_generations=int(d["num_generations"]),
            mutation_rate_initial=float(d["mutation_rate_initial"]),
            crossover_rate=float(d["crossover_rate"]),
            elitism_count=int(d["elitism_count"]),
            random_injection_rate=float(d["random_injection_rate"]),
            tournament_size=int(d["tournament_size"]),
            enable_hill_climb=bool(d["enable_hill_climb"]),
            hill_climb_period=int(d["hill_climb_period"]),
            hill_climb_junc=int(d["hill_climb_junc"]),
            hill_climb_cands=int(d["hill_climb_cands"]),
            allowed_motif_scale=float(d.get("allowed_motif_scale", 1.0)),
            workers=workers,
            batch_chunk_size=max(1, chunk),
        )


# --------------------------
# Top-level global config
# --------------------------

@dataclass(frozen=True)
class GlobalConfig:
    _ga_raw: Dict[str, Any]
    _cpu_workers_fraction: float

    # Assembler knobs
    max_total_adjustment_rounds: int
    division_adjust_attempts_per_junction: int
    division_adjust_step_bp: int

    # Thermo method + parameters
    tm_method: str
    tm: Dict[str, float]

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "GlobalConfig":
        ga_raw = dict(data.get("ga", {}))
        cpu_frac = float(data.get("cpu_workers_fraction", 0.9))
        return cls(
            _ga_raw=ga_raw,
            _cpu_workers_fraction=cpu_frac,
            max_total_adjustment_rounds=int(data.get("max_total_adjustment_rounds", 2)),
            division_adjust_attempts_per_junction=int(
                data.get("division_adjust_attempts_per_junction", 4)
            ),
            division_adjust_step_bp=int(data.get("division_adjust_step_bp", 5)),
            tm_method=str(data.get("tm_method", "PRIMER3")),
            tm={k: float(v) for k, v in dict(data.get("tm", {})).items()},
        )

    @classmethod
    def from_json_file(cls, path: str) -> "GlobalConfig":
        with open(path, "r", encoding="utf-8") as f:
            return cls.from_dict(json.load(f))

    # ---- Views used by the assembler / GA ----

    def ga_params(self) -> GAParamsView:
        """
        Return a resolved GA params view:
          - workers >0 honored; 0/missing derived from cpu_workers_fraction
          - batch_chunk_size defaulted to 128 if missing
        """
        return GAParamsView.from_dict(self._ga_raw, cpu_workers_fraction=self._cpu_workers_fraction)

    @property
    def cpu_workers_fraction(self) -> float:
        return self._cpu_workers_fraction

    # ---- Back-compat shims expected by assembler ----

    def tm_params_dict(self) -> Dict[str, float]:
        """Legacy shim: return the melting temperature parameters dict."""
        return dict(self.tm)

    def tm_method_str(self) -> str:
        """Legacy shim: return the selected TM method string."""
        return str(self.tm_method)


# --------------------------
# Legacy shim for CLI import
# --------------------------

def load_global_config(path: str) -> GlobalConfig:
    """
    Legacy helper used by `assemble_cli.py`:
    Equivalent to `GlobalConfig.from_json_file(path)`.
    """
    return GlobalConfig.from_json_file(path)


__all__ = [
    "GAParamsView",
    "GlobalConfig",
    "load_global_config",
]
