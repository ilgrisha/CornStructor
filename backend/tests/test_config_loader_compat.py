# File: backend/tests/test_config_loader_compat.py
# Version: v0.1.0

"""
Ensure legacy import path used by assemble_cli keeps working.
"""

from backend.app.config.config_global import load_global_config, GlobalConfig


def test_load_global_config_shim(tmp_path):
    p = tmp_path / "g.json"
    p.write_text(
        """{
          "ga": {
            "population_size": 10,
            "num_generations": 1,
            "mutation_rate_initial": 0.3,
            "crossover_rate": 0.6,
            "elitism_count": 1,
            "random_injection_rate": 0.05,
            "tournament_size": 3,
            "enable_hill_climb": false,
            "hill_climb_period": 5,
            "hill_climb_junc": 1,
            "hill_climb_cands": 3,
            "allowed_motif_scale": 1.0
          },
          "cpu_workers_fraction": 0.9,
          "tm_method": "NN",
          "tm": {}
        }""",
        encoding="utf-8",
    )
    cfg = load_global_config(str(p))
    assert isinstance(cfg, GlobalConfig)
    ga = cfg.ga_params()
    assert ga.population_size == 10
