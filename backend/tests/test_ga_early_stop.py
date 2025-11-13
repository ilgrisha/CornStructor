# File: backend/tests/test_ga_early_stop.py
# Version: v0.1.0

"""
Early-stop smoke test:
Force a nearly static GA (no crossover/mutation/immigration) so the best fitness
doesn't improve; verify that early-stop triggers before the requested number of
generations and that the run summary reflects it.
"""

from backend.app.core.assembly.ga_overlap_selector import GAOverlapSelector, GAParameters

def test_early_stop_triggers_and_summary():
    # Simple two-oligo case with plenty of candidates but no operators to change population
    full = "ACGT" * 50
    positions = [(0, 80), (80, 200)]
    bodies = [full[s:e] for s, e in positions]

    params = GAParameters(
        population_size=20,
        num_generations=50,          # large cap – we expect to stop much earlier
        mutation_rate_initial=0.0,   # freeze mutation
        crossover_rate=0.0,          # disable crossover
        elitism_count=1,
        random_injection_rate=0.0,   # no immigrants
        tournament_size=3,
        enable_hill_climb=False,

        # Stagnation / early-stop knobs
        stagnation_patience=2,
        early_stop_multiplier=1.0,   # stop after 2 gens without improvement
        enable_early_stop=True,

        # Keep defaults for anti-stagnation boosts off (they won't fire with mutation=0 anyway)
        workers=0,
        batch_chunk_size=64,
    )

    sel = GAOverlapSelector(
        full_sequence=full,
        oligo_positions=positions,
        oligo_seqs=bodies,
        overlap_min_size=8,
        overlap_max_size=12,
        ga_params=params,
        tm_method="Wallace",  # make Tm trivial for test speed
    )

    best, log = sel.evolve()

    # Expect we did not run all 50 generations due to early-stop
    assert len(log) < 50
    assert sel.last_run_summary is not None
    assert sel.last_run_summary.early_stopped is True
    # Should have a valid best fitness
    assert isinstance(sel.last_run_summary.best_fitness, float)
