# File: backend/tests/test_ga_selector_backcompat.py
# Version: v0.1.0
"""
Backwards-compat test: GAOverlapSelector accepts GA-like objects lacking new fields.
"""

from backend.app.core.assembly.ga_overlap_selector import GAOverlapSelector

class GAViewLike:
    # Intentionally missing 'workers' and 'batch_chunk_size'
    population_size = 12
    num_generations = 1
    mutation_rate_initial = 0.2
    crossover_rate = 0.5
    elitism_count = 1
    random_injection_rate = 0.1
    tournament_size = 3
    enable_hill_climb = False
    hill_climb_period = 5
    hill_climb_junc = 1
    hill_climb_cands = 5
    allowed_motif_scale = 1.0

def test_backcompat_ga_params():
    full = "ACGT" * 25
    positions = [(0, 30), (30, 50)]
    bodies = [full[s:e] for s, e in positions]

    sel = GAOverlapSelector(
        full_sequence=full,
        oligo_positions=positions,
        oligo_seqs=bodies,
        overlap_min_size=8,
        overlap_max_size=12,
        ga_params=GAViewLike(),  # ← pass a view-like object (no workers attr)
    )
    best, log = sel.evolve()
    assert len(log) == 1
    assert best.fitness is not None
