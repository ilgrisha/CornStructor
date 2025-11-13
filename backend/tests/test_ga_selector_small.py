# File: backend/tests/test_ga_selector_small.py
# Version: v0.1.0

"""
Smoke test for GAOverlapSelector using the threaded, C-accelerated hot paths.

This ensures:
- Candidate generation works without span requirement.
- Threaded evaluation (ThreadPool) runs and returns a best chromosome.
- Fitness is assigned and progress log length == num_generations.
"""

from __future__ import annotations

from backend.app.core.assembly.ga_overlap_selector import (
    GAOverlapSelector,
    GAParameters,
)


def test_ga_evolve_small():
    # Simple synthetic sequence with 3 oligos → 2 junctions
    full = "ACGT" * 30  # 120 bp
    # Oligo bodies (non-overlapping bodies; overlaps are chosen around edges)
    positions = [(0, 40), (40, 80), (80, 120)]
    bodies = [full[s:e] for s, e in positions]

    ga = GAParameters(
        population_size=20,
        num_generations=2,
        workers=4,           # exercise threaded evaluation
        batch_chunk_size=8,  # small chunks
    )

    sel = GAOverlapSelector(
        full_sequence=full,
        oligo_positions=positions,
        oligo_seqs=bodies,
        overlap_min_size=8,
        overlap_max_size=16,
        tm_min=None,
        tm_max=None,
        gc_min=20.0,
        gc_max=80.0,
        run_min=1,
        run_max=6,
        disallowed_motifs=["AAAAAA"],
        allowed_motifs={"CG": 0.5},
        tm_method="PRIMER3",
        tm_params={},  # defaults in compute_tm
        ga_params=ga,
    )

    best, log = sel.evolve()
    assert len(log) == ga.num_generations
    assert best.fitness is not None
    assert isinstance(best.overlaps, list)
    assert len(best.overlaps) == len(positions) - 1
