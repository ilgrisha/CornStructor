# File: backend/tests/test_run_index.py
# Version: v0.1.0
"""Tests for the run_index HTML generator."""
from __future__ import annotations

from pathlib import Path
import re

from backend.app.core.visualization.run_index import write_run_index

def test_write_run_index_creates_links(tmp_path: Path):
    # Arrange: create fake artifacts
    outdir = tmp_path / 'jobabc123'
    outdir.mkdir(parents=True)
    (outdir / 'seq1_analysis.html').write_text('<html>analysis</html>', encoding='utf-8')
    (outdir / 'seq1_ga_progress.html').write_text('<html>ga</html>', encoding='utf-8')
    (outdir / 'seq1_cluster_reports').mkdir()
    # Act
    idx = write_run_index(outdir, job_id='jobabc123', reports_public_base='/reports')
    # Assert
    assert idx is not None and idx.exists()
    text = idx.read_text(encoding='utf-8')
    assert 'CornStructor results — jobabc123' in text
    assert re.search(r"href='seq1_analysis.html'", text)
    assert re.search(r"href='seq1_ga_progress.html'", text)
    assert re.search(r"href='seq1_cluster_reports/'", text)
