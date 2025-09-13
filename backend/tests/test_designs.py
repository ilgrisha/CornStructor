# File: backend/tests/test_designs.py
# Version: v0.1.0
"""Smoke tests for Design persistence endpoints."""
from fastapi.testclient import TestClient

from backend.app.main import app

client = TestClient(app)

def test_start_persists_design_and_fetch_by_run(monkeypatch):
    # Start a design
    payload = {"sequence": "ACGT" * 20, "params": {"gcWin": 50}}
    res = client.post("/api/design/start", json=payload)
    assert res.status_code == 200
    job_id = res.json()["jobId"]
    # Try to fetch design by run
    res2 = client.get(f"/api/designs/by-run/{job_id}")
    assert res2.status_code == 200
    data = res2.json()
    assert data["job_id"] == job_id
    assert data["sequence_len"] == 80
    assert "sequence" in data
