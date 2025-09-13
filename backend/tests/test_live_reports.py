# File: backend/tests/test_live_reports.py
# Version: v0.1.0
"""Smoke tests for live (on-the-fly) report endpoints."""
from fastapi.testclient import TestClient
import json

from backend.app.main import app

client = TestClient(app)


def test_live_ga_progress_endpoints(monkeypatch):
    # Start a design so pipeline runs and persists GA progress
    seq = "ACGT" * 50
    res = client.post("/api/design/start", json={"sequence": seq, "params": {"gcWin": 50}})
    assert res.status_code == 200
    job_id = res.json()["jobId"]

    # Drain logs to completion to allow pipeline to finish in test env
    with client.stream("GET", f"/api/design/{job_id}/logs") as s:
        for line in s.iter_text():
            if "Log stream completed" in line:
                break

    # GA progress JSON should be available
    res2 = client.get(f"/api/live/by-run/{job_id}/ga_progress.json")
    assert res2.status_code == 200
    data = json.loads(res2.text)
    assert isinstance(data, dict) or isinstance(data, list)

    # HTML should render
    res3 = client.get(f"/api/live/by-run/{job_id}/ga_progress.html")
    assert res3.status_code == 200
    assert "<html" in res3.text.lower() or "<!doctype" in res3.text.lower()

    # Index should include links
    res4 = client.get(f"/api/live/by-run/{job_id}/index.html")
    assert res4.status_code == 200
    assert "GA progress" in res4.text
