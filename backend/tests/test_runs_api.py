# File: backend/tests/test_runs_api.py
# Version: v0.1.0
"""
Simple tests for the Runs API.
"""
from fastapi.testclient import TestClient
from backend.app.main import app
from backend.app.db.session import SessionLocal
from backend.app.services.run_store import record_run_start, record_run_completion

client = TestClient(app)

def test_runs_list_and_detail_roundtrip():
    # create a run
    db = SessionLocal()
    try:
        record_run_start(db, job_id="abcd1234ef56", sequence_len=123, note="Integration test run")
        record_run_completion(db, job_id="abcd1234ef56", report_url="/reports/abcd1234ef56/index.html", exit_code=0)
    finally:
        db.close()

    # list
    r = client.get("/api/runs")
    assert r.status_code == 200
    data = r.json()
    assert data["total"] >= 1
    target = next(item for item in data["items"] if item["job_id"] == "abcd1234ef56")
    assert target["note"] == "Integration test run"

    # detail
    r2 = client.get("/api/runs/abcd1234ef56")
    assert r2.status_code == 200
    d2 = r2.json()
    assert d2["job_id"] == "abcd1234ef56"
    assert d2["status"] in ("completed", "failed")
    assert d2["report_url"] == "/reports/abcd1234ef56/index.html"

def test_delete_run():
    db = SessionLocal()
    try:
        record_run_start(db, job_id="to-delete-0001")
    finally:
        db.close()
    r = client.delete("/api/runs/to-delete-0001")
    assert r.status_code == 204
    r2 = client.get("/api/runs/to-delete-0001")
    assert r2.status_code == 404
