# File: backend/tests/test_api_health.py
# Version: v0.1.0
"""
Basic smoke test for health endpoint.
"""
from fastapi.testclient import TestClient
from backend.app.main import app


def test_health():
    client = TestClient(app)
    r = client.get("/api/health")
    assert r.status_code == 200
    assert r.json() == {"status": "ok"}
