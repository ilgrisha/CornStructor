# File: backend/tests/test_secondary_structure.py
# Version: v0.1.0
"""
Tests for the secondary structure analysis API.

These tests use httpx.AsyncClient to hit the FastAPI app in-memory.
"""

import pytest
from httpx import AsyncClient
from backend.app.main import app

pytestmark = pytest.mark.asyncio


async def test_stems_basic_ok():
    async with AsyncClient(app=app, base_url="http://test") as ac:
        payload = {
            "sequence": "ACGT" * 30,  # 120bp
            "min_stem_len": 4,
            "merge_max_gap": 2,
        }
        resp = await ac.post("/api/v1/analysis/stems", json=payload)
        assert resp.status_code == 200
        data = resp.json()
        assert "regions" in data
        assert data["length"] == 120
        # Regions may be empty depending on MFE, but schema should be correct
        assert isinstance(data["regions"], list)
        for r in data["regions"]:
            assert r["kind"] == "stems"
            assert 0 <= r["start"] <= r["end"] <= data["length"]


async def test_stems_validation():
    async with AsyncClient(app=app, base_url="http://test") as ac:
        bad_payload = {"sequence": "", "min_stem_len": 0, "merge_max_gap": -1}
        resp = await ac.post("/api/v1/analysis/stems", json=bad_payload)
        # Pydantic validation error
        assert resp.status_code == 422
