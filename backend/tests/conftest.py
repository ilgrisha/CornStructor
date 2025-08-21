# File: backend/tests/conftest.py
# Version: v0.1.0
"""
Test bootstrap: ensure project root is on sys.path so 'backend.*' imports work.

This avoids requiring editable installs or extra plugins. It keeps tests hermetic
to the repo layout (works in CI and locally).
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]  # repo root (../.. from this file)
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
