"""Convenience launcher — run from the repo root.

    python scripts/start_studio.py

Then open http://127.0.0.1:5000 in a browser.
"""
import sys
from pathlib import Path

# The script lives in <repo>/scripts/, so repo root is parents[1]
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from bioscreen.studio.app import app, _load_engine

if __name__ == "__main__":
    print("Loading engine (ESM-C 600M + reference pickle)…")
    _load_engine()
    print("Studio running at http://127.0.0.1:5000")
    app.run(host="127.0.0.1", port=5000, debug=False, use_reloader=False)
