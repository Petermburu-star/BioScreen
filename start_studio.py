"""Convenience launcher — run from the repo root."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "src"))

from bioscreen.studio.app import app, _load_engine

if __name__ == "__main__":
    print("Loading engine (ESM-C 600M + reference pickle)...")
    _load_engine()
    print("Studio running at http://127.0.0.1:5000")
    app.run(host="127.0.0.1", port=5000, debug=False, use_reloader=False)
