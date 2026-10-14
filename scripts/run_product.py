"""BioScreen product UI — Flask wrapper around kernel.screen_sequence()."""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]   # BioScreen/ — fixes the launcher bug
sys.path.insert(0, str(ROOT / "src"))
import json
from flask import Flask, request, jsonify, render_template_string
from bioscreen.kernel import screen_sequence, load_db

BANK = json.loads((ROOT / "data/product_db/demo_bank.json").read_text())
HTML = (ROOT / "src/bioscreen/product_ui/index.html").read_text(encoding="utf-8")
app = Flask(__name__)

@app.get("/")
def home(): return render_template_string(HTML)

@app.get("/api/demo")
def demos():
    return jsonify([{k: b[k] for k in ("id", "label", "expected")} for b in BANK])

@app.get("/api/meta")
def meta():
    m = load_db()["meta"]
    return jsonify({"db_version": m["db_version"], "n": m["n_entries"]})

@app.post("/api/screen")
def screen():
    d = request.get_json(force=True, silent=True) or {}
    if d.get("demo_id"):
        b = next((x for x in BANK if x["id"] == d["demo_id"]), None)
        if b is None:
            return jsonify({"error": "unknown demo id"}), 400
        seq, label = b["sequence"], b["label"]
    else:
        seq, label = d.get("sequence", ""), None
    try:
        v = screen_sequence(seq)
    except ValueError as e:
        return jsonify({"error": str(e)}), 400
    v["demo_label"] = label
    return jsonify(v)

if __name__ == "__main__":
    import threading, webbrowser
    threading.Timer(1.5, lambda: webbrowser.open("http://127.0.0.1:5055")).start()
    app.run(host="127.0.0.1", port=5055, debug=False)
