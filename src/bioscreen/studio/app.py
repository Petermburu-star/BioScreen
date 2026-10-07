"""
BioScreen Studio — Flask backend.

Loads ESM-C, the reference embeddings, and the fingerprint pickle once,
then serves the Studio UI. Every result is computed live.
"""
import os
import threading
from pathlib import Path

import numpy as np
from flask import Flask, jsonify, render_template, request

try:
    from dotenv import load_dotenv
    load_dotenv(Path(__file__).resolve().parents[3] / ".env")
except ImportError:
    pass

ROOT = Path(__file__).resolve().parents[3]

from bioscreen.screening_v2 import BioScreenV2, _local_alignment
from bioscreen.fingerprint import FingerprintEngine
from bioscreen.pricing import ActuarialPricer, OrderRisk

app = Flask(__name__, static_folder="static", template_folder="templates")

_engine = {
    "screener": None,
    "fingerprinter": None,
    "pricer": None,
    "status": "idle",
    "error": None,
}
_lock = threading.Lock()


def _load_engine():
    with _lock:
        if _engine["screener"] is not None:
            return
        if _engine["status"] == "loading":
            return
        _engine["status"] = "loading"
        try:
            import torch
            from esm.models.esmc import EsmcForMaskedLM, EsmcTokenizer

            ref_path = Path(os.environ.get(
                "BIOSCREEN_REF_PATH",
                ROOT / "data" / "processed" / "per_residue_embeddings.pkl",
            ))
            fp_path = Path(os.environ.get(
                "BIOSCREEN_FP_PATH",
                ROOT / "data" / "processed" / "toxin_fingerprints.pkl",
            ))

            if not ref_path.exists():
                raise FileNotFoundError(
                    f"Reference pickle missing: {ref_path}\n"
                    f"Run: python scripts/build_references.py"
                )
            if not fp_path.exists():
                raise FileNotFoundError(
                    f"Fingerprint pickle missing: {fp_path}\n"
                    f"Run: python scripts/build_fingerprints.py"
                )

            model = EsmcForMaskedLM.from_pretrained("biohub/ESMC-600M")
            model.eval()
            tokenizer = EsmcTokenizer()

            _engine["screener"] = BioScreenV2(
                model=model,
                tokenizer=tokenizer,
                per_residue_path=ref_path,
                hazard_threshold=0.02,
                top_k=50,
                toxin_similarity_floor=0.90,
            )
            _engine["fingerprinter"] = FingerprintEngine(
                model=model,
                tokenizer=tokenizer,
                fingerprint_path=fp_path,
            )
            _engine["pricer"] = ActuarialPricer()
            _engine["status"] = "ready"
        except Exception as e:
            _engine["status"] = "error"
            _engine["error"] = str(e)
            raise


# ----------------------------------------------------------------------
# Routes
# ----------------------------------------------------------------------

@app.route("/")
def index():
    return render_template("studio.html")


@app.route("/api/status")
def api_status():
    s = _engine["screener"]
    f = _engine["fingerprinter"]
    return jsonify({
        "status": _engine["status"],
        "error": _engine["error"],
        "toxins": len(s.toxin_accessions) if s else 0,
        "safes": len(s.safe_accessions) if s else 0,
        "fingerprints": len(f.fingerprints) if f else 0,
    })


@app.route("/api/warmup", methods=["POST"])
def api_warmup():
    _load_engine()
    return jsonify({"status": _engine["status"]})


@app.route("/api/screen", methods=["POST"])
def api_screen():
    _load_engine()
    data = request.get_json(silent=True) or {}
    sequence = (data.get("sequence") or "").strip().upper()

    if len(sequence) < 20:
        return jsonify({"error": "Sequence must be at least 20 amino acids"}), 400
    if len(sequence) > 2048:
        sequence = sequence[:2048]

    s = _engine["screener"]
    f = _engine["fingerprinter"]

    # --- Contrastive screening (existing method) ---
    query_emb = s.embed_per_residue(sequence)

    toxin_scores = np.array([
        _local_alignment(query_emb, s.toxin_embeddings[a], top_k=s.top_k)
        for a in s.toxin_accessions
    ])
    ti = int(np.argmax(toxin_scores))
    best_toxin_acc = s.toxin_accessions[ti]
    best_toxin_name = s.toxin_names[ti]
    best_toxin_score = float(toxin_scores[ti])

    safe_scores = np.array([
        _local_alignment(query_emb, s.safe_embeddings[a], top_k=s.top_k)
        for a in s.safe_accessions
    ])
    si = int(np.argmax(safe_scores))
    best_safe_acc = s.safe_accessions[si]
    best_safe_name = s.safe_names[si]
    best_safe_score = float(safe_scores[si])

    delta = best_toxin_score - best_safe_score
    gate_delta = delta >= s.hazard_threshold
    gate_floor = best_toxin_score >= s.toxin_similarity_floor
    flagged = bool(gate_delta and gate_floor)

    # --- Fingerprint screening (new method) ---
    fp_result = f.screen(sequence)

    # --- Per-residue profiles for ribbons ---
    def profile(q, ref):
        qn = q / (np.linalg.norm(q, axis=1, keepdims=True) + 1e-8)
        rn = ref / (np.linalg.norm(ref, axis=1, keepdims=True) + 1e-8)
        return (qn @ rn.T).max(axis=1).astype(float).tolist()

    return jsonify({
        "sequence_length": len(sequence),
        # Contrastive
        "best_toxin_score": round(best_toxin_score, 4),
        "best_toxin_name": best_toxin_name,
        "best_toxin_acc": best_toxin_acc,
        "best_safe_score": round(best_safe_score, 4),
        "best_safe_name": best_safe_name,
        "best_safe_acc": best_safe_acc,
        "hazard_delta": round(delta, 4),
        "flagged": flagged,
        "gate_delta": bool(gate_delta),
        "gate_floor": bool(gate_floor),
        "threshold_delta": s.hazard_threshold,
        "threshold_floor": s.toxin_similarity_floor,
        "toxin_profile": profile(query_emb, s.toxin_embeddings[best_toxin_acc]),
        "safe_profile": profile(query_emb, s.safe_embeddings[best_safe_acc]),
        # Fingerprint
        **fp_result,
    })


@app.route("/api/price", methods=["POST"])
def api_price():
    _load_engine()
    d = request.get_json(silent=True) or {}
    order = OrderRisk(
        customer_verified=bool(d.get("verified", False)),
        order_size_bp=int(d.get("size_bp", 2000)),
        organism=str(d.get("organism", "human")),
        screening_hazard_delta=float(d.get("hazard_delta", 0.0)),
        screening_best_toxin_sim=float(d.get("toxin_sim", 0.0)),
        customer_history=float(d.get("history", 0.5)),
    )
    result = _engine["pricer"].compute_premium(
        order, order_cost_usd=float(d.get("order_cost", 1000.0))
    )
    return jsonify(result)


@app.route("/api/example/<name>")
def api_example(name):
    import requests
    accessions = {
        "ricin": "P02879",
        "botulinum_a": "P0DPI1",
        "diphtheria": "P00588",
        "gfp": "P42212",
        "insulin": "P01308",
        "lysozyme": "P00698",
    }
    acc = accessions.get(name.lower())
    if not acc:
        return jsonify({"error": "Unknown example"}), 404
    r = requests.get(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", timeout=20)
    if r.status_code != 200:
        return jsonify({"error": f"UniProt returned {r.status_code}"}), 502
    lines = r.text.strip().split("\n")
    return jsonify({"header": lines[0][1:], "sequence": "".join(lines[1:]), "accession": acc})


if __name__ == "__main__":
    print("Pre-loading engine...")
    _load_engine()
    print(f"Engine status: {_engine['status']}")
    print("Studio running at http://127.0.0.1:5000")
    app.run(host="127.0.0.1", port=5000, debug=False, use_reloader=False)
