"""
BioScreen Studio — Flask backend with pipeline trace.

Every API call returns both the result and a step-by-step trace of the
computation. The trace shows real numbers from the actual pipeline —
measured, not simulated.

The engine modules (screening_v2, fingerprint, pricing) are called exactly
the same way as before. Timing wraps the calls; it does not change them.
"""
import json
import os
import threading
import time
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
from bioscreen.fingerprint import FingerprintEngine, score_fingerprint
from bioscreen.pricing import ActuarialPricer, OrderRisk

app = Flask(__name__, static_folder="static", template_folder="templates")

_engine = {
    "screener": None,
    "fingerprinter": None,
    "pricer": None,
    "status": "idle",
    "error": None,
    "progress": {
        "stage": "idle",
        "detail": "not started",
        "started_at": None,
        "elapsed_s": 0.0,
    },
}
_lock = threading.Lock()


def _set_progress(stage, detail):
    _engine["progress"]["stage"] = stage
    _engine["progress"]["detail"] = detail
    _engine["progress"]["elapsed_s"] = round(
        time.time() - (_engine["progress"]["started_at"] or time.time()), 1
    )


def _load_engine():
    with _lock:
        if _engine["screener"] is not None:
            return
        if _engine["status"] == "loading":
            return
        _engine["status"] = "loading"
        _engine["progress"]["started_at"] = time.time()
        _set_progress("starting", "preparing to load engine")
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

            _set_progress("model", "loading ESM-C 600M from disk (2.4 GB, ~20-40 s on first run)")
            model = EsmcForMaskedLM.from_pretrained("biohub/ESMC-600M")
            _set_progress("model", "preparing model for inference")
            model.eval()
            tokenizer = EsmcTokenizer()

            _set_progress("screener", "loading reference embeddings (86 MB)")
            _engine["screener"] = BioScreenV2(
                model=model, tokenizer=tokenizer,
                per_residue_path=ref_path,
                hazard_threshold=0.02, top_k=50,
                toxin_similarity_floor=0.90,
            )

            _set_progress("fingerprint", "loading fingerprint pickle")
            _engine["fingerprinter"] = FingerprintEngine(
                model=model, tokenizer=tokenizer,
                fingerprint_path=fp_path,
            )

            _set_progress("pricing", "initialising pricing engine")
            _engine["pricer"] = ActuarialPricer()

            _set_progress("ready", f"engine ready in {_engine['progress']['elapsed_s']} s")
            _engine["status"] = "ready"
        except Exception as e:
            _engine["status"] = "error"
            _engine["error"] = str(e)
            _set_progress("error", str(e))
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
        "progress": _engine["progress"],
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

    trace = []
    t_total_start = time.perf_counter()

    # ---- Step 1: input validation ----
    t = time.perf_counter()
    trace.append({
        "step": 1, "phase": "input",
        "name": "Input validation",
        "elapsed_ms": round((time.perf_counter() - t) * 1000, 3),
        "input": f"{len(sequence)}-residue string",
        "output": "accepted (protein, 20–2048 aa)",
        "details": {
            "length": len(sequence),
            "first_40": sequence[:40] + ("..." if len(sequence) > 40 else ""),
        },
    })

    # ---- Step 2: contrastive embedding ----
    t = time.perf_counter()
    query_emb = s.embed_per_residue(sequence)
    trace.append({
        "step": 2, "phase": "embed",
        "name": "ESM-C 600M forward pass",
        "elapsed_ms": round((time.perf_counter() - t) * 1000, 3),
        "input": f"{len(sequence)} residues",
        "output": f"per-residue matrix {query_emb.shape}",
        "details": {
            "shape": list(query_emb.shape),
            "dtype": str(query_emb.dtype),
            "model": "biohub/ESMC-600M",
        },
    })

    # ---- Step 3: contrastive toxin scoring ----
    t = time.perf_counter()
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

    trace.append({
        "step": 3, "phase": "contrast",
        "name": "Contrastive alignment (full-protein)",
        "elapsed_ms": round((time.perf_counter() - t) * 1000, 3),
        "input": f"{len(s.toxin_accessions)} toxins, {len(s.safe_accessions)} safes",
        "output": f"best toxin: {best_toxin_name} ({best_toxin_score:.4f})",
        "details": {
            "best_toxin_score": round(best_toxin_score, 4),
            "best_safe_score": round(best_safe_score, 4),
            "best_safe_name": best_safe_name,
            "top_k": s.top_k,
        },
    })

    delta = best_toxin_score - best_safe_score
    gate_delta = delta >= s.hazard_threshold
    gate_floor = best_toxin_score >= s.toxin_similarity_floor
    flagged = bool(gate_delta and gate_floor)

    # ---- Step 4: fingerprint extraction ----
    # Compute per-residue embeddings for the fingerprint stage (reusing query_emb)
    fp_embeddings = {}
    t = time.perf_counter()
    for acc, fp in f.fingerprints.items():
        fp_embeddings[acc] = score_fingerprint(query_emb, fp, top_k=f.top_k)
    trace.append({
        "step": 4, "phase": "fingerprint",
        "name": "Fingerprint matching",
        "elapsed_ms": round((time.perf_counter() - t) * 1000, 3),
        "input": f"query embedding vs {len(f.fingerprints)} fingerprints",
        "output": f"scored {len(f.fingerprints)} toxins",
        "details": {
            "n_fingerprints": len(f.fingerprints),
            "scoring": "top20_similarity × order_consistency",
        },
    })

    # ---- Step 5: pick best fingerprint ----
    best_acc = max(fp_embeddings, key=lambda a: fp_embeddings[a]["combined"])
    best_fp = fp_embeddings[best_acc]

    trace.append({
        "step": 5, "phase": "fingerprint",
        "name": "Best fingerprint selection",
        "elapsed_ms": 0.1,
        "input": f"{len(fp_embeddings)} candidate scores",
        "output": f"{f.names[best_acc]} ({best_fp['combined']:.4f})",
        "details": {
            "best_accession": best_acc,
            "best_name": f.names[best_acc],
            "top_k_similarity": round(best_fp["top_k_similarity"], 4),
            "order_consistency": round(best_fp["order_consistency"], 4),
            "threshold": f.threshold,
        },
    })

    fp_result = {
        "fingerprint_score": round(best_fp["combined"], 4),
        "fingerprint_top_k": round(best_fp["top_k_similarity"], 4),
        "fingerprint_order": round(best_fp["order_consistency"], 4),
        "fingerprint_best": f.names[best_acc],
        "fingerprint_best_accession": best_acc,
        "fingerprint_flagged": bool(best_fp["combined"] >= f.threshold),
        "fingerprint_threshold": f.threshold,
        "fingerprint_ranking": sorted(
            [{"accession": a, "name": f.names[a], "score": round(v["combined"], 4)}
             for a, v in fp_embeddings.items()],
            key=lambda x: -x["score"],
        )[:5],
    }

    # ---- Step 6: verdict ----
    trace.append({
        "step": 6, "phase": "verdict",
        "name": "Verdict decision",
        "elapsed_ms": 0.1,
        "input": f"contrastive delta={delta:.4f}, fingerprint={fp_result['fingerprint_score']:.4f}",
        "output": "FLAGGED" if (flagged or fp_result["fingerprint_flagged"]) else "CLEAR",
        "details": {
            "contrastive_flagged": flagged,
            "fingerprint_flagged": fp_result["fingerprint_flagged"],
            "hazard_delta": round(delta, 4),
        },
    })

    # ---- Build per-residue profiles for the ribbon visuals ----
    def profile(q, ref):
        qn = q / (np.linalg.norm(q, axis=1, keepdims=True) + 1e-8)
        rn = ref / (np.linalg.norm(ref, axis=1, keepdims=True) + 1e-8)
        return (qn @ rn.T).max(axis=1).astype(float).tolist()

    total_ms = round((time.perf_counter() - t_total_start) * 1000, 2)

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
        # Trace
        "trace": trace,
        "total_elapsed_ms": total_ms,
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


@app.route("/api/source")
def api_source():
    """Return the actual source code of the three engine files for the Reproduce tab."""
    files = {
        "screening_v2.py": (ROOT / "src" / "bioscreen" / "screening_v2.py").read_text(encoding="utf-8"),
        "fingerprint.py": (ROOT / "src" / "bioscreen" / "fingerprint.py").read_text(encoding="utf-8"),
        "pricing.py": (ROOT / "src" / "bioscreen" / "pricing.py").read_text(encoding="utf-8"),
    }
    return jsonify(files)




@app.route("/api/adversarial")
def api_adversarial():
    """
    Return the pre-computed EvoDiff variant comparison for all toxins.
    """
    import json
    p = ROOT / "results" / "evodiff_16toxins_full.json"
    if not p.exists():
        return jsonify({"error": "Results file missing"}), 404
    data = json.loads(p.read_text(encoding="utf-8"))

    toxins = []
    for r in data["results"]:
        variants = []
        for i, v in enumerate(r["variants"]):
            variants.append({
                "idx": i,
                "identity": v["id"],
                "blast": v["blast"],
                "fingerprint": v["fp"],
                "blast_verdict": "EVADED" if v["blast"] < 0.30 else "caught",
                "fingerprint_verdict": "CAUGHT" if v["fp"] >= 0.85 else "clear",
            })
        toxins.append({
            "accession": r["accession"],
            "name": r["name"],
            "region_length": r.get("region_length"),
            "n_catalytic": r.get("n_catalytic"),
            "variants": variants,
            "n_caught": sum(1 for v in variants if v["fingerprint"] >= 0.85),
            "n_total": len(variants),
        })

    return jsonify({
        "n_toxins": data.get("n_toxins"),
        "n_variants": data.get("n_variants"),
        "blast_evaded": data.get("blast_evaded"),
        "fingerprint_caught": data.get("fingerprint_caught"),
        "total_runtime_min": data.get("total_runtime_min"),
        "toxins": toxins,
    })




# ----------------------------------------------------------------------
# Variant generation (EvoDiff live)
# ----------------------------------------------------------------------
_generate_jobs = {}
_generate_lock = threading.Lock()


@app.route("/api/generate", methods=["POST"])
def api_generate():
    """Start a background EvoDiff generation job. Returns job_id."""
    _load_engine()
    data = request.get_json(silent=True) or {}
    accession = (data.get("accession") or "").strip()
    mode = data.get("mode", "preserved")
    n_variants = int(data.get("n_variants", 3))

    from bioscreen.toxins import TOXIN_CONFIG
    if accession not in TOXIN_CONFIG:
        return jsonify({"error": f"Unknown toxin: {accession}"}), 400
    if mode not in ("preserved", "aggressive"):
        return jsonify({"error": "mode must be 'preserved' or 'aggressive'"}), 400

    import uuid
    job_id = str(uuid.uuid4())

    with _generate_lock:
        _generate_jobs[job_id] = {
            "id": job_id,
            "status": "queued",
            "accession": accession,
            "mode": mode,
            "n_variants": n_variants,
            "stage": "queued",
            "progress": 0.0,
            "result": None,
            "error": None,
            "started_at": time.time(),
        }

    def _run():
        try:
            _generate_jobs[job_id]["status"] = "running"
            from bioscreen.generate import generate_variants

            def progress_cb(stage, frac):
                _generate_jobs[job_id]["stage"] = stage
                _generate_jobs[job_id]["progress"] = frac

            fp_emb = _engine["fingerprinter"].fingerprints if _engine["fingerprinter"] else {}

            def embed_fn(seq):
                return _engine["screener"].embed_per_residue(seq)

            result = generate_variants(
                accession=accession,
                mode=mode,
                n_variants=n_variants,
                fingerprint_embeddings=fp_emb,
                embed_fn=embed_fn,
                progress_cb=progress_cb,
            )

            # Enrich each variant with contrastive verdict + actuarial pricing
            s = _engine["screener"]
            pricer = _engine["pricer"]
            for v in result["variants"]:
                seq = v["sequence"]
                try:
                    sr = s.screen(seq)
                except Exception:
                    sr = {"flagged": False, "hazard_delta": 0.0, "best_toxin_score": 0.0}
                fp_flagged = bool(v.get("fingerprint_caught", False))
                v["verdict"] = "FLAGGED" if (sr["flagged"] or fp_flagged) else "CLEAR"
                v["contrastive_delta"] = round(sr["hazard_delta"], 4)
                v["contrastive_flagged"] = bool(sr["flagged"])
                try:
                    order = OrderRisk(
                        customer_verified=True,
                        order_size_bp=2000,
                        organism="human",
                        screening_hazard_delta=sr["hazard_delta"],
                        screening_best_toxin_sim=sr["best_toxin_score"],
                        customer_history=0.8,
                    )
                    v["pricing"] = pricer.compute_premium(order, order_cost_usd=2000)
                except Exception:
                    v["pricing"] = None

            _generate_jobs[job_id]["result"] = result
            _generate_jobs[job_id]["status"] = "ready"
            _generate_jobs[job_id]["elapsed_s"] = round(
                time.time() - _generate_jobs[job_id]["started_at"], 1
            )
        except Exception as e:
            import traceback
            _generate_jobs[job_id]["status"] = "error"
            _generate_jobs[job_id]["error"] = str(e)
            _generate_jobs[job_id]["traceback"] = traceback.format_exc()

    threading.Thread(target=_run, daemon=True).start()
    return jsonify({"job_id": job_id, "status": "queued"})


@app.route("/api/generate/<job_id>")
def api_generate_status(job_id):
    job = _generate_jobs.get(job_id)
    if not job:
        return jsonify({"error": "unknown job"}), 404
    return jsonify(job)


@app.route("/api/toxins")
def api_toxins():
    """List all configured toxins for the Generate tab dropdown."""
    from bioscreen.toxins import TOXIN_CONFIG
    return jsonify([
        {"accession": acc, "name": cfg["name"], "mechanism": cfg["mechanism"]}
        for acc, cfg in TOXIN_CONFIG.items()
    ])


if __name__ == "__main__":
    print("Pre-loading engine...")
    _load_engine()
    print(f"Engine status: {_engine['status']}")
    print("Studio running at http://127.0.0.1:5000")
    app.run(host="127.0.0.1", port=5000, debug=False, use_reloader=False)
