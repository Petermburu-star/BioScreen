"""
Re-run the ricin benchmark end-to-end.

Loads ESM-C 600M, loads the fingerprint pickle, generates 3 EvoDiff
variants of the ricin A-chain region in preserved mode, and reports:

    identity  — sequence identity to the native region
    BLAST     — sequence-similarity score (the incumbent's signal)
    FP        — ordered catalytic fingerprint score

EvoDiff is stochastic. Re-runs produce new variants with similar
qualitative behaviour, not identical scores. One recorded run is in
results/rerun_ricin_benchmark.json.

Run from repo root:

    python scripts/rerun_ricin_benchmark.py

First run downloads ESM-C weights (~2.4 GB) from Hugging Face.
"""
import json
import pickle
import sys
import time
from pathlib import Path

# --- repo root ---
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "external" / "evodiff"))

import torch
from dotenv import load_dotenv

load_dotenv(ROOT / ".env")

from esm.models.esmc import EsmcForMaskedLM, EsmcTokenizer
from bioscreen.generate import generate_variants


ACCESSION = "P02879"   # ricin
MODE = "preserved"
N_VARIANTS = 3


def embed_isolated(seq, model, tokenizer):
    inputs = tokenizer([seq], return_tensors="pt", padding=True)
    with torch.inference_mode():
        out = model(**inputs, output_hidden_states=True)
    return out.hidden_states[-1][0, 1:-1].cpu().numpy()


def fmt_fp(v):
    return f"{v:>9.4f}" if isinstance(v, (int, float)) else f"{'—':>9}"


def main():
    print("Loading ESM-C 600M...")
    model = EsmcForMaskedLM.from_pretrained("biohub/ESMC-600M")
    model.eval()
    tokenizer = EsmcTokenizer()

    fp_path = ROOT / "data" / "processed" / "toxin_fingerprints.pkl"
    if not fp_path.exists():
        raise SystemExit(f"Missing {fp_path}. Run scripts/build_fingerprints.py first.")
    with open(fp_path, "rb") as f:
        fp_data = pickle.load(f)
    fingerprints = fp_data["fingerprints"]
    print(f"  fingerprints : {len(fingerprints)} toxins loaded")

    def embed(seq):
        return embed_isolated(seq, model, tokenizer)

    print(f"\nGenerating {N_VARIANTS} {MODE}-mode variants of {ACCESSION} (ricin)...")
    t0 = time.time()
    result = generate_variants(
        accession=ACCESSION,
        mode=MODE,
        n_variants=N_VARIANTS,
        fingerprint_embeddings=fingerprints,
        embed_fn=embed,
    )
    elapsed = time.time() - t0

    print(f"\nDone in {elapsed:.0f}s")
    print(f"  source        : {result.get('source')}")
    print(f"  region_length : {result.get('region_length')}")
    print(f"  n_catalytic   : {result.get('n_catalytic')}")

    variants = result.get("variants", [])
    print()
    print("=" * 80)
    print(f"  {'v':<4} {'identity':>10} {'BLAST':>9} {'FP':>9}  "
          f"{'BLAST verdict':>14}  {'FP verdict':>12}")
    print("=" * 80)
    for v in variants:
        bl_v = "EVADED" if v.get("blast_evaded") else "caught"
        fp_val = v.get("fingerprint")
        fp_v = ("CAUGHT" if v.get("fingerprint_caught") else "clear") \
               if fp_val is not None else "—"
        ident = v.get("identity", float("nan"))
        blast = v.get("blast", float("nan"))
        print(f"  v{v.get('idx','?'):<3} {ident*100:>9.2f}% {blast:>9.4f} "
              f"{fmt_fp(fp_val)}  {bl_v:>14}  {fp_v:>12}")
    print("=" * 80)

    caught = sum(1 for v in variants if v.get("fingerprint_caught"))
    evaded = sum(1 for v in variants if v.get("blast_evaded"))
    print(f"\n  BLAST evaded      : {evaded} / {len(variants)}")
    print(f"  Fingerprint caught: {caught} / {len(variants)}")

    out = ROOT / "results" / "rerun_ricin_benchmark.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({
        "accession": ACCESSION,
        "name": "Ricin",
        "mode": MODE,
        "n_variants": N_VARIANTS,
        "elapsed_seconds": round(elapsed, 1),
        "result": result,
    }, indent=2, default=str), encoding="utf-8")
    print(f"\n  Saved: {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
