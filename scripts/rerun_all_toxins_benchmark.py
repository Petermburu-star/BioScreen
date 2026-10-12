"""
Re-run the full EvoDiff benchmark with corrected TOXIN_CONFIG.

Loads ESM-C 600M once, then for each toxin in TOXIN_CONFIG:
  - generates 3 preserved-mode EvoDiff variants
  - scores BLAST and the ESM fingerprint
  - writes incremental results so a crash doesn't lose everything

Writes: results/evodiff_16toxins_v2.json

The v2 file replaces the pre-correction benchmark. Original (now stale)
is kept as results/evodiff_16toxins_full.json for the record.

Run from repo root:

    python scripts/rerun_all_toxins_benchmark.py

Runtime: ~2 minutes per toxin, ~25 minutes total on the reference machine.
"""
import os

# Force UTF-8 before any print (Windows cp1252 crashes on box-drawing chars)
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
os.environ.setdefault("PYTHONUTF8", "1")

import argparse
import json
import pickle
import sys
import time
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "external" / "evodiff"))

import torch
from dotenv import load_dotenv

load_dotenv(ROOT / ".env")

from esm.models.esmc import EsmcForMaskedLM, EsmcTokenizer
from bioscreen.generate import generate_variants
from bioscreen.toxins import TOXIN_CONFIG


def embed_isolated(seq, model, tokenizer):
    inputs = tokenizer([seq], return_tensors="pt", padding=True)
    with torch.inference_mode():
        out = model(**inputs, output_hidden_states=True)
    return out.hidden_states[-1][0, 1:-1].cpu().numpy()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", default="preserved",
                        choices=["preserved", "aggressive"])
    args = parser.parse_args()
    MODE = args.mode

    print("=" * 78)
    print(f"BioScreen full-toxin EvoDiff benchmark (v2, corrected fingerprints)")
    print(f"mode: {MODE}")
    print("=" * 78)

    # --- Load ESM-C ---
    print()
    print("Loading ESM-C 600M...")
    model = EsmcForMaskedLM.from_pretrained("biohub/ESMC-600M")
    model.eval()
    tokenizer = EsmcTokenizer()
    print("  ready")

    def embed(seq):
        return embed_isolated(seq, model, tokenizer)

    # --- Load fingerprints ---
    fp_path = ROOT / "data" / "processed" / "toxin_fingerprints.pkl"
    with open(fp_path, "rb") as f:
        fp_data = pickle.load(f)
    fingerprints = fp_data["fingerprints"]
    print(f"  {len(fingerprints)} fingerprints loaded")

    # --- Resume support ---
    out_path = ROOT / "results" / f"evodiff_16toxins_v2_{MODE}.json"
    partial_path = ROOT / "results" / f"evodiff_16toxins_v2_{MODE}.partial.json"

    if partial_path.exists():
        print(f"  found partial results, resuming from {partial_path.name}")
        all_results = json.loads(partial_path.read_text(encoding="utf-8"))
        print(f"  {len(all_results)} toxins already complete")
    else:
        all_results = {}

    # --- Run each toxin ---
    items = list(TOXIN_CONFIG.items())
    total = len(items)

    for i, (acc, cfg) in enumerate(items, 1):
        if acc in all_results and "error" not in all_results[acc]:
            print(f"\n  [{i}/{total}] {acc} {cfg['name']} -- already done, skipping")
            continue

        print(f"\n  [{i}/{total}] {acc} {cfg['name']}")
        print(f"    region {cfg['region']} | {len(cfg['catalytic'])} catalytic residues")

        t0 = time.time()
        try:
            result = generate_variants(
                accession=acc,
                mode=MODE,
                n_variants=3,
                fingerprint_embeddings=fingerprints,
                embed_fn=embed,
            )
            elapsed = time.time() - t0

            all_results[acc] = {
                "accession": acc,
                "name": cfg["name"],
                "elapsed_seconds": round(elapsed, 1),
                "result": result,
            }
            # Incremental save
            partial_path.write_text(
                json.dumps(all_results, indent=2, default=str),
                encoding="utf-8",
            )

            caught = sum(1 for v in result["variants"] if v.get("fingerprint_caught"))
            evaded = sum(1 for v in result["variants"] if v.get("blast_evaded"))
            print(f"    done in {elapsed:.0f}s | BLAST evaded {evaded}/3 | FP caught {caught}/3")
        except Exception as e:
            print(f"    ERROR: {type(e).__name__}: {e}")
            traceback.print_exc()
            all_results[acc] = {
                "accession": acc,
                "name": cfg["name"],
                "error": f"{type(e).__name__}: {e}",
            }
            partial_path.write_text(
                json.dumps(all_results, indent=2, default=str),
                encoding="utf-8",
            )

    # --- Summary table ---
    print()
    print("=" * 78)
    print("SUMMARY")
    print("=" * 78)

    total_caught = 0
    total_evaded = 0
    total_variants = 0
    failures = []

    print(f"  {'toxin':<24} {'BLAST evaded':>14}  {'FP caught':>10}")
    print("  " + "-" * 55)
    for acc, entry in sorted(all_results.items()):
        if "error" in entry:
            print(f"  {entry['name']:<24}  ERROR")
            failures.append(entry["name"])
            continue
        variants = entry["result"]["variants"]
        caught = sum(1 for v in variants if v.get("fingerprint_caught"))
        evaded = sum(1 for v in variants if v.get("blast_evaded"))
        total_caught += caught
        total_evaded += evaded
        total_variants += len(variants)
        print(f"  {entry['name']:<24} {evaded:>10}/{len(variants)}  {caught:>10}/{len(variants)}")

    print("  " + "-" * 55)
    print(f"  {'TOTAL':<24} {total_evaded:>10}/{total_variants}  "
          f"{total_caught:>10}/{total_variants}")

    if failures:
        print(f"\n  FAILED TOXINS: {failures}")

    # --- Final save ---
    final = {
        "experiment": "evodiff_full_benchmark_v2",
        "date": "2026-10-12",
        "mode": MODE,
        "n_toxins": len(all_results),
        "n_variants": total_variants,
        "blast_evaded": total_evaded,
        "fingerprint_caught": total_caught,
        "note": ("Computed with corrected TOXIN_CONFIG "
                 "(session 5, commit 927154c). Supersedes "
                 "evodiff_16toxins_full.json."),
        "results": all_results,
    }
    out_path.write_text(json.dumps(final, indent=2, default=str), encoding="utf-8")
    print(f"\n  saved final: {out_path.relative_to(ROOT)}")

    # Clean up partial
    if partial_path.exists():
        partial_path.unlink()
        print("  removed partial file")

    print()
    print("=" * 78)
    print(f"  BLAST evaded:       {total_evaded} / {total_variants}")
    print(f"  Fingerprint caught: {total_caught} / {total_variants}")
    print("=" * 78)


if __name__ == "__main__":
    main()
