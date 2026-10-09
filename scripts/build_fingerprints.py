"""
Build the fingerprint pickle — catalytic-window definition.

For each toxin, extracts the fingerprint as the residues within ±6 of
each catalytic residue (the definition used in PROTOCOL.md and validated
against EvoDiff variants).

Writes: data/processed/toxin_fingerprints.pkl
"""
import pickle
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
import pandas as pd
import torch

from esm.models.esmc import EsmcForMaskedLM, EsmcTokenizer
from bioscreen.toxins import TOXIN_CONFIG


WINDOW = 6   # ±6 residues around each catalytic site


def embed_isolated(seq, model, tokenizer):
    inputs = tokenizer([seq], return_tensors="pt", padding=True)
    with torch.inference_mode():
        out = model(**inputs, output_hidden_states=True)
    return out.hidden_states[-1][0, 1:-1].cpu().numpy()


def main():
    print("Loading ESM-C 600M...")
    model = EsmcForMaskedLM.from_pretrained("biohub/ESMC-600M")
    model.eval()
    tokenizer = EsmcTokenizer()

    tox_df = pd.read_csv(ROOT / "data" / "raw" / "toxin_references.csv")
    seqs = {r["accession"]: r["sequence"] for _, r in tox_df.iterrows()}

    fingerprints = {}
    names = {}
    regions = {}

    print("\nExtracting catalytic-window fingerprints:")
    for acc, cfg in TOXIN_CONFIG.items():
        if acc not in seqs:
            print(f"  [SKIP] {acc} not in reference CSV")
            continue

        lo, hi = cfg["region"]
        region_seq = seqs[acc][lo:hi]
        L = len(region_seq)

        # Catalytic residues relative to region
        cat_rel = [c - lo for c in cfg["catalytic"] if lo <= c < hi]
        if not cat_rel:
            print(f"  [SKIP] {acc} ({cfg['name']}): no catalytic residues in region")
            continue

        # Fingerprint = catalytic windows ±WINDOW
        fp_idx = sorted(set(
            p + off for p in cat_rel for off in range(-WINDOW, WINDOW + 1)
            if 0 <= p + off < L
        ))

        # Embed region in ISOLATED context (critical — matches the protocol)
        region_emb = embed_isolated(region_seq, model, tokenizer)
        fingerprint = region_emb[fp_idx]

        fingerprints[acc] = fingerprint
        names[acc] = cfg["name"]
        regions[acc] = (lo, hi)

        print(f"  [OK] {acc:<10} {cfg['name'][:34]:<34}  "
              f"region={L:>4}  catalytic={len(cat_rel)}  fp={len(fp_idx)}")

    out_path = ROOT / "data" / "processed" / "toxin_fingerprints.pkl"
    out_path.write_bytes(pickle.dumps({
        "fingerprints": fingerprints,
        "names": names,
        "regions": regions,
        "definition": "catalytic_window",
        "window": WINDOW,
    }))
    print(f"\n[OK] wrote {out_path}  ({out_path.stat().st_size:,} bytes)")
    print(f"     {len(fingerprints)} fingerprints")
    print(f"     definition: catalytic windows ±{WINDOW}")


if __name__ == "__main__":
    main()
