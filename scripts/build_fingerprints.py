"""
Build the fingerprint pickle.

Extracts an ordered catalytic fingerprint for each toxin from the
reference embeddings. Writes:

    data/processed/toxin_fingerprints.pkl

Run from repo root after build_references.py:

    python scripts/build_fingerprints.py
"""
import pickle
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np

from bioscreen.fingerprint import extract_fingerprint


N_RESIDUES = 42


def main():
    ref_path = ROOT / "data" / "processed" / "per_residue_embeddings.pkl"
    out_path = ROOT / "data" / "processed" / "toxin_fingerprints.pkl"

    assert ref_path.exists(), f"Missing {ref_path}. Run build_references.py first."

    data = pickle.loads(ref_path.read_bytes())
    toxin_embeddings = data["toxins"]["embeddings"]
    toxin_names = data["toxins"]["names"]
    toxin_accessions = data["toxins"]["accessions"]
    safe_embeddings = data["safes"]["embeddings"]

    # --- Safe pool ---
    safe_pool = np.vstack([safe_embeddings[a] for a in safe_embeddings])
    print(f"Safe pool: {safe_pool.shape}")

    # --- Extract per toxin ---
    fingerprints = {}
    names = {}

    for acc in toxin_accessions:
        emb = toxin_embeddings[acc]
        idx = extract_fingerprint(emb, safe_pool, n_residues=N_RESIDUES)
        fingerprints[acc] = emb[idx]
        names[acc] = toxin_names[toxin_accessions.index(acc)]
        print(f"  [OK] {acc:<10} {names[acc][:40]:<40}  fp shape={fingerprints[acc].shape}")

    payload = {
        "fingerprints": fingerprints,
        "names": names,
        "n_residues": N_RESIDUES,
    }
    out_path.write_bytes(pickle.dumps(payload))
    print(f"\n[OK] wrote {out_path}  ({out_path.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
