"""
Build the reference embeddings pickle.

Downloads the 16 controlled toxins and 25 safe reference proteins from
UniProt, embeds each with ESM-C 600M, and writes:

    data/processed/per_residue_embeddings.pkl

Run from repo root:

    python scripts/build_references.py

First run downloads ~2.4 GB of ESM-C weights from Hugging Face. Subsequent
runs load from cache in ~10 seconds.
"""
import os
import pickle
import sys
from pathlib import Path

# --- repo root ---
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
import requests
import torch
from dotenv import load_dotenv

load_dotenv(ROOT / ".env")

from esm.models.esmc import EsmcForMaskedLM, EsmcTokenizer


TOXIN_ACCESSIONS = {
    "P02879": "Ricin",
    "P0DPI1": "Botulinum neurotoxin type A",
    "P10844": "Botulinum neurotoxin type B",
    "P0DPI0": "Botulinum neurotoxin type E",
    "P10149": "Shiga-like toxin 1 subunit A (Stx1A)",
    "P69179": "Shiga-like toxin 1 subunit B (Stx1B)",
    "P09385": "Shiga-like toxin 2 subunit A (Stx2A)",
    "P09386": "Shiga-like toxin 2 subunit B (Stx2B)",
    "Q8X531": "Shiga toxin 2B variant",
    "P00588": "Diphtheria toxin",
    "P11140": "Abrin-a chain A",
    "P33183": "Modeccin",
    "P13423": "Anthrax protective antigen",
    "P15917": "Anthrax lethal factor",
    "P40136": "Anthrax edema factor",
    "P04958": "Tetanus toxin",
}

SAFE_ACCESSIONS = {
    # Structural / housekeeping
    "P69905": "Hemoglobin subunit alpha (human)",
    "P68871": "Hemoglobin subunit beta (human)",
    "P0DP23": "Calmodulin (human)",
    "P68133": "Actin, alpha skeletal muscle (human)",
    "P07437": "Tubulin beta chain (human)",
    "P60709": "Actin, cytoplasmic 1 (human)",
    "P06733": "Alpha-enolase (human)",
    "P00558": "Phosphoglycerate kinase 1 (human)",
    "P04406": "GAPDH (human)",
    "P68104": "Elongation factor 1-alpha 1 (human)",
    # Lab standards
    "P00698": "Lysozyme C (chicken)",
    "P00918": "Carbonic anhydrase 2 (human)",
    "P61626": "Lysozyme C (human)",
    "P02769": "Serum albumin (bovine)",
    "P42212": "GFP (Aequorea victoria)",
    "P00722": "Beta-galactosidase (E. coli)",
    "P01857": "Ig gamma-1 chain C region (human)",
    "P0DTC2": "Spike glycoprotein (SARS-CoV-2)",
    # Peptide hormones
    "P01317": "Insulin (human) mature",
    "P01308": "Preproinsulin (human)",
    "P01315": "Insulin A chain (human)",
    "P01343": "IGF-1 (human)",
    "P01344": "IGF-2 (human)",
    "P01275": "Glucagon (human)",
    "P0DTC2": "Spike glycoprotein (SARS-CoV-2)",
}


def fetch_fasta(accession: str) -> str:
    r = requests.get(
        f"https://rest.uniprot.org/uniprotkb/{accession}.fasta",
        timeout=30,
    )
    r.raise_for_status()
    lines = r.text.strip().split("\n")
    return "".join(lines[1:])


def embed_per_residue(sequence, model, tokenizer):
    inputs = tokenizer([sequence], return_tensors="pt", padding=True)
    with torch.inference_mode():
        out = model(**inputs, output_hidden_states=True)
    return out.hidden_states[-1][0, 1:-1].cpu().numpy()


def main():
    out_path = ROOT / "data" / "processed" / "per_residue_embeddings.pkl"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print("Loading ESM-C 600M...")
    model = EsmcForMaskedLM.from_pretrained("biohub/ESMC-600M")
    model.eval()
    tokenizer = EsmcTokenizer()
    print("  ready\n")

    def build_set(accessions, label):
        result = {"accessions": [], "names": [], "embeddings": {}}
        for acc, name in accessions.items():
            try:
                seq = fetch_fasta(acc)
            except Exception as e:
                print(f"  ✗ {label} {acc} ({name}): {e}")
                continue
            emb = embed_per_residue(seq, model, tokenizer)
            result["accessions"].append(acc)
            result["names"].append(name)
            result["embeddings"][acc] = emb
            print(f"  ✓ {label} {acc}  {name[:40]:<40}  shape={emb.shape}")
        return result

    print(f"Fetching {len(TOXIN_ACCESSIONS)} toxins...")
    toxins = build_set(TOXIN_ACCESSIONS, "toxin")
    print()

    print(f"Fetching {len(SAFE_ACCESSIONS)} safe references...")
    safes = build_set(SAFE_ACCESSIONS, "safe")
    print()

    payload = {"toxins": toxins, "safes": safes}
    out_path.write_bytes(pickle.dumps(payload))
    print(f"✓ wrote {out_path}  ({out_path.stat().st_size:,} bytes)")
    print(f"  toxins: {len(toxins['accessions'])}  safes: {len(safes['accessions'])}")


if __name__ == "__main__":
    main()
