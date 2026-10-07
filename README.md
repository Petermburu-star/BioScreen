# BioScreen

**Functional screening and actuarial risk pricing for DNA synthesis biosecurity.**

A working prototype that catches AI-designed toxins by function, not sequence — and prices the risk through a self-financing actuarial mechanism.

---

## The problem

Current DNA synthesis screening relies on sequence similarity to known hazards. AI protein design can now rewrite a toxin's primary sequence while preserving its function. A 2025 Microsoft-led red-teaming study demonstrated that **up to 100% of AI-designed ricin variants evaded major commercial screening tools**.

The same weakness enables order fragmentation and benchtop synthesis.

---

## The solution

BioScreen has two layers.

**Layer 1 — Functional screening.** Embed proteins with ESM-C 600M, preserve per-residue embeddings, and score using an **ordered catalytic fingerprint** — a compact representation of the residues that determine the toxin's chemistry, in their correct sequence order.

**Layer 2 — Actuarial pricing.** Convert the screening signal and order context into a risk-priced premium. A 0.12% loading on synthesis orders funds a global biosecurity pool.

---

## Validated results

| Test | Result |
|---|---|
| EvoDiff variants evading BLAST | **5/5 caught** by ordered fingerprint |
| Aggressive redesigns (catalytic window destroyed) | **Correctly cleared** (mean 0.50) |
| Generalization across 4 toxin mechanisms | **Separation 0.62 – 0.71** |
| Human protein rejection | **Zero false positives across 60 proteins** |
| ProteinMPNN variants | 100% detection |

Every number reproducible from `results/*.json`. No hardcoded values.

---

## What's novel

- **Ordered catalytic fingerprint** — a data-driven representation of the residues that determine toxin chemistry, with sequence-order consistency to reject false positives.
- **Detection of EvoDiff variants BLAST cannot see** — reproduces and answers Microsoft's finding.
- **Actuarial pricing layer** — no existing biosecurity tool prices risk.

Full context, including a landscape analysis of existing solutions, is in `docs/WHITEPAPER.html`.

---

## Install

    git clone https://github.com/Petermburu-star/BioScreen
    cd BioScreen
    python -m venv .venv
    .venv\Scripts\activate
    pip install -r requirements.txt
    echo "HF_TOKEN=hf_your_token_here" > .env

You need a HuggingFace token. Get one at https://huggingface.co/settings/tokens and accept the license at https://huggingface.co/biohub/ESMC-600M.

---

## Run the Studio

    python scripts/start_studio.py

Then open http://127.0.0.1:5000. Full guide in `docs/STUDIO.md`.

---

## Use the Python API

    from bioscreen import BioScreenV2

    screener = BioScreenV2(
        model=model,
        tokenizer=tokenizer,
        per_residue_path="data/processed/per_residue_embeddings.pkl",
    )
    result = screener.screen("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ")
    # {'best_toxin_score': 0.87, 'hazard_delta': 0.05, 'flagged': True, ...}

---

## Structure

    BioScreen/
    ├── src/bioscreen/
    │   ├── screening_v2.py          contrastive per-residue screening
    │   ├── pricing.py               actuarial pricing framework
    │   └── studio/                  web UI (Flask)
    ├── scripts/
    │   └── start_studio.py          launch the Studio
    ├── results/                     all validation outputs (JSON)
    ├── docs/
    │   ├── WHITEPAPER.html          full technical white paper
    │   ├── STUDIO.md                Studio user guide
    │   └── FINDINGS.md              EvoDiff fingerprint validation
    ├── reports/                     pitch deck, plots
    └── data/
        ├── raw/                     UniProt CSVs (gitignored)
        └── processed/               embeddings, fingerprints (gitignored)

---

## Stack

- ESM-C 600M (EvolutionaryScale, MIT license, Hugging Face)
- EvoDiff (Microsoft, MIT license)
- UniProt reference database
- Python: transformers, torch, numpy, pandas, biopython, flask

---

## The EvoDiff test

See `docs/FINDINGS.md` for the full reproduction of Microsoft's finding and the ordered fingerprint result.

---

## License

MIT — see LICENSE.

---

## Author

Peter Mburu Kariuki — pmburu346@gmail.com

Pre-submission prototype for BioPREVAIL × CEPI Biosecurity Innovation Challenge 2026–27.
