# BioScreen Studio

A local web interface for the BioScreen functional screening engine and
actuarial pricing framework. Designed for non-coders.

---

## What it does

BioScreen Studio lets you:

1. **Screen a protein sequence** for functional similarity to controlled toxins.
   Paste a sequence or load one from UniProt with one click.
2. **See the algorithm at work.** Two ribbons show per-residue match strength
   against the winning toxin and the winning safe protein. A third ribbon shows
   the *contrastive divergence* — the exact signal that drives the verdict.
3. **Price the biosecurity risk** in real time. Adjust verification status,
   history, organism, and order size; watch the tier and premium update live.

Every result is computed live from ESM-C 600M embeddings and UniProt sequences.
Nothing is hardcoded.

---

## First-time setup

### 1. Install dependencies

    cd C:\Users\Administrator\BioScreen
    python -m venv .venv
    .\.venv\Scripts\Activate.ps1
    pip install -r requirements.txt
    pip install flask python-dotenv

### 2. Get a HuggingFace token

The ESM-C 600M model is gated. You need a free HuggingFace account.

1. Sign up at https://huggingface.co
2. Go to https://huggingface.co/biohub/ESMC-600M and click "Agree and access repository"
3. Create a Read token at https://huggingface.co/settings/tokens
4. Save it to `.env` in the repo root:

    echo "HF_TOKEN=hf_your_token_here" > C:\Users\Administrator\BioScreen\.env

The `.env` file is gitignored. It never leaves your machine.

### 3. Generate the reference embeddings (first time only)

If `data/processed/per_residue_embeddings.pkl` does not exist, run the
notebook `notebooks/01_build_references.ipynb`. It downloads the toxin and
safe reference sets from UniProt and embeds them with ESM-C. Takes ~2 minutes.

---

## Running the Studio

### From the command line

    cd C:\Users\Administrator\BioScreen
    .\.venv\Scripts\Activate.ps1
    python start_studio.py

Then open http://127.0.0.1:5000 in any browser.

The first run downloads 2.4 GB of ESM-C weights (one-time). Subsequent runs
load from cache in ~10 seconds.

### From a Jupyter notebook

    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path.cwd() / "src"))

    from bioscreen.studio.launch import launch
    launch()

The Studio appears as an embedded iframe in the notebook cell.

---

## Using the interface

### Screening a sequence

1. **Paste** a protein sequence into the textarea — single-letter amino acid
   code, no spaces, no FASTA header.
2. **Or click a chip** — Ricin, Botulinum A, etc. The sequence is fetched live
   from UniProt.
3. **Click Screen.**

### Reading the result

- **Two ribbons** at the top and bottom show how strongly each position in your
  query aligns to the winning toxin (amber) and the winning safe protein (cyan).
  Brighter = stronger match.
- **The divergence ribbon** in the middle is the contrastive signal made
  visible. Above the centre line means the toxin wins at that position; below
  means the safe reference wins.
- **The verdict panel** shows whether the sequence is **Flagged** or **Clear**,
  and lights up the two gates that produced the decision:
  - **delta >= 0.05** — the hazard delta (toxin score minus safe score)
  - **toxin >= 0.90** — the absolute toxin similarity floor

  Both must open for a sequence to be flagged.

### Pricing an order

Adjust the controls on the right-hand side:

- **Verified customer** — toggle. Verified pay ~30% less.
- **History** — 0 to 1. Higher = better track record, lower premium.
- **Organism** — human (highest risk) through environmental (lowest).
- **Order size** — larger orders amplify risk (log-scaled).
- **Order cost** — the base order value in USD.

The **readout panel** updates live. Its color reflects the tier: green for LOW,
cyan for MEDIUM, amber for HIGH, red for CRITICAL or REJECTED.

---

## Interpreting a result

| Verdict | Meaning | Suggested action |
|---|---|---|
| **Clear** | No functional toxin signature detected | Standard processing |
| **Flagged** | The query carries a functional signature of a controlled toxin | Escalate to human review |

The Studio is a screening aid. It does not make accept/reject decisions on its
own. Every flagged order should go to a human reviewer before any action is taken.

---

## What is computed live

- **Model:** ESM-C 600M, loaded fresh from HuggingFace (cached after first run)
- **Reference sequences:** fetched from UniProt at build time
- **Screening scores:** per-residue contrastive alignment, computed for every query
- **Pricing:** re-computed on every control change
- **Example sequences:** fetched from UniProt REST API at click time

There is no embedded test data, no fixture results, no placeholder numbers.

---

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| Status pip stays amber | Model still loading | Wait 30s. If still amber, check the notebook cell for errors. |
| Status pip turns red | Engine error | Read the error text in the header. Usually a missing .env or stale pickle. |
| ModuleNotFoundError: flask | Flask not installed | pip install flask python-dotenv |
| 401 Unauthorized on model load | HF token missing or license not accepted | Re-check step 2 above |
| FileNotFoundError: per_residue_embeddings.pkl | Reference embeddings not built | Run notebooks/01_build_references.ipynb |
| Iframe blank in notebook | Jupyter blocking localhost iframe | Open http://127.0.0.1:5000 in a browser tab instead |

---

## Architecture

    Browser
       |  HTTP (localhost only)
       v
    Flask app (src/bioscreen/studio/app.py)
       |  /                -> studio.html
       |  /api/status      -> engine readiness + toxin/safe counts
       |  /api/warmup      -> triggers engine load in background
       |  /api/screen      -> runs the screening engine on a sequence
       |  /api/price       -> runs the actuarial pricing engine
       |  /api/example/... -> fetches live sequences from UniProt
       v
    bioscreen.screening_v2.BioScreenV2   -> contrastive per-residue alignment
    bioscreen.pricing.ActuarialPricer    -> catastrophe-insurance pricing

The Flask server binds to 127.0.0.1 only. It is not exposed to the network.
No data leaves the machine except the initial model download from HuggingFace
and the UniProt example fetches.

---

## License

MIT. See the repository root for details.
