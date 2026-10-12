# BioScreen Studio

A local web interface for the BioScreen functional screening engine and
actuarial pricing framework.

---

## What it does

BioScreen Studio has five tabs:

| Tab | Purpose |
|---|---|
| **Screen** | Paste a protein sequence (or load a reference). Get the contrastive verdict, per-residue ribbons, fingerprint score, tier, and premium. |
| **Pipeline trace** | Step-by-step timeline of the last Screen run — every intermediate value, every timing. |
| **Generate** | Run EvoDiff on a curated toxin, or on any custom sequence (paste or fetch from UniProt by accession). Preserved and aggressive modes. |
| **Adversarial** | Pre-computed 13-toxin EvoDiff comparison. BLAST evaded 39/39. ESM fingerprint 20/39 preserved, context-free 39/39 preserved / 0/39 aggressive. |
| **Re-run** | Clone/build commands and a live source viewer for the engine modules. |

Every number is computed live from ESM-C 600M embeddings and UniProt sequences.

---

## First-time setup

### 1. Install dependencies

    git clone https://github.com/Petermburu-star/BioScreen
    cd BioScreen
    python -m venv .venv
    .venv\Scripts\activate          # Windows
    source .venv/bin/activate        # macOS / Linux
    pip install -r requirements.txt

EvoDiff is not pip-installable. Clone it separately:

    git clone https://github.com/microsoft/evodiff external/evodiff

### 2. Get a HuggingFace token

The ESM-C 600M model is gated. You need a free HuggingFace account.

1. Sign up at https://huggingface.co
2. Go to https://huggingface.co/biohub/ESMC-600M and click "Agree and access repository"
3. Create a Read token at https://huggingface.co/settings/tokens
4. Save it to `.env` in the repo root:

    echo "HF_TOKEN=hf_your_token_here" > .env

The `.env` file is gitignored. It never leaves your machine.

### 3. Generate the reference data (first time only)

Two build steps, in order:

    python scripts/build_references.py

Downloads the 16 controlled toxins and 25 safe reference proteins from
UniProt, embeds each with ESM-C 600M, and writes
`data/processed/per_residue_embeddings.pkl`. Takes a few minutes; the
ESM-C weights download once (~2.4 GB) and cache after that.

    python scripts/build_fingerprints.py

Derives the ordered catalytic fingerprints from those embeddings and
`src/bioscreen/toxins.py`. Writes `data/processed/toxin_fingerprints.pkl`.

---

## Running the Studio

    python scripts/start_studio.py

Then open http://127.0.0.1:5000 in any browser.

The first run downloads ESM-C weights (one-time). Subsequent runs load from
cache in seconds.

---

## Using the interface

### Screening a sequence

1. **Paste** a protein sequence into the textarea — single-letter amino acid
   code, no spaces, no FASTA header.
2. **Or pick a reference** from the dropdown — Ricin, Botulinum A, etc. The
   sequence is fetched live from UniProt.
3. **Click Screen.**

### Reading the result

- **Two ribbons** show how strongly each position in your query aligns to the
  winning toxin (amber) and the winning safe protein (cyan).
- **The divergence ribbon** makes the contrastive signal visible. Above the
  centre line means the toxin wins at that position; below means the safe
  reference wins.
- **The verdict panel** shows whether the sequence is **Flagged** or **Clear**,
  and lights up the gates that produced the decision:

  | Layer | Gate |
  |---|---|
  | Contrastive | `hazard_delta >= 0.02` **and** `best_toxin_score >= 0.90` |
  | ESM fingerprint | `fingerprint_score >= 0.85` |
  | Context-free | `motif_passes` **and** `per_site_max >= 0.60` |

  A query is flagged if **either** layer fires.

### Pricing an order

Adjust the controls on the right-hand side:

- **Verified customer** — toggle. Verified customers pay 30% less.
- **History** — 0 to 1. Higher = better track record, up to 30% further discount.
- **Organism** — human (1.0×) through animal (0.6×), plant (0.3×), environmental (0.2×).
- **Order size** — larger orders amplify risk (log-scaled).
- **Order cost** — the base order value in USD.

The **readout panel** updates live. Its color reflects the tier: green for LOW,
cyan for MEDIUM, amber for HIGH, red for CRITICAL or REJECTED.

---

## Interpreting a result

| Verdict | Meaning | Suggested action |
|---|---|---|
| **Clear** | Neither detection layer fired | Standard processing |
| **Flagged** | Contrastive or fingerprint layer fired | Escalate to human review |

The Studio is a screening aid. It does not make accept/reject decisions on its
own. Every flagged order should go to a human reviewer before any action is taken.

---

## What is computed live

- **Model:** ESM-C 600M, loaded fresh from HuggingFace (cached after first run)
- **Screening scores:** contrastive alignment and fingerprint scoring, per query
- **Pricing:** re-computed on every control change
- **Example sequences:** fetched from UniProt REST API at click time

The `results/*.json` files are recorded outputs from past EvoDiff runs — they
document what the engine produced, not fixtures the UI reads from.

---

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| Status pip stays amber | Model still loading | Wait 30s. If still amber, check the terminal for errors. |
| Status pip turns red | Engine error | Read the error text in the header. Usually a missing `.env` or stale pickle. |
| `ModuleNotFoundError: flask` | Flask not installed | `pip install -r requirements.txt` |
| 401 Unauthorized on model load | HF token missing or license not accepted | Re-check step 2 above |
| `FileNotFoundError: per_residue_embeddings.pkl` | Reference embeddings not built | `python scripts/build_references.py` |
| `FileNotFoundError: toxin_fingerprints.pkl` | Fingerprints not built | `python scripts/build_fingerprints.py` |
| `ModuleNotFoundError: evodiff` | EvoDiff not cloned | `git clone https://github.com/microsoft/evodiff external/evodiff` |

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
       |  /api/source      -> serves engine source for the Re-run tab
       |  /api/example/... -> fetches live sequences from UniProt
       v
    bioscreen.screening_v2.BioScreenV2   -> contrastive per-residue alignment
    bioscreen.fingerprint                -> ordered catalytic fingerprint (ESM-based)
    bioscreen.contextfree                -> context-free fingerprint (motif + per-site)
    bioscreen.pricing.ActuarialPricer    -> catastrophe-insurance pricing

The Flask server binds to 127.0.0.1 only. It is not exposed to the network.
No data leaves the machine except the initial model download from HuggingFace
and the UniProt example fetches.

---

## License

MIT. See the repository root for details.
