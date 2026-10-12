# BioScreen — Re-run Protocol

Exact commands to rebuild the engine from scratch and re-run each validation.
Everything below assumes you are at the repo root:

    cd BioScreen

For first-time environment setup (venv, dependencies, HuggingFace token,
EvoDiff clone), see `docs/STUDIO.md` § First-time setup.

---

## 1. Build the reference embeddings

    python scripts/build_references.py

Fetches 16 controlled toxins + 25 safe reference proteins from UniProt,
embeds each with ESM-C 600M in isolated context, and writes:

    data/processed/per_residue_embeddings.pkl

One-time cost: ~2.4 GB model download from HuggingFace (cached after).
Subsequent runs: ~10 seconds.

---

## 2. Build the fingerprints

    python scripts/build_fingerprints.py

For each catalytic toxin, extracts the residues within ±6 of each
literature-curated catalytic residue (`src/bioscreen/toxins.py`), reads the
corresponding rows from the reference pickle, and writes:

    data/processed/toxin_fingerprints.pkl

Deterministic — same input, same output. Re-running produces an identical
pickle.

---

## 3. Re-run the ricin benchmark

    python scripts/rerun_ricin_benchmark.py

Generates 3 EvoDiff variants of the ricin A-chain region in preserved mode,
scores them with BLAST and the ordered catalytic fingerprint, and reports:

- sequence identity to native
- BLAST similarity (the incumbent's signal)
- fingerprint score

Writes:

    results/rerun_ricin_benchmark.json

**EvoDiff is stochastic.** Re-runs produce new variants each time. Expect
similar *qualitative* behaviour (BLAST evades, fingerprint catches) but not
identical scores. The current recorded run:

| Variant | Identity | BLAST | FP | Verdict |
|---|---|---|---|---|
| v0 | 1.50% | 0.1866 | 0.9094 | EVADED / CAUGHT |
| v1 | 1.50% | 0.1829 | 0.9104 | EVADED / CAUGHT |
| v2 | 1.12% | 0.1924 | 0.9312 | EVADED / CAUGHT |

---

## 4. Re-run the full 16-toxin benchmark

The 12 toxins with catalytic residues each get 3 EvoDiff variants (36 total).
Preserved mode: catalytic windows protected, everything else redesigned.

Expected outcomes (recorded in `results/evodiff_16toxins_full.json`):

- BLAST evaded: 36 / 36 (100%)
- Fingerprint caught: 26 / 36 (72%)

Per-toxin breakdown is in `docs/FINDINGS.md`.

**Runtime:** ~2 minutes per toxin × 12 = ~25 minutes on the reference machine.

---

## 5. Deterministic checks (no EvoDiff)

These reproduce exactly — no stochastic component:

### 5a. Specificity on 60 human proteins

    results/specificity_60human.json

60 diverse human proteins scored against 4 fingerprints (ricin, botulinum A,
shiga 1A, diphtheria). Zero false positives; smallest margin below the 0.85
threshold is +0.39.

### 5b. Circularity control

    results/fingerprint_discrimination_test.json

Preserved variants score ~0.925; aggressive variants (catalytic window
destroyed) score ~0.495 — same level as safe human proteins. Gap: +0.43.

### 5c. Generalization across mechanisms

    results/generalization_4toxins.json

Ricin, botulinum A, shiga 1A, diphtheria. Separation from best human match:
0.62 – 0.71.

### 5d. Pricing mechanism

    results/pricing_validation.json

5,000 simulated orders → 0.120% aggregate loading → $14.9M/year global pool
projection.

---

## 6. Environment versions

Recorded on the development machine:

    python  3.13
    torch   2.11.0
    transformers 4.57.6
    esm     3.4.1.post1
    numpy, pandas, biopython, flask, python-dotenv, tqdm (see requirements.txt)

EvoDiff is cloned from https://github.com/microsoft/evodiff (not pip-installed).
Exact commit hash not pinned — re-run behaviour may shift with upstream changes.

---

## 7. Known limitations

- **Auto-detection on custom sequences is not yet validated.** The current
  method picks the 42 residues least similar to the safe pool, which scatter
  across the whole region and produce an unusably large fingerprint on
  regions smaller than ~250 aa.
- **Botulinum and anthrax LF fingerprints are incomplete** versus literature
  (missing third zinc ligand / catalytic tyrosine). Fix pending.
- **Tetanus catalytic set** is literature-verified (H233, E234, H237, E271,
  Y375), but only catches 1 / 3 EvoDiff variants — likely a sensitivity issue,
  not an annotation issue.
- **No wet-lab validation.** All results are computational.
