# BioScreen

**Functional screening and actuarial risk pricing for DNA synthesis biosecurity.**

A working detector for AI-designed toxin variants that sequence screening cannot
see — validated against Microsoft's own EvoDiff model — plus a self-financing
economic mechanism that funds global biosecurity capacity.

---

## The problem

Current DNA synthesis screening relies on sequence similarity to known hazards.
AI protein design can now rewrite a toxin's primary sequence while preserving
its function. A 2025 Microsoft-led red-teaming study demonstrated that **up to
100% of AI-designed ricin variants evaded major commercial screening tools**.

As of October 2026, no working countermeasure has been deployed. Microsoft's
patches are documented stopgaps. DeepMind's watermarking only works for
cooperative designers. The ADAPT programme is unfunded. Every deployed tool
remains sequence-based.

---

## The solution

BioScreen has two layers.

**Layer 1 — Functional screening.** Embed proteins with ESM-C 600M and run
**two independent detectors**:

1. **Contrastive** — the difference between best-match similarity to 16
   controlled toxins and best-match similarity to 25 safe reference proteins.
2. **Ordered catalytic fingerprint** — a novel method that compares only the
   residues within ±6 of each toxin's literature-curated catalytic residues,
   with a sequence-order consistency requirement to reject false positives.

A query is flagged if **either** detector fires.

**Layer 2 — Actuarial pricing.** Convert the screening signal and order context
into a risk-priced premium. A 0.12% loading on synthesis orders funds a global
biosecurity pool.

---

## Validated results

| Test | Result |
|---|---|
| EvoDiff ricin variants, BLAST evaded | 5/5 |
| EvoDiff ricin variants, fingerprint caught | **5/5** |
| Full 16-toxin EvoDiff benchmark, BLAST evaded | 36/36 (100%) |
| Full 16-toxin EvoDiff benchmark, fingerprint caught | **26/36 (72%)** |
| Circularity control (aggressive variants) | Correctly cleared at 0.50 vs 0.92 |
| 60 human proteins vs 4 fingerprints | **Zero false positives** |
| Generalization across 4 toxin mechanisms | Separation 0.62 – 0.71 |
| Pricing mechanism | 0.12% loading, $14.9M/year projected |

Every number is derived from `results/*.json`. EvoDiff is stochastic — re-runs produce new variants with similar qualitative behaviour, not identical scores. The deterministic parts (fingerprint scoring, specificity on 60 human proteins, pricing) reproduce exactly.

### Per-toxin EvoDiff detection

| Toxin | Caught | FP range |
|---|---|---|
| Ricin | 3/3 | 0.89 – 0.95 |
| Botulinum A | 3/3 | 0.88 – 0.93 |
| Botulinum B | 2/3 | 0.84 – 0.89 |
| Botulinum E | 2/3 | 0.48 – 0.90 |
| Shiga 1A | 2/3 | 0.85 – 0.94 |
| Shiga 2A | 2/3 | 0.84 – 0.91 |
| Diphtheria | 2/3 | 0.75 – 0.92 |
| Abrin | 1/3 | 0.80 – 0.90 |
| Modeccin | 2/3 | 0.78 – 0.91 |
| Anthrax LF | 3/3 | 0.87 – 0.92 |
| Anthrax EF | 3/3 | 0.88 – 0.96 |
| Tetanus | 1/3 | 0.49 – 0.90 |

---

## What's novel

- **Ordered catalytic fingerprint** — a data-driven representation of the
  structural signature of a toxin's active site, with a sequence-order
  consistency constraint that rejects safe proteins while keeping functional
  variants.
- **Detection of EvoDiff variants BLAST cannot see** — reproduces and answers
  Microsoft's finding.
- **Two independent detection layers OR-combined** — contrastive for
  interpretability, fingerprint for sharp discrimination. Defense in depth
  inside one tool.
- **Actuarial pricing layer** — no existing biosecurity tool prices risk.
- **Live 5-tab Studio** — every number computed live, every step traced.

---

## Install

```bash
git clone https://github.com/Petermburu-star/BioScreen
cd BioScreen
python -m venv .venv
.venv\Scripts\activate           # Windows
# source .venv/bin/activate      # macOS / Linux

pip install -r requirements.txt
echo "HF_TOKEN=hf_your_token_here" > .env
