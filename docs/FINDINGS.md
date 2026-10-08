# BioScreen — Findings

Technical results from the EvoDiff validation. Full reproduction steps below.

---

## Reproducing Microsoft's finding

Microsoft's 2025 red-teaming study showed that AI-designed ricin variants evade commercial sequence-similarity screening. We reproduced this with the public EvoDiff model (Microsoft, MIT license).

**Setup:**

- 12 toxins with literature-curated catalytic residues
- Ricin A-chain (267 aa), botulinum light chains (437 aa), shiga A subunits (315 aa), and others
- EvoDiff `inpaint_simple` redesigns scaffold regions
- Catalytic residues preserved (±8 window)
- 3 independent variants per toxin (36 variants total)

**BLAST detection:** 36/36 evaded (all scored below 0.30 local alignment threshold).

---

## The detection method

A **fingerprint** built from each toxin's catalytic window (residues within ±6 of each catalytic residue) and scored with two constraints:

1. **Per-residue similarity** — how well do fingerprint residues match the query?
2. **Order consistency** — do consecutive fingerprint residues land in the correct relative positions?

    score = top20_similarity × order_consistency

The order term eliminates false positives. GFP scored 0.92 on per-residue similarity but dropped to 0.45 once order was required.

---

## Validation results — 16 toxins

| Toxin | Variants caught | FP range |
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

**Aggregate: 26/36 caught (72%). BLAST evaded 36/36 (100%).**

---

## The circularity test

The EvoDiff variants preserved the catalytic residues used to build the fingerprint. Is the match trivial?

**The test:** generate "aggressive" variants where EvoDiff redesigns the whole A-chain including the catalytic window. If the fingerprint is only checking that we preserved what we preserved, these should also score high. If the fingerprint detects catalytic residues specifically, these should score low.

**The result:**

| Class | Mean fingerprint score |
|---|---|
| Preserved variants (catalytic window intact) | 0.925 |
| Aggressive variants (catalytic window destroyed) | 0.495 |
| **Gap** | **+0.429** |

Aggressive variants score at the same level as safe human proteins. The fingerprint is discriminating on catalytic residues.

---

## Generalization across mechanisms

Tested on four mechanistically distinct toxins:

| Toxin | Mechanism | Self | Best human | Separation |
|---|---|---|---|---|
| Ricin | N-glycosidase | 1.00 | 0.29 | +0.71 |
| Botulinum A | Zinc protease | 1.00 | 0.34 | +0.66 |
| Shiga 1A | rRNA depurination | 1.00 | 0.38 | +0.62 |
| Diphtheria | ADP-ribosylation | 1.00 | 0.32 | +0.68 |

All four pass. The method is not ricin-specific.

---

## Specificity — 60 human proteins

| Toxin | Self | Human mean | Human max | Margin below threshold |
|---|---|---|---|---|
| Ricin | 1.00 | 0.25 | 0.32 | +0.53 |
| Botulinum A | 1.00 | 0.27 | 0.39 | +0.46 |
| Shiga 1A | 1.00 | 0.25 | 0.46 | +0.39 |
| Diphtheria | 1.00 | 0.28 | 0.38 | +0.47 |

Zero false positives across 60 human proteins. Every fingerprint leaves at least 0.39 headroom below the 0.85 threshold.

---

## What's uncertain

1. **10 of 36 variants were missed (28%).** Most landed at 0.78–0.84 — one point below the 0.85 threshold. Threshold tuning may recover some without sacrificing specificity.

2. **Abrin is the weakest fingerprint** (1/3). Abrin and ricin are both plant RIPs; the abrin fingerprint may be partially shadowed by ricin's.

3. **Tetanus needs UniProt annotation verification.** Currently 1/3 with a best-guess catalytic set. UniProt active-site annotations will fix it.

4. **Adversarial robustness untested.** An attacker who knows about the fingerprint could attempt chemistry-preserving substitutions (E→D, R→K, Y→F). Whether the fingerprint catches those is unverified.

5. **Function of EvoDiff variants not experimentally validated.** The design preserves the catalytic window, but we haven't confirmed the variants fold and catalyze correctly.

6. **Only one design tool tested.** EvoDiff is the specific tool Microsoft used. Other tools (ProteinMPNN, RFdiffusion) may produce variants with different characteristics.

---

## Interpretation

The ordered catalytic fingerprint detects AI-designed variants of known toxins that BLAST cannot see. It generalizes across four mechanistically distinct toxins. It cleanly rejects safe proteins when the order-consistency term is included. It catches 72% of EvoDiff variants on average across 12 toxins — enough to prove the method, not enough to be production-grade.

This is not a general solution to AI-designed biothreats. It is a specific solution to the threat Microsoft identified: function-preserving, sequence-diverging variants of toxins already in a reference set.

---

## Reproduction

All results reproducible from `results/evodiff_16toxins_full.json`.

Full pipeline:

1. Fetch toxin sequences from UniProt
2. Install EvoDiff from GitHub (microsoft/evodiff)
3. For each toxin: extract region, build fingerprint, run `inpaint_simple` with catalytic window protected
4. Score with BLAST and fingerprint

See `docs/PROTOCOL.md` for the exact recipe.
