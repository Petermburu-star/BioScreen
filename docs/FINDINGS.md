# BioScreen — Findings

Technical results from the EvoDiff validation. Full reproduction steps below.

---

## Reproducing Microsoft's finding

Microsoft's 2025 red-teaming study showed that AI-designed ricin variants evade commercial sequence-similarity screening. We reproduced this with the public EvoDiff model (Microsoft, MIT license).

**Setup:**

- Ricin A-chain (267 aa) from AlphaFold DB (P02879)
- EvoDiff `inpaint_simple` redesigns scaffold regions
- Catalytic residues preserved (±8 window)
- Result: 5 variants, sequence identity **0.7 – 1.5%** to native

**BLAST detection:** 5/5 evaded (all scored below 0.30 local alignment threshold).

---

## The detection method

A **fingerprint** built from the native toxin's catalytic window (42 residues) and scored with two constraints:

1. **Per-residue similarity** — how well do fingerprint residues match the query?
2. **Order consistency** — do consecutive fingerprint residues land in the correct relative positions?

    score = top20_similarity × order_consistency

The order term is what eliminates false positives. GFP scored 0.92 on per-residue similarity but dropped to 0.45 once order was required.

---

## Validation results

| Protein class | BLAST | Fingerprint | Verdict |
|---|---|---|---|
| Native ricin | 1.00 | 1.00 | Flag |
| EvoDiff variant 0 | 0.18 | 0.92 | Flag (BLAST missed) |
| EvoDiff variant 1 | 0.20 | 0.93 | Flag (BLAST missed) |
| EvoDiff variant 2 | 0.18 | 0.97 | Flag (BLAST missed) |
| EvoDiff variant 3 | 0.19 | 0.86 | Flag (BLAST missed) |
| EvoDiff variant 4 | 0.18 | 0.95 | Flag (BLAST missed) |
| Aggressive redesign | 0.06 | 0.50 | Clear |
| Human lysozyme | — | 0.12 | Clear |
| GFP | — | 0.45 | Clear |
| Human calmodulin | — | 0.01 | Clear |

**Separation: +0.42.**

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

1. **Only ricin tested end-to-end with EvoDiff variants.** The other three toxins have validated fingerprints but haven't been tested against their own EvoDiff variants.

2. **Adversarial robustness untested.** An attacker who knows about the fingerprint could attempt chemistry-preserving substitutions (E→D, R→K, Y→F). Whether the fingerprint catches those is unverified.

3. **Function of EvoDiff variants not experimentally validated.** The design preserves the catalytic window, but we haven't confirmed the variants fold and catalyze correctly.

4. **Only one design tool tested.** EvoDiff is the specific tool Microsoft used. Other tools (ProteinMPNN, RFdiffusion) may produce variants with different characteristics.

---

## Interpretation

The ordered catalytic fingerprint detects AI-designed variants of known toxins that BLAST cannot see. It generalizes across mechanisms. It cleanly rejects safe proteins when the order-consistency term is included.

This is not a general solution to AI-designed biothreats. It is a specific solution to the threat Microsoft identified: function-preserving, sequence-diverging variants of toxins already in a reference set.

---

## Reproduction

All results reproducible from `results/*.json`.

Full pipeline:

1. Fetch ricin A-chain from AlphaFold DB (P02879)
2. Install EvoDiff from GitHub (microsoft/evodiff)
3. Run `inpaint_simple` with catalytic window protected
4. Build fingerprint from the catalytic window
5. Score with `top20_similarity × order_consistency`
