# BioScreen — Findings

Honest record of what we found, what we can claim, and what we cannot.

Last updated: 2026-10-12 (session 5)

---

## The headline result

**Full EvoDiff benchmark — 13 toxins, 39 variants per mode.**

| Mode | BLAST evaded | ESM fingerprint | Context-free |
|---|---|---|---|
| **Preserved** | 39/39 (100%) | **20/39 (51.3%)** | **39/39 (100%)** |
| **Aggressive** | 39/39 (100%) | 0/39 | **0/39** |

**Two-layer separation: 39/39 preserved vs 0/39 aggressive (100 pp).**

The ESM fingerprint alone catches about half of preserved-mode variants.
The context-free layer catches all of them. Combined, the detector works.

---

## What "preserved" and "aggressive" mean

**Preserved mode:** EvoDiff is told to protect ±8 residues around each
catalytic position, and to redesign everything else. The catalytic window
is preserved by construction.

**Aggressive mode:** EvoDiff redesigns straight through the catalytic
window. The site is destroyed.

These are designed extremes. Real AI-designed threats could fall anywhere
between them.

---

## Reproducing the 51.3% (correction to submitted claim)

The BioPREVAIL × CEPI submission reported **26/36 (72%)** for the ESM
fingerprint. That number was inflated by a config bug:

- 7 of 13 `TOXIN_CONFIG` catalytic positions were using mature-chain
  numbering where the region slicing assumed precursor coordinates.
- The wrong positions happened to sit in regions EvoDiff preserved well,
  so the fingerprint caught variants by proxy — not by detecting the
  actual active site.
- After correcting positions against UniProt active-site annotations,
  the same benchmark drops to **20/39 (51.3%)**.

The corrected number is the honest one.

---

## Why the ESM fingerprint drops

ESM-C embeddings are **context-dependent**. A variant that preserves
every fingerprint residue can still score below threshold if the
surrounding scaffold is redesigned.

**Direct evidence:** P11439 EvoDiff variant v1 preserved all 50
fingerprint residues verbatim. It scored 0.6748 — below the 0.85
threshold — because redesigning the 158 non-fingerprint residues shifted
ESM-C's embedding at every fingerprint position.

Preserving residues does not guarantee preserving the embedding
neighborhood. That is a fundamental property of contextual protein
language models, not a bug.

---

## The context-free layer — what it does

A residue-identity + physicochemical check that does not depend on any
protein language model.

**Signals:**
- Residue identity at catalytic positions
- BLOSUM62 substitution scores
- AAindex physicochemical profile correlation over the fingerprint window
- Mechanism-specific catalytic motif patterns (HExxH, YYER, ER)

**Verdict form:** identity acts as a multiplicative gate on the inner
score. If no catalytic residue is preserved, the score is near zero
regardless of other signals.

---

## The context-free layer — honest assessment

**The components are not new.** BLOSUM62 (1992), AAindex scales
(1968–1982), catalytic motifs (standard enzymology). What is new is
the specific combination and its application to biosecurity screening.

**Validated on the designed test:**
- Preserved: 39/39 fires
- Aggressive: 0/39 fires
- Separation: 100 percentage points

**Validated for specificity:**
- 25 safe human proteins, ~2,500 sliding windows
- 0 false positives

**Not validated:**
- Adversarial variants (preserve residues, break geometry)
- Cross-mechanism transfer (does a ricin fingerprint false-positive on
  a botulinum variant?)
- Real AI-designed threats beyond EvoDiff
- Misfolded proteins with preserved sequence (will false-positive)

---

## Frustration analysis — preliminary

Ran FrustrAI-Seq (Leusch et al., 2026) on 4 toxins: ricin, botulinum A,
shiga 1A, Pseudomonas ExoA.

**Finding:** 4 of 13 catalytic residues fall in the top-5% of frustration
(6.15× enrichment over random, p ~ 0.003).

**Stronger pattern:** all 4 are **catalytic glutamates** — the
proton-transfer residue in each mechanism. Tyr, His, and Arg catalytic
residues do not show the signal.

**Mechanistic interpretation:** the catalytic glutamate must be
deprotonated (charged) at physiological pH, often in a buried active
site. That conflicts with folding stability. The frustration is the
energetic cost of that conflict — consistent with Galpern et al. (2026)
"dark energy" theory.

**Caveats:**
- n=13 catalytic residues across 4 toxins. Small sample.
- Post-hoc observation, not pre-registered.
- Untested on AI-designed variants.
- Confirmation on all 13 toxins required.

**Status:** documented as a promising direction, not yet integrated into
the detection pipeline.

---

## Reference set

**17 controlled toxins, 13 with catalytic fingerprints:**

| Accession | Name | Mechanism | FP size |
|---|---|---|---|
| P00588 | Diphtheria | ADP-ribosylation | 26 |
| P02879 | Ricin | N-glycosidase | 42 |
| P04958 | Tetanus | Zinc protease | 43 |
| P09385 | Shiga 2A | rRNA depurination | 16 |
| P0DPI0 | Botulinum E | Zinc protease | 17 |
| P0DPI1 | Botulinum A | Zinc protease | 17 |
| P10149 | Shiga 1A | rRNA depurination | 16 |
| P10844 | Botulinum B | Zinc protease | 30 |
| P11140 | Abrin | N-glycosidase | 42 |
| P11439 | Pseudomonas ExoA | ADP-ribosylation | 50 |
| P15917 | Anthrax LF | Zinc protease | 30 |
| P33183 | Modeccin | N-glycosidase | 42 (unresolved) |
| P40136 | Anthrax EF | Adenylate cyclase | 13 |

**Non-catalytic (no fingerprint):** Shiga 1B (P69179), Shiga 2B
(P09386), Shiga 2B variant (Q8X531), Anthrax PA (P13423). These are
binding or pore-forming subunits with no catalytic residues.

---

## Config corrections (session 5)

7 of 13 catalytic configurations were corrected against UniProt
active-site annotations:

| Toxin | Was | Now |
|---|---|---|
| Diphtheria | [21, 148] | [53, 180] |
| Ricin | [80, 123, 177, 180] | [115, 158, 212, 215] |
| Shiga 1A | [167, 170] | [189, 192] |
| Shiga 2A | [167, 170] | [189, 192] |
| Anthrax LF | [686, 690, 728, 735] | [719, 720, 723, 768] |
| Botulinum B | [223, 224, 227] | [230, 231, 234, 268] |
| Modeccin | unresolved | unresolved (flagged) |

Verified correct: tetanus, botulinum A/E, abrin, P11439, anthrax EF.

---

## What we cannot claim

- **Not "functional screening."** The context-free layer scores
  sequence-level features. It cannot verify 3D active-site geometry.
  Misfolded proteins with preserved sequences will false-positive.
- **Not a general-purpose detector.** Validated on 13 toxins from
  5 mechanistic classes. Novel mechanisms with no reference are not
  covered.
- **Not regulatory-grade.** Proof-of-concept.
- **Not wet-lab validated.** All results computational.

---

## Open questions

1. **Adversarial variants.** Preserve residues, break geometry. Does the
   context-free layer still fire?
2. **Cross-mechanism transfer.** Does a ricin fingerprint false-positive
   on botulinum variants?
3. **60-protein specificity panel** for the context-free layer (25 done).
4. **Frustration enrichment on all 13 toxins** (4 done).
5. **Structural layer** (ESMFold + active-site RMSD) — the only way to
   honestly claim functional detection.

---

## Files

- `results/evodiff_16toxins_v2_preserved.json` — preserved-mode benchmark
- `results/evodiff_16toxins_v2_aggressive.json` — aggressive-mode benchmark
- `results/contextfree_preserved_vs_aggressive.json` — separation test
- `results/contextfree_specificity_v2.json` — 0/13 false positives
- `results/frustration_enrichment.json` — preliminary frustration result
- `results/p11439_evodiff_benchmark.json` — the case that motivated the CF layer
