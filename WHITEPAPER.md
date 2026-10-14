# BioScreen: Functional-Motif Screening of AI-Redesigned Toxin Sequences on Commodity Hardware

*Technical white paper — BioPREVAIL Biosecurity Innovation Challenge submission*

## 1. Abstract

AI protein-design tools can generate sequences that preserve a toxin's
functional motif while redesigning the surrounding scaffold, dropping
sequence identity to levels where similarity screening (BLAST and commercial
equivalents) loses detection. BioScreen screens for the conserved catalytic
*chemistry* instead: a gated, alignment-light functional fingerprint over
48 toxins in 4 mechanistic classes, scoring in ~0.02 s on a 16 GB laptop
with no network. We report a fully measured detection envelope, a
 UniProt-verified residue database, zero false positives across ~15,000
safe-protein windows, and a live, self-invalidating proof protocol for
reviewers.

## 2. Threat model

DNA synthesis orders are screened by sequence similarity. A designer with
EvoDiff-class tools can hold a toxin's catalytic window fixed and regenerate
the rest: our benchmark variants reach 1–25% identity to the native toxin.
Similarity screeners were not built for this. Additionally, many synthesis
facilities operate air-gapped; cloud screeners are unavailable there.
BioScreen targets both gaps: motif-level detection, offline deployment.

## 3. Positioning

The scoring core is deliberately classical — conserved-window identity,
BLOSUM62, and AAindex physicochemical correlation, combined in a **gated
product form**. The contribution is not a new algorithm class; it is (a) the
validation discipline — every claim pre-registered and measured, including
the failure modes; (b) false-positive engineering (the gate form dropped
safe-protein maxima from 0.61 to 0.41, enabling a 0.60 threshold with zero
false positives); and (c) a deployment model that runs during a synthesis
order review with no internet.

## 4. Database and verification

48 toxins across 4 functional classes (13 enzymatic, 12 superantigen,
16 pore-forming, 7 binding-only). For the enzymatic layer we:

1. Read `Active site` / `Binding site` features from UniProt (1-based
   precursor coordinates — the coordinate-convention failure mode is
   itself documented in our audit trail).
2. Gate on literature residue patterns (RIP `YYER`, zincin `HE..H`,
   diphtheria `HE`, ExoA `HYYE`…). A config that does not read back its
   literature residues cannot enter the database.
3. Verify arithmetically: under the gated form, a native self-score has a
   closed form `0.5 + 0.3·(mean self-BLOSUM + 3)/11 + 0.2`. All 12 enzymatic
   natives match prediction to 4 decimals (e.g., ExoA = 0.9659 exactly).

One toxin (Modeccin, P33183) failed verification and is **excluded** — the
tool does not vouch for what it cannot verify. This refusal is demonstrated
live: variants of the excluded toxin pass through unscored, by design.

## 5. The detection envelope (measured, pre-registered)

**Axis A — scaffold redesign + conservative window mutation.** With
catalytic positions protected, all 12 enzymatic fingerprints retain ≥90%
detection through 14 conservative substitutions in the ±6 window
(median k90 = 14, the top of our ladder — censored, not exhausted);
detection at 6 substitutions is 1.00 across the board. Small-window
fingerprints (≤16 residues) degrade by k = 10; large ones do not.

**Axis B — the single-codon adversary (exact, no sampling).** All 62
conservative single substitutions of catalytic residues, enumerated: 56/62
(90%) are flagged. Every multi-residue core (8/8 toxins) catches *every*
single substitution. The evasion set is precisely the small-core cluster:
diphtheria H>K (0.4754) and H>R (0.5514), shiga R>H (0.5435), anthrax EF
entirely (0.1188–0.2463).

**Full-core substitution:** replacing every catalytic residue conservatively
evades 12/12 fingerprints. The biological counterweight: these are the
substitutions that abolish catalysis (Y→F deletes the reactive hydroxyl;
E→D translocates the carboxylate by one methylene; HExxH→KExxR destroys
zinc coordination). The adversary who evades BioScreen has, with high
probability, not built a working toxin. We state this as a measured
boundary with a biology argument, not as a guarantee.

**Specificity.** Zero false positives across ~15,000 sliding windows of
safe proteins (25- and 60-protein panels, length-unbiased variant
included). Product-path margins: human TNF 0.4745, albumin 0.3769 vs
threshold 0.60.

**Real adversaries.** 39 EvoDiff preserved-mode variants across 13 toxins:
36/36 flagged through the production kernel (the 3 unflagged are variants
of the excluded toxin — coherent by construction). A live, reviewer-seeded
run reproduced the pipeline end-to-end in 73 s: preserved variant at 25.5%
identity flagged at 0.9455 with catalytic identity exactly 1.00; the
aggressive control cleared at 0.3777.

## 6. Live proof protocol

`scripts/live_proof.py` generates variants at a reviewer-chosen seed,
hashes them (SHA-256), confirms absence from all stored data, scores with
full component decomposition, and enforces its own controls — preserved
must flag, aggressive must clear — printing DEMO INVALIDATED otherwise.
We retain a deliberately mis-wired run in which the harness refused to
validate (`live_proof_miswired_refusal_example.json`). The scoring engine
also reproduces the published exact map live with |Δ| = 0.00000 and is
deterministic across re-scores.

## 7. Supplementary ESM layer

A per-residue embedding fingerprint (ESM-C) adds an independent signal.
On corrected configurations it flags 20/39 preserved variants (51.3%).
Its documented weakness motivates the CF layer: ESM embeddings are
context-sensitive — a variant preserving 50/50 fingerprint residues scored
0.6748 (below threshold) because surrounding redesign shifted every
embedding. Verdict logic is OR across layers; disagreement is reported.
Known issue, scheduled: the ESM builder applies a +1 offset to catalytic
positions; it is internally consistent (reference and query shifted
together) so published numbers stand as measured, but it will be corrected
with the next ESM re-validation.

## 8. Limitations

Sequence-level (fold-blind); cross-genus <40% identity out of scope;
full-core conservative substitution evades (Section 5 counterweight);
botulinum A/E attributable only at family level (identical catalytic
windows — reported honestly as a tie, not resolved arbitrarily);
Modeccin excluded; not wet-lab validated; not regulatory-grade.

## 9. Roadmap

1. **Structural layer** (acceleration-funded): ESMFold + pLDDT-gated
   geometric catalytic specs — already prototyped in this project — to
   close the core-substitution boundary and add fold evidence.
2. ESM wiring into the live-proof protocol.
3. New-class EvoDiff benchmark to full 48-toxin coverage.
4. Comparison against IBBIS commec.

## 10. Reproducibility

Every number in this paper is a JSON file: `uniprot_catalytic_audit.json`,
`conservative_substitution_envelope.json` (seed 20261013),
`adversarial_envelope_final.json`, `contextfree_specificity_*.json`,
`evodiff_16toxins_v2_*.json`, `product_freeze_report.json`,
`live_proof_*.json`. Protocol recipes: `PROTOCOL.md`.
