# BioScreen — Offline Biosecurity Screening

**Paste a protein sequence → FLAGGED / CLEAR in ~0.02 s. Entirely on a 16 GB laptop. No internet, no GPU, no cloud.**

BioScreen screens protein sequences for toxin functional motifs that survive
AI-driven sequence redesign. Where sequence-similarity screening (BLAST-class)
loses the signal — AI-designed variants at 1–25% identity — BioScreen detects
the conserved catalytic chemistry and names the toxin family.

## 30-second start

    double-click START_BIOSCREEN.bat      → browser opens at http://127.0.0.1:5055
    click the 10 demo buttons top to bottom

No internet required. The verdict card shows the score, the matched toxin
family, and the verification provenance of the fingerprint that fired.

## Live reviewer proof

    python scripts/live_proof.py --accession P02879 --seed <YOUR_NUMBER> --n 1

You choose the toxin and the seed. Variants are generated live (EvoDiff),
hashed, confirmed absent from every stored dataset, and scored. Preserved
variants must FLAG; catalytic-destroyed controls must CLEAR — or the script
prints DEMO INVALIDATED on itself. Details: REVIEWER_PROTOCOL.md.

## How it works

    sequence ──► 47 functional fingerprints (4 toxin classes)
                 enzymatic layer: catalytic residues verified against
                 UniProt feature annotations, residue by residue
             ──► gated CF score = identity × (0.50·fp + 0.30·BLOSUM62 + 0.20·AAindex)
                 catalytic identity is a *precondition*, not a weight
             ──► FLAGGED if ≥ 0.60  ·  attribution: catalytic-consistent match only

Detection runs offline from a 55 KB JSON database (`data/product_db/`).
A supplementary ESM-embedding layer and an actuarial pricing module
(0.12% pool loading) live in `src/bioscreen/`.

## Measured evidence (all in results/)

| Claim | Number | Source |
|---|---|---|
| Real AI variants caught (motif preserved, scaffold redesigned) | 36/36 | evodiff_16toxins_v2_preserved.json |
| Robust to conservative window mutation | k90 ≥ 8–14 (median 14, censored) | conservative_substitution_envelope.json |
| Single-codon catalytic adversaries caught | 56/62 (exact map) | adversarial_envelope_final.json |
| Zero false positives on safe proteins | 0 across ~15,000 windows | contextfree_specificity_*.json |
| Database provenance | 12/12 natives match first-principles prediction to 4 decimals | uniprot_catalytic_audit.json |
| Live reproduction of published map | \|Δ\| = 0.00000, deterministic | live_proof_*.json |
| Per-screen latency | 0.01–0.02 s | product_freeze_report.json |

## What BioScreen is NOT (measured, not asserted)

- **Not functional confirmation.** Sequence-level, fold-blind — appropriate
  for DNA-synthesis screening, not wet-lab claims.
- **Documented evasion boundary:** full conservative substitution of the
  catalytic core evades CF — and typically destroys catalysis
  (Y→F removes the reactive hydroxyl; HExxH→KExxR kills metal coordination).
  Single-substitution evaders are mapped exactly in the envelope file.
- Cross-genus variants <40% identity: out of scope.
- Not regulatory-grade; not wet-lab validated.

## Repo layout

    START_BIOSCREEN.bat          the product (double-click)
    scripts/run_product.py       Flask wrapper → http://127.0.0.1:5055
    scripts/live_proof.py        reviewer-driven live validation
    src/bioscreen/kernel.py      screen_sequence(): the entire product brain
    src/bioscreen/contextfree.py gated fingerprint engine
    src/bioscreen/fingerprint.py ESM supplementary layer
    src/bioscreen/pricing.py     actuarial pricing model
    data/product_db/             frozen offline database + demo bank
    results/                     every number above, as JSON
    WHITEPAPER.md                technical detail · REVIEWER_PROTOCOL.md
    DEMO_SCRIPT.md               2-minute judge walkthrough

**Author:** Peter Mburu Kariuki · BioPREVAIL Biosecurity Innovation Challenge submission
