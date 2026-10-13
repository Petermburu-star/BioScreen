# BioScreen

**Two-layer screening and actuarial risk pricing for DNA synthesis biosecurity.**

A working detector for AI-designed toxin variants that sequence screening cannot
see — validated against Microsoft's own EvoDiff model — plus a self-financing
economic mechanism that funds global biosecurity capacity.

---

## The problem

Current DNA synthesis screening relies on sequence similarity to known hazards.
AI protein design can now rewrite a toxin's primary sequence while preserving
its function. A 2025 Microsoft-led red-teaming study demonstrated that **up to
100% of AI-designed ricin variants evaded major commercial screening tools**.

As of October 2026, no working countermeasure has been deployed.

---

## The solution

BioScreen runs **two independent detection layers**.

**Layer 1 — ESM ordered catalytic fingerprint.** ESM-C 600M per-residue
embeddings, scored at the residues surrounding each toxin's catalytic sites,
with a sequence-order consistency constraint that rejects safe proteins.

**Layer 2 — Context-free fingerprint.** A fallback layer that does not depend
on any protein language model. Scores four signals:

- residue identity at catalytic positions
- BLOSUM62 substitution scores
- AAindex physicochemical profile correlation over the fingerprint window
- mechanism-specific catalytic motifs (HExxH, YYER, ER, etc.)

A query is flagged if **either** layer fires.

**Layer 3 (pricing).** Converts the screening signal and order context into
a risk-priced premium. A 0.12% loading on synthesis orders funds a global
biosecurity pool.

---

## Validated results

### Full EvoDiff benchmark — 13 enzymatic toxins, 39 variants

| Mode | BLAST evaded | ESM fingerprint | Context-free |
|---|---|---|---|
| **Preserved** | 39/39 (100%) | **20/39 (51.3%)** | **39/39 (100%)** |
| **Aggressive** | 39/39 (100%) | 0/39 | **0/39** |

**Separation on the designed test: 39/39 preserved vs 0/39 aggressive.**

**EvoDiff benchmark for superantigen and pore-forming classes is running.**

The ESM fingerprint alone catches 51% of preserved-mode variants. The
context-free layer catches all of them. Combined, the two-layer detector
catches everything on this benchmark.

### Specificity across all 48 toxins

| Test | Result |
|---|---|
| Length-unbiased panel (26 fingerprints × 25 safes, 150-aa windows) | **0 / 9,984** |
| Original 60-protein panel (13 enzymatic, full region) | **0 / 8,154** |
| CDC specificity (8 fingerprints × 25 safes) | **0 / 717** |
| Alpha-hemolysin specificity (8 fingerprints × 25 safes) | **0 / 2,019** |
| Shiga B specificity (7 fingerprints × 25 safes) | **0 / 3,815** |
| Cross-class off-diagonal (4 classes × 4) | **0 / 6** |

### Why we added the context-free layer

The ESM fingerprint has a **context-dependence limitation**. A P11439 variant
preserved all 50 fingerprint residues verbatim but scored 0.6748 — below the
0.85 threshold — because redesigning the 158 surrounding residues shifted
ESM-C's embedding at every fingerprint position. Preserving residues does not
guarantee preserving the embedding neighborhood.

The context-free layer does not have this failure mode. It scores residue
identity and physicochemical properties directly, which are invariant to
scaffold redesign.

### Reference set

**48 controlled toxins across 4 functional classes.**

| Class | Count | Families |
|---|---|---|
| **Enzymatic** | 13 | Ricin, botulinum A/B/E, tetanus, shiga 1A/2A, diphtheria, abrin, modeccin, anthrax LF/EF, Pseudomonas ExoA |
| **Superantigen** | 12 | Staphylococcal enterotoxins (SEA–SEI, TS3T-1) + streptococcal pyrogenic exotoxins (SPEA, SPEG, SPEJ) |
| **Pore-forming** | 16 | Cholesterol-dependent cytolysins (SLO, PFO, LLO, pneumolysin, ILY) + α-hemolysin family (Hla, LukF, LukS, LukD, LukE, HlgA, HlgB, HlgC) |
| **Binding-only** | 7 | Shiga toxin B subunits (receptor-binding, no catalytic activity) |

Every class uses the same CF algorithm. Only the functional residue library differs.

**Discovery method:** Functional residues for superantigen, pore-forming, and
binding-only classes were discovered by **unsupervised conservation** across
each family's aligned members. The method found:
- CDC undecapeptide `ECTGLAWEWWR` (membrane-inserting motif) at 11/11 across 7/8 members
- Shiga B receptor-binding pocket (7/8 literature residues)
- Alpha-hemolysin subfamily split (Hla vs bicomponent leukocidins)

Validated against published mutagenesis. No manual residue annotation for the new classes.

---

## What's novel — honestly

- **The ESM ordered catalytic fingerprint** with the sequence-order consistency
  constraint is the primary innovation.
- **The context-free fallback layer** is built from standard bioinformatics
  tools (BLOSUM62, AAindex, catalytic motif patterns). It is a well-engineered
  second signal, not a new algorithm.
- **Unsupervised discovery of functional residues across 4 toxin classes.**
  The same conservation analysis works for catalytic, binding-interface,
  pore-forming, and receptor-binding residues. Validated against published
  mutagenesis for three independent families.
- **The finding that ESM fingerprints are context-sensitive** on AI-designed
  variants is a real contribution — a specific diagnosis and a specific fix.
- **The actuarial pricing layer** applying catastrophe-insurance methods to
  DNA synthesis screening is a novel combination.
- **Preliminary frustration analysis** (catalytic glutamate 5.24× enrichment
  in top-5% frustration) documented but not integrated.

---

## What BioScreen is not

- **Not a replacement for BLAST.** Both layers should run in parallel.
- **Not a detector of novel toxins** with no relative in the reference set.
- **Not regulatory-grade.** It is a proof-of-concept for accelerator work.
- **Not a functional detector.** It flags sequences that preserve functional
  signatures. It does not verify 3D geometry. Misfolded proteins with preserved
  sequences will false-positive.
- **Not a cross-genus detector.** Verified against three alignment methods —
  cross-genus variants at <40% identity are outside CF's design. See
  `docs/PROTOCOL.md` §10 for the boundary documentation.
- **Not validated on adversarial variants** designed to preserve residues but
  break geometry.

---

## Install

```bash
git clone https://github.com/Petermburu-star/BioScreen
cd BioScreen
python -m venv .venv
.venv\Scripts\activate           # Windows
# source .venv/bin/activate      # macOS / Linux

pip install -r requirements.txt

# EvoDiff is not pip-installable — clone separately:
git clone https://github.com/microsoft/evodiff external/evodiff

echo "HF_TOKEN=hf_your_token_here" > .env
