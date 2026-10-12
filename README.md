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

### Full EvoDiff benchmark — 13 toxins, 39 variants

| Mode | BLAST evaded | ESM fingerprint | Context-free |
|---|---|---|---|
| **Preserved** | 39/39 (100%) | **20/39 (51.3%)** | **39/39 (100%)** |
| **Aggressive** | 39/39 (100%) | 0/39 | **0/39** |

**Separation on the designed test: 39/39 preserved vs 0/39 aggressive.**

The ESM fingerprint alone catches 51% of preserved-mode variants. The
context-free layer catches all of them. Combined, the two-layer detector
catches everything on this benchmark.

### Specificity

| Panel | Result |
|---|---|
| Safe human proteins (25, ~2,500 windows) | **0 false positives** |
| Synthetic destroyed fingerprint (all residues → Ala) | correctly rejected |

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

**17 controlled toxins** (13 with catalytic fingerprints).

Recently expanded with **P11439 (Pseudomonas exotoxin A)** and corrected 7 of
13 catalytic configurations against UniProt active-site annotations.

---

## What's novel — honestly

- **The ESM ordered catalytic fingerprint** with the sequence-order consistency
  constraint is the primary innovation.
- **The context-free fallback layer** is built from standard bioinformatics
  tools (BLOSUM62, AAindex, catalytic motif patterns). It is a well-engineered
  second signal, not a new algorithm.
- **The finding that ESM fingerprints are context-sensitive** on AI-designed
  variants is a real contribution — a specific diagnosis and a specific fix.
- **The actuarial pricing layer** applying catastrophe-insurance methods to
  DNA synthesis screening is a novel combination.
- **Preliminary frustration analysis** (6.15× enrichment of catalytic
  glutamates in top-5% frustration) is documented as a promising direction,
  not yet integrated.

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
