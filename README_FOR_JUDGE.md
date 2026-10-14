# BioScreen — Offline Biosecurity Screening

**Run it:** double-click `START_BIOSCREEN.bat` → browser opens at http://127.0.0.1:5055.
No internet needed. No GPU. Runs on a 16 GB laptop.

**What it does:** paste (or click) a protein sequence → FLAGGED / CLEAR verdict
in under a second, with the matched toxin family, score decomposition, and the
verification provenance of the fingerprint that fired.

**The 10 demo buttons show the full measured envelope:**
- Real AI-designed variants (EvoDiff, motif preserved, scaffold redesigned) → FLAGGED
- AI variants with the catalytic window destroyed (negative control) → CLEAR
- A single-codon adversary (ricin Y>F) → FLAGGED
- A known evasion boundary (diphtheria catalytic H>K) → CLEAR, honestly labeled
- Native toxins → FLAGGED; safe human proteins → CLEAR

**Under the hood (all measured, results/ has the JSON):**
- 47 fingerprints across 4 toxin functional classes; enzymatic layer verified
  residue-by-residue against UniProt feature annotations.
- Gated scoring: catalytic identity is a precondition, not a weight —
  0 false positives across ~15,000 safe-protein windows.
- Adversarial envelope mapped exactly: 36/36 preserved AI variants caught;
  56/62 single catalytic substitutions caught; full-core conservative
  substitution concedes (and typically breaks catalysis).

**Limitations:** sequence-level (fold-blind), cross-genus <40% identity out of
scope, not wet-lab validated, not regulatory-grade. See the in-app panel.
