# BioScreen — Live Proof Protocol for Reviewers

**Purpose:** prove, in one session, that BioScreen's detection is real —
generated live, scored live, reproducible, consistent with published evidence.

## Run

    python scripts/live_proof.py --accession P02879 --seed <YOUR_NUMBER> --n 2

You choose the toxin and the seed. Generation: a few minutes per variant (CPU).
Everything else is sub-second. An `--with-esm` ESM layer can be added once the
embedding call is adapted to screening_v2's API (see Cell B output).

## Why it is real (not bugs, not bias)

1. **Freshness** — SHA-256 per variant, low identity to native, exists in no
   stored dataset.
2. **Controls** — same seed and generator: preserved variants must FLAG,
   aggressive (catalytic-destroyed) variants must CLEAR. Any failure prints
   DEMO INVALIDATED.
3. **Exact-map reproduction** — hand-apply any substitution from
   results/adversarial_envelope_final.json; the live system returns the
   published score to 4 decimals (verified in Cell B).
4. **Full decomposition** — every CF score prints its components:
   catalytic identity gate, BLOSUM62, AAindex, window identity.
5. **Provenance** — catalytic residues verified against UniProt feature
   annotations (results/uniprot_catalytic_audit.json).

## Known boundaries (volunteered, not hidden)

- Preserved mode protects the catalytic window by instruction. The controls,
  the single-codon test, and the published failure mode (full conservative
  core substitution evades — and typically breaks catalysis) address circularity.
- Sequence-level detection is fold-blind; cross-genus variants below 40%
  identity are out of scope; not wet-lab validated; not regulatory-grade.
