# BioScreen — Validation Protocols (v2, post-verification)

## 1. Database verification (the gate that replaced trust)
- Fetch UniProt JSON per toxin; read `Active site`/`Binding site` features
  (1-based precursor coordinates).
- Convert to region-relative 0-based as `j = c - 1 - lo` (the c−lo convention
  is a documented failure mode — do not use it).
- Pattern gate per family: RIP `YYER`, zincin `HE..H`, DT `HE`, ExoA `HYYE`,
  Stx `ER`. Failure = exclusion (Modeccin excluded on this basis).
- Arithmetic gate: native self-score must equal
  `0.5 + 0.3·(mean self-BLOSUM+3)/11 + 0.2` to 4 decimals. 12/12 pass.

## 2. Adversarial envelope
- Axis A (window drift): conservative substitutions sampled from ±6 window,
  catalytic positions excluded; 25 draws/k; seed 20261013; threshold 0.60.
- Axis B (single-codon): exact enumeration, every catalytic residue × every
  same-group alternative. 62 combos; results in
  `results/adversarial_envelope_final.json`.
- Full-core: all catalytic residues substituted conservatively; 0/12 detected
  (documented boundary; biology counterweight in WHITEPAPER §5).

## 3. Specificity
- Safe panels (25 + 60 proteins, human), sliding windows, threshold 0.60:
  0 false positives (~15,000 windows). Product path: TNF 0.4745, albumin 0.3769.

## 4. Live proof
- `python scripts/live_proof.py --accession P02879 --seed N --n 1`
- Frame gate → generation → SHA-256 freshness → scored controls →
  PROTOCOL VALID / DEMO INVALIDATED. Mis-wired refusal retained as evidence.

## 5. Known technical notes
- ESM builder (+1 offset on catalytic positions): internally consistent,
  correction scheduled with ESM re-validation.
- Botulinum A/E: identical catalytic windows → family-level attribution (tie reported).
- Custom preserved-mode auto-detection: not a validated feature; product does not ship it.
