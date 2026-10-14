# Claim -> Evidence -> Demo map

| Pitch claim | Evidence (results/) | Demo button |
|---|---|---|
| Catches AI-designed variants of known toxins | evodiff_16toxins_v2_preserved.json (36/36) | ai_real_* |
| Survives scaffold redesign | adversarial_envelope_final.json Axis A | ai_real_* |
| Catches single-codon adversaries | adversarial_envelope_final.json Axis B (56/62) | ai_singlecodon |
| Zero false positives | contextfree_specificity_*.json (~15k windows) | safe_* |
| Knows its own limits | uniprot_catalytic_audit.json, adv_boundary | adv_boundary |
| Laptop-scale, offline | product_freeze_report.json | every click |

**Do not claim:** functional detection (fold-blind), cross-genus <40%,
regulatory-grade, wet-lab validated. Modeccin excluded (unverified config).
BoNT A/E attributed at family level (identical catalytic windows).
