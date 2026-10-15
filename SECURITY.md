# Security Policy — BioScreen

## Our Commitments

1. **No hazardous sequences are synthesized.** All screening is in silico.
2. **Toxin reference sequences are from public databases** (UniProt).
3. **No customer data is stored in this repository.**
4. **Access control on any deployed version** (RBAC + encryption).

## Ethics Statement

This project develops defensive biosecurity tooling. It does not create,
synthesize, or distribute hazardous biological materials. All validation
is computational. Any deployed version would operate under institutional
biosafety oversight.

## Reporting

Contact: pmburu346@gmail.com

## Dual-use disclosure (biosecurity hygiene)

This repository deliberately does NOT distribute AI-generated toxin-variant
sequences. Evidence files store SHA-256 hashes of every variant: claims are
verifiable by regenerating deterministically via scripts/live_proof.py and
matching the hash. Native reference sequences are UniProt public data with
per-residue provenance. Variant sequences briefly tracked in git history were
removed by history rewrite (this commit range). Documented detection boundary:
conservative substitution of an entire catalytic core evades detection and
typically breaks catalysis (WHITEPAPER.md section 5).
