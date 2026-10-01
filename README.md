# BioScreen

**Functional screening and actuarial risk pricing for DNA synthesis biosecurity.**

## The Problem

Current DNA synthesis screening relies on sequence similarity to known hazards.
AI-designed proteins can retain function while evading sequence-based detection —
up to 100% evasion has been demonstrated for some AI-generated variants of
controlled proteins (Microsoft red-teaming study, 2025).

No practical, deployed system uses function-based screening.
No economic mechanism prices biosecurity risk to incentivise transparency.

## The Innovation

BioScreen has two layers:

**Layer 1 — Functional screening.** Embed query protein sequences using ESM-C
(EvolutionaryScale's protein language model). Compare embeddings to a curated
reference database of controlled toxins. Detect functional similarity even when
sequence similarity is low.

**Layer 2 — Actuarial risk pricing.** Convert screening signals and order context
into a risk score. Price biosecurity risk using catastrophe-insurance methods.
Create incentive-compatible premiums that reward transparency.

## Why This Matters

- Catches AI-designed evasion attempts that sequence screening misses
- Creates economic incentives for labs to submit to enhanced screening
- Subsidises screening for resource-limited settings (LMIC labs)
- Scales to pandemic-response conditions (100 Days Mission)

## Status

Working prototype. Every result computed from live predictions — no hardcoded data.

## Stack

- ESM-C 600M (protein language model, MIT license, Hugging Face)
- UniProt toxin reference database
- Python: transformers, torch, numpy, pandas, scikit-learn

## Author

Peter Mburu Kariuki — pmburu346@gmail.com

## Status

Pre-submission prototype for BioPREVAIL × CEPI Biosecurity Innovation Challenge 2026-27.