"""
BioScreen — functional residue acquisition.

Fetches what an AI-designed variant must preserve for the toxin to remain
toxic, from three open sources:

    1. UniProt REST  — active sites, binding sites, metal ligands
    2. InterPro API  — Pfam domain boundaries, full Stockholm alignments
    3. pyhmmer       — HMM built from the alignment, per-position conservation

Everything is cached under data/reference/ so repeated calls are offline.

No torch, no ESM. Only requests, pyhmmer, numpy.

See docs/HANDOFF.md session 4 for the discovery notes that led here:
the Pfam *seed* HMMs (nseq <= 17) are too sparse for conservation work.
The *full* alignments, exposed as `annotation=alignment:full`, carry
hundreds to thousands of sequences and produce usable signal.
"""
from __future__ import annotations

import gzip
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import requests

ROOT = Path(__file__).resolve().parents[2]
CACHE_DIR = ROOT / "data" / "reference"
CACHE_DIR.mkdir(parents=True, exist_ok=True)

UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"
INTERPRO_BASE = "https://www.ebi.ac.uk/interpro/api/entry/pfam"
INTERPRO_HMM = "https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam"

# --------------------------------------------------------------------------
# Data classes
# --------------------------------------------------------------------------
@dataclass
class FunctionalResidues:
    """A toxin's functionally obligate residues, grouped by role."""
    accession: str
    sequence: str
    active_site: list[int] = field(default_factory=list)
    binding_site: list[int] = field(default_factory=list)
    metal_binding: list[int] = field(default_factory=list)
    site: list[int] = field(default_factory=list)

    @property
    def all_positions(self) -> list[int]:
        """Sorted, deduplicated union of every functional position."""
        s = set(self.active_site) | set(self.binding_site) \
            | set(self.metal_binding) | set(self.site)
        return sorted(s)


@dataclass
class PfamDomain:
    pfam_id: str
    name: str
    start: int
    end: int


# --------------------------------------------------------------------------
# UniProt
# --------------------------------------------------------------------------
def fetch_uniprot_sequence(accession: str) -> str:
    """Fetch a UniProt entry's sequence as a plain string."""
    cache = CACHE_DIR / f"seq_{accession}.txt"
    if cache.exists():
        return cache.read_text(encoding="utf-8").strip()

    r = requests.get(f"{UNIPROT_BASE}/{accession}.fasta", timeout=30)
    r.raise_for_status()
    lines = r.text.strip().split("\n")
    seq = "".join(lines[1:])
    cache.write_text(seq, encoding="utf-8")
    return seq


def fetch_uniprot_features(accession: str) -> FunctionalResidues:
    """Fetch a UniProt entry's functional-site features.

    Returns active site, binding site, metal-binding, and generic site
    positions (1-based, in UniProt coordinates).
    """
    cache = CACHE_DIR / f"features_{accession}.json"
    if cache.exists():
        data = json.loads(cache.read_text(encoding="utf-8"))
    else:
        r = requests.get(
            f"{UNIPROT_BASE}/{accession}.json",
            timeout=30,
            headers={"Accept": "application/json"},
        )
        r.raise_for_status()
        data = r.json()
        cache.write_text(json.dumps(data), encoding="utf-8")

    seq = data["sequence"]["value"]
    out = FunctionalResidues(accession=accession, sequence=seq)

    for feat in data.get("features", []):
        ft = feat.get("type", "")
        loc = feat.get("location", {})
        try:
            start = int(loc["start"]["value"])
            end = int(loc["end"]["value"])
        except (KeyError, TypeError):
            continue
        positions = list(range(start, end + 1))
        if ft == "Active site":
            out.active_site.extend(positions)
        elif ft == "Binding site":
            out.binding_site.extend(positions)
        elif ft == "Metal binding":
            out.metal_binding.extend(positions)
        elif ft == "Site":
            out.site.extend(positions)

    # dedupe preserving order
    for attr in ("active_site", "binding_site", "metal_binding", "site"):
        seen = set()
        vals = []
        for v in getattr(out, attr):
            if v not in seen:
                seen.add(v)
                vals.append(v)
        setattr(out, attr, vals)
    return out


# --------------------------------------------------------------------------
# InterPro / Pfam
# --------------------------------------------------------------------------
def fetch_pfam_domains(accession: str) -> list[PfamDomain]:
    """Return Pfam domains on a UniProt entry."""
    cache = CACHE_DIR / f"pfam_{accession}.json"
    if cache.exists():
        data = json.loads(cache.read_text(encoding="utf-8"))
    else:
        r = requests.get(
            f"{INTERPRO_BASE}/protein/uniprot/{accession}/",
            timeout=30,
            headers={"Accept": "application/json"},
        )
        r.raise_for_status()
        data = r.json()
        cache.write_text(json.dumps(data), encoding="utf-8")

    out: list[PfamDomain] = []
    for entry in data.get("results", []):
        md = entry.get("metadata", {})
        acc = md.get("accession")
        nm = md.get("name", {})
        name = nm.get("name") if isinstance(nm, dict) else str(nm)
        for prot in entry.get("proteins", []):
            for loc in prot.get("entry_protein_locations", []):
                for frag in loc.get("fragments", []):
                    out.append(PfamDomain(
                        pfam_id=acc,
                        name=name or "",
                        start=int(frag["start"]),
                        end=int(frag["end"]),
                    ))
                break
            break
    return out


def fetch_pfam_alignment(pfam_id: str) -> Path:
    """Download the full Stockholm alignment for a Pfam family.

    Cached at data/reference/{pfam_id}.stockholm. Returns the path.
    """
    out = CACHE_DIR / f"{pfam_id}.stockholm"
    if out.exists() and out.stat().st_size > 5_000:
        return out

    url = f"{INTERPRO_HMM}/{pfam_id}?annotation=alignment:full"
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    if len(r.content) < 5_000:
        raise RuntimeError(f"{pfam_id}: alignment too small ({len(r.content)} B)")
    out.write_bytes(r.content)
    return out


def fetch_pfam_hmm(pfam_id: str) -> Path:
    """Download the Pfam *seed* HMM (gzipped) and decompress it.

    NOTE: seed HMMs are sparse (nseq <= 17 for most families) and are NOT
    suitable for conservation work. Use build_hmm_from_alignment() instead
    for anything that needs real per-position signal.
    """
    out = CACHE_DIR / f"{pfam_id}.seed.hmm"
    if out.exists():
        return out

    url = f"{INTERPRO_HMM}/{pfam_id}?annotation=hmm"
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    raw = r.content
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    out.write_bytes(raw)
    return out


# --------------------------------------------------------------------------
# HMM construction + conservation
# --------------------------------------------------------------------------
def build_hmm_from_alignment(stockholm_path: Path | str):
    """Build a pyhmmer HMM from a Stockholm alignment.

    Returns the HMM object (not the profile). Uses the full alignment, so
    nseq_effective typically lands in the hundreds — unlike the seed HMMs.
    """
    from pyhmmer.easel import Alphabet, MSAFile
    from pyhmmer.plan7 import Background, Builder

    alphabet = Alphabet.amino()
    with MSAFile(str(stockholm_path), format="stockholm") as f:
        text_msa = next(iter(f))

    if not getattr(text_msa, "name", None):
        text_msa.name = stockholm_path.stem.encode("utf-8")

    digital = text_msa.digitize(alphabet)
    builder = Builder(alphabet)
    bg = Background(alphabet)
    hmm, _profile, _opt = builder.build_msa(digital, background=bg)
    return hmm


def conservation_from_hmm(hmm) -> np.ndarray:
    """Per-match-state conservation, 0 = uniform, 1 = fully conserved.

    Computed as 1 - entropy/max_entropy over the HMM's match emissions.
    Array length equals hmm.M (number of match states).
    """
    em = np.asarray(hmm.match_emissions)
    em_match = em[1:]  # drop the null/begin row
    sums = em_match.sum(axis=1, keepdims=True)
    probs = em_match / np.where(sums > 0, sums, 1.0)
    eps = 1e-12
    entropy = -np.sum(probs * np.log2(probs + eps), axis=1)
    return 1.0 - (entropy / np.log2(20.0))


# --------------------------------------------------------------------------
# Convenience
# --------------------------------------------------------------------------
def fetch_functional_residues(accession: str) -> FunctionalResidues:
    """Alias for fetch_uniprot_features."""
    return fetch_uniprot_features(accession)


__all__ = [
    "FunctionalResidues",
    "PfamDomain",
    "fetch_uniprot_sequence",
    "fetch_uniprot_features",
    "fetch_functional_residues",
    "fetch_pfam_domains",
    "fetch_pfam_alignment",
    "fetch_pfam_hmm",
    "build_hmm_from_alignment",
    "conservation_from_hmm",
]
