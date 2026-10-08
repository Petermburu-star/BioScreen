"""
BioScreen — EvoDiff generation wrapper.

Given a toxin accession, generates N EvoDiff variants using either:
  - "preserved" mode: catalytic windows protected (functional variants)
  - "aggressive" mode: entire region redesigned (broken variants)

Each variant is scored with BLAST against the native toxin and with the
ordered fingerprint. Returns results ready for the Studio to display.
"""
import sys
from pathlib import Path
from difflib import SequenceMatcher
from typing import Callable, Optional

ROOT = Path(__file__).resolve().parents[2]
EVO_DIR = ROOT / "external" / "evodiff"
if str(EVO_DIR) not in sys.path:
    sys.path.insert(0, str(EVO_DIR))

import numpy as np
import pandas as pd
import torch

from Bio.Align import PairwiseAligner, substitution_matrices

from bioscreen.toxins import TOXIN_CONFIG


# ----------------------------------------------------------------------
# Lazy model loading
# ----------------------------------------------------------------------
_evo = {"model": None, "tokenizer": None, "loaded": False}


def _load_evo():
    if _evo["loaded"]:
        return _evo["model"], _evo["tokenizer"]
    from evodiff.pretrained import OA_DM_38M
    model, _, tokenizer, _ = OA_DM_38M()
    model.eval()
    _evo["model"] = model
    _evo["tokenizer"] = tokenizer
    _evo["loaded"] = True
    return model, tokenizer


# ----------------------------------------------------------------------
# BLAST
# ----------------------------------------------------------------------
_blosum = substitution_matrices.load("BLOSUM62")
_aligner = PairwiseAligner()
_aligner.substitution_matrix = _blosum
_aligner.open_gap_score = -11
_aligner.extend_gap_score = -1
_aligner.mode = "local"


def blast_score(query: str, reference: str) -> float:
    try:
        return _aligner.score(query, reference) / _aligner.score(reference, reference)
    except Exception:
        return 0.0


# ----------------------------------------------------------------------
# Core generation
# ----------------------------------------------------------------------
def _region_seq(accession: str) -> str:
    """Load the sequence of a toxin's functional region from the raw CSV."""
    tox_path = ROOT / "data" / "raw" / "toxin_references.csv"
    df = pd.read_csv(tox_path)
    row = df[df["accession"] == accession]
    if row.empty:
        raise ValueError(f"Unknown toxin: {accession}")
    full = row.iloc[0]["sequence"]
    lo, hi = TOXIN_CONFIG[accession]["region"]
    return full[lo:hi]


def _build_chunks(L: int, cat_rel: list, mode: str) -> list:
    """
    Returns list of (start, end) redesign chunks.
    preserved mode: chunks avoid the catalytic windows (±8)
    aggressive mode: single chunk covering almost the whole region
    """
    if mode == "aggressive":
        return [(5, L - 5)] if L > 15 else []

    protected = set()
    for p in cat_rel:
        for off in range(-8, 9):
            if 0 <= p + off < L:
                protected.add(p + off)
    chunks = []
    prev = 0
    for p in sorted(protected):
        if p > prev:
            chunks.append((prev, p))
        prev = p + 1
    if prev < L:
        chunks.append((prev, L))
    return [(a, b) for a, b in chunks if b - a >= 10]


def generate_variants(
    accession: str,
    mode: str = "preserved",
    n_variants: int = 3,
    fingerprint_embeddings: Optional[dict] = None,
    embed_fn: Optional[Callable] = None,
    progress_cb: Optional[Callable[[str, float], None]] = None,
) -> dict:
    """
    Generate N variants of a toxin. Returns scores for each.

    fingerprint_embeddings: {accession: ndarray(42, 1152)}
    embed_fn: callable(sequence) -> ndarray(L, 1152)
    progress_cb: optional callback(stage, fraction_done)
    """
    if accession not in TOXIN_CONFIG:
        raise ValueError(f"Unknown toxin: {accession}")
    cfg = TOXIN_CONFIG[accession]

    model_evo, tokenizer_evo = _load_evo()
    region = _region_seq(accession)
    L = len(region)

    lo, hi = cfg["region"]
    cat_rel = [c - lo for c in cfg["catalytic"] if lo <= c < hi]

    if mode == "preserved":
        chunks = _build_chunks(L, cat_rel, "preserved")
        fp_idx = sorted(set(
            p + off for p in cat_rel for off in range(-6, 7)
            if 0 <= p + off < L
        ))
    else:
        chunks = _build_chunks(L, cat_rel, "aggressive")
        fp_idx = None

    if not chunks:
        raise ValueError(f"No redesign region for {accession} in mode {mode}")

    # Build fingerprint for this toxin (fresh, isolated context)
    own_fp = None
    if embed_fn is not None:
        region_emb = embed_fn(region)
        if fp_idx:
            own_fp = region_emb[fp_idx]
        else:
            if fingerprint_embeddings and accession in fingerprint_embeddings:
                own_fp = fingerprint_embeddings[accession]

    results = []
    for i in range(n_variants):
        if progress_cb:
            progress_cb(f"variant {i+1}/{n_variants}", i / max(n_variants, 1))

        torch.manual_seed(abs(hash(accession)) % 100000 + i * 7)
        np.random.seed(abs(hash(accession)) % 100000 + i * 7)

        current = region
        for (a, b) in chunks:
            from evodiff.conditional_generation import inpaint_simple
            raw = inpaint_simple(
                model=model_evo, sequence=current,
                start_idx=a, end_idx=b,
                tokenizer=tokenizer_evo, device="cpu",
            )
            current = raw[1]
        variant = current

        identity = SequenceMatcher(None, variant, region).ratio()
        bl = blast_score(variant, region)

        fp_own = None
        fp_ranking = []
        if own_fp is not None and embed_fn is not None:
            v_emb = embed_fn(variant)
            q = v_emb / (np.linalg.norm(v_emb, axis=1, keepdims=True) + 1e-8)
            f = own_fp / (np.linalg.norm(own_fp, axis=1, keepdims=True) + 1e-8)
            sims = f @ q.T
            pos = sims.argmax(axis=1)
            val = sims.max(axis=1)
            order = float(np.mean(pos[1:] > pos[:-1])) if len(pos) > 1 else 0.0
            k = min(20, len(val))
            fp_own = float(np.sort(val)[-k:].mean()) * order

            # Cross-match against all other fingerprints
            if fingerprint_embeddings:
                for acc, ref_fp in fingerprint_embeddings.items():
                    if acc == accession:
                        continue
                    f2 = ref_fp / (np.linalg.norm(ref_fp, axis=1, keepdims=True) + 1e-8)
                    s2 = f2 @ q.T
                    p2 = s2.argmax(axis=1)
                    v2 = s2.max(axis=1)
                    o2 = float(np.mean(p2[1:] > p2[:-1])) if len(p2) > 1 else 0.0
                    k2 = min(20, len(v2))
                    fp_ranking.append((acc, float(np.sort(v2)[-k2:].mean()) * o2))

        results.append({
            "idx": i,
            "sequence": variant,
            "length": len(variant),
            "identity": round(identity, 4),
            "blast": round(bl, 4),
            "blast_evaded": bl < 0.30,
            "fingerprint": round(fp_own, 4) if fp_own is not None else None,
            "fingerprint_caught": (fp_own is not None and fp_own >= 0.85),
            "cross_ranking": sorted(fp_ranking, key=lambda x: -x[1])[:5] if fp_ranking else [],
        })

    if progress_cb:
        progress_cb("done", 1.0)

    return {
        "accession": accession,
        "name": cfg["name"],
        "mode": mode,
        "region_length": L,
        "n_catalytic": len(cat_rel),
        "chunks": chunks,
        "variants": results,
    }
