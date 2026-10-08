"""
BioScreen — EvoDiff generation wrapper.

Given a toxin accession, generates N EvoDiff variants using either:
  - "preserved" mode: catalytic windows protected (functional variants)
  - "aggressive" mode: entire region redesigned (broken variants)

Each variant is scored with BLAST against the native toxin and with the
ordered fingerprint. Every step is timed and recorded in a `trace` list,
which the Studio renders in the Pipeline trace tab.
"""
import sys
import time
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
# Trace helper
# ----------------------------------------------------------------------
class _Tracer:
    """Small helper to keep step counting and timing in one place."""

    def __init__(self):
        self.steps = []
        self._n = 0

    def add(self, phase, name, elapsed_s, variant, input_str, output_str, details=None):
        self._n += 1
        self.steps.append({
            "step": self._n,
            "phase": phase,
            "variant": variant,
            "name": name,
            "elapsed_ms": round(elapsed_s * 1000, 2),
            "input": input_str,
            "output": output_str,
            "details": details or {},
        })

    def time(self, phase, name, variant, input_str, output_str, fn, details_fn=None):
        """Time a callable and record a trace step."""
        t0 = time.perf_counter()
        result = fn()
        elapsed = time.perf_counter() - t0
        details = details_fn(result) if details_fn else {}
        self.add(phase, name, elapsed, variant, input_str, output_str, details)
        return result


# ----------------------------------------------------------------------
# Core generation
# ----------------------------------------------------------------------
def _region_seq(accession: str) -> str:
    tox_path = ROOT / "data" / "raw" / "toxin_references.csv"
    df = pd.read_csv(tox_path)
    row = df[df["accession"] == accession]
    if row.empty:
        raise ValueError(f"Unknown toxin: {accession}")
    full = row.iloc[0]["sequence"]
    lo, hi = TOXIN_CONFIG[accession]["region"]
    return full[lo:hi]


def _build_chunks(L: int, cat_rel: list, mode: str) -> list:
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
    if accession not in TOXIN_CONFIG:
        raise ValueError(f"Unknown toxin: {accession}")
    cfg = TOXIN_CONFIG[accession]

    tracer = _Tracer()

    # ---- Setup ----
    def _do_setup():
        model_evo, tokenizer_evo = _load_evo()
        region = _region_seq(accession)
        return model_evo, tokenizer_evo, region

    model_evo, tokenizer_evo, region = tracer.time(
        "setup", "Load EvoDiff model + region sequence", None,
        accession, "", _do_setup,
        details_fn=lambda r: {
            "region_length": len(r[2]),
            "model": "OA_DM_38M",
        },
    )

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

    # Setup trace: chunk layout
    tracer.add(
        "setup", "Configure redesign chunks", 0.0, None,
        f"mode={mode}, {len(cat_rel)} catalytic residues",
        f"{len(chunks)} chunks: {chunks}",
        details={
            "mode": mode,
            "catalytic_residues_relative": cat_rel,
            "chunks": chunks,
            "total_redesign_residues": sum(b - a for a, b in chunks),
        },
    )

    # ---- Build own fingerprint ----
    own_fp = None
    if embed_fn is not None:
        region_emb = tracer.time(
            "fingerprint", "Embed region for fingerprint extraction", None,
            f"{L} residues", "", lambda: embed_fn(region),
            details_fn=lambda e: {"embedding_shape": list(e.shape)},
        )
        if fp_idx:
            own_fp = region_emb[fp_idx]
            tracer.add(
                "fingerprint", "Extract fingerprint from catalytic windows", 0.0, None,
                f"fingerprint indices: {len(fp_idx)} residues",
                f"fingerprint shape: {own_fp.shape}",
                details={"n_fingerprint_residues": len(fp_idx)},
            )
        else:
            if fingerprint_embeddings and accession in fingerprint_embeddings:
                own_fp = fingerprint_embeddings[accession]
                tracer.add(
                    "fingerprint", "Reuse reference fingerprint (aggressive mode)", 0.0, None,
                    accession, f"shape: {own_fp.shape}", {},
                )

    # ---- Generate variants ----
    results = []
    for i in range(n_variants):
        if progress_cb:
            progress_cb(f"variant {i+1}/{n_variants}", i / max(n_variants, 1))

        torch.manual_seed(abs(hash(accession)) % 100000 + i * 7)
        np.random.seed(abs(hash(accession)) % 100000 + i * 7)

        # --- EvoDiff generation, per chunk ---
        current = region
        for ci, (a, b) in enumerate(chunks):
            def _do_inpaint():
                from evodiff.conditional_generation import inpaint_simple
                raw = inpaint_simple(
                    model=model_evo, sequence=current,
                    start_idx=a, end_idx=b,
                    tokenizer=tokenizer_evo, device="cpu",
                )
                return raw[1]

            next_seq = tracer.time(
                "evodiff", f"EvoDiff inpaint — chunk {ci+1}/{len(chunks)}", i,
                f"redesign residues [{a}:{b}]",
                "", _do_inpaint,
                details_fn=lambda s: {
                    "chunk_start": a,
                    "chunk_end": b,
                    "chunk_length": b - a,
                    "output_length": len(s),
                },
            )
            current = next_seq
        variant = current

        # --- Identity ---
        def _do_identity():
            return SequenceMatcher(None, variant, region).ratio()
        identity = tracer.time(
            "score", "Sequence identity to native", i,
            f"{len(variant)} residues", "", _do_identity,
            details_fn=lambda v: {"identity_pct": round(v * 100, 2)},
        )

        # --- BLAST ---
        def _do_blast():
            return blast_score(variant, region)
        bl = tracer.time(
            "score", "BLAST local alignment", i,
            "variant vs native region", "", _do_blast,
            details_fn=lambda v: {
                "normalized_score": round(v, 4),
                "threshold": 0.30,
                "verdict": "EVADED" if v < 0.30 else "caught",
            },
        )

        # --- Fingerprint ---
        fp_own = None
        fp_ranking = []
        if own_fp is not None and embed_fn is not None:
            def _do_embed():
                return embed_fn(variant)
            v_emb = tracer.time(
                "embed", "ESM-C forward pass (variant)", i,
                f"{len(variant)} residues", "", _do_embed,
                details_fn=lambda e: {"shape": list(e.shape)},
            )

            def _do_fp():
                q = v_emb / (np.linalg.norm(v_emb, axis=1, keepdims=True) + 1e-8)
                f = own_fp / (np.linalg.norm(own_fp, axis=1, keepdims=True) + 1e-8)
                sims = f @ q.T
                pos = sims.argmax(axis=1)
                val = sims.max(axis=1)
                order = float(np.mean(pos[1:] > pos[:-1])) if len(pos) > 1 else 0.0
                k = min(20, len(val))
                top_k = float(np.sort(val)[-k:].mean())
                return top_k * order, top_k, order
            fp_own, fp_topk, fp_order = tracer.time(
                "fingerprint", "Ordered fingerprint score", i,
                f"fingerprint vs variant", "", _do_fp,
                details_fn=lambda r: {
                    "combined": round(r[0], 4),
                    "top_k_similarity": round(r[1], 4),
                    "order_consistency": round(r[2], 4),
                    "threshold": 0.85,
                    "verdict": "CAUGHT" if r[0] >= 0.85 else "clear",
                },
            )

            # --- Cross-match ---
            if fingerprint_embeddings:
                def _do_cross():
                    out = []
                    q = v_emb / (np.linalg.norm(v_emb, axis=1, keepdims=True) + 1e-8)
                    for acc, ref_fp in fingerprint_embeddings.items():
                        if acc == accession:
                            continue
                        f2 = ref_fp / (np.linalg.norm(ref_fp, axis=1, keepdims=True) + 1e-8)
                        s2 = f2 @ q.T
                        p2 = s2.argmax(axis=1)
                        v2 = s2.max(axis=1)
                        o2 = float(np.mean(p2[1:] > p2[:-1])) if len(p2) > 1 else 0.0
                        k2 = min(20, len(v2))
                        out.append((acc, float(np.sort(v2)[-k2:].mean()) * o2))
                    return out
                fp_ranking = tracer.time(
                    "cross-match", "Cross-match vs 15 other fingerprints", i,
                    "", "", _do_cross,
                    details_fn=lambda r: {"n_compared": len(r)},
                )

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

    # Final verdict per variant
    for i, r in enumerate(results):
        tracer.add(
            "verdict", f"Variant v{i} verdict", 0.0, i,
            f"BLAST={r['blast']:.3f}, FP={r['fingerprint']:.3f}" if r['fingerprint'] is not None else f"BLAST={r['blast']:.3f}",
            f"{'EVADED' if r['blast_evaded'] else 'caught'} by BLAST / "
            f"{'CAUGHT' if r['fingerprint_caught'] else 'clear'} by fingerprint",
            {},
        )

    return {
        "accession": accession,
        "name": cfg["name"],
        "mode": mode,
        "region_length": L,
        "n_catalytic": len(cat_rel),
        "chunks": chunks,
        "variants": results,
        "trace": tracer.steps,
    }
