"""
BioScreen — EvoDiff generation wrapper.

Two entry points:

  generate_variants(accession, mode, n_variants, ...)
      For the 12 curated toxins. Uses pre-configured regions and catalytic residues.

  generate_variants_custom(sequence, region_start, region_end,
                            catalytic_residues, auto_catalytic, mode,
                            n_variants, safe_pool, ...)
      For arbitrary user-provided sequences. Auto-detects catalytic residues
      if not supplied.

Both return the same result structure. Both instrument every step with trace.
"""
import hashlib
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


def stable_seed(accession: str, variant_index: int = 0) -> int:
    """Deterministic seed from accession string (PEP 456-safe)."""
    h = hashlib.md5(accession.encode("utf-8")).hexdigest()
    base = int(h[:8], 16) % 100000
    return base + variant_index * 7919


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
    def __init__(self):
        self.steps = []
        self._n = 0

    def add(self, phase, name, elapsed_s, variant, input_str, output_str, details=None):
        self._n += 1
        self.steps.append({
            "step": self._n, "phase": phase, "variant": variant,
            "name": name,
            "elapsed_ms": round(elapsed_s * 1000, 2),
            "input": input_str, "output": output_str,
            "details": details or {},
        })

    def time(self, phase, name, variant, input_str, output_str, fn, details_fn=None):
        t0 = time.perf_counter()
        result = fn()
        elapsed = time.perf_counter() - t0
        details = details_fn(result) if details_fn else {}
        self.add(phase, name, elapsed, variant, input_str, output_str, details)
        return result


# ----------------------------------------------------------------------
# Reference data loading (for curated path)
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


# ----------------------------------------------------------------------
# Chunk builder
# ----------------------------------------------------------------------
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


# ----------------------------------------------------------------------
# Auto-detect catalytic residues (fallback for custom sequences)
# ----------------------------------------------------------------------
def _auto_catalytic_positions(region_emb, safe_pool, n_residues=42):
    """
    Fallback: pick the residues least similar to the safe-protein background.
    This is the heuristic used in the earliest fingerprint iterations.
    It won't be as sharp as literature-curated catalytic windows, but it's
    a reasonable default when no annotations exist.
    """
    r = region_emb / (np.linalg.norm(region_emb, axis=1, keepdims=True) + 1e-8)
    s = safe_pool / (np.linalg.norm(safe_pool, axis=1, keepdims=True) + 1e-8)
    max_sim = (r @ s.T).max(axis=1)
    idx = np.argsort(max_sim)[:n_residues]
    return sorted([int(i) for i in idx])


# ----------------------------------------------------------------------
# Shared variant scoring
# ----------------------------------------------------------------------
def _score_variant(variant, region, own_fp, fingerprint_embeddings, embed_fn):
    """Score one variant. Returns (identity, blast, fp_own, cross_ranking)."""
    identity = SequenceMatcher(None, variant, region).ratio()
    bl = blast_score(variant, region)

    fp_own = None
    cross_ranking = []
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

        if fingerprint_embeddings:
            for acc, ref_fp in fingerprint_embeddings.items():
                f2 = ref_fp / (np.linalg.norm(ref_fp, axis=1, keepdims=True) + 1e-8)
                s2 = f2 @ q.T
                p2 = s2.argmax(axis=1)
                v2 = s2.max(axis=1)
                o2 = float(np.mean(p2[1:] > p2[:-1])) if len(p2) > 1 else 0.0
                k2 = min(20, len(v2))
                cross_ranking.append((acc, float(np.sort(v2)[-k2:].mean()) * o2))

    return identity, bl, fp_own, sorted(cross_ranking, key=lambda x: -x[1])[:5]


# ======================================================================
# PATH A — Curated toxin
# ======================================================================
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

    def _do_setup():
        model_evo, tokenizer_evo = _load_evo()
        region = _region_seq(accession)
        return model_evo, tokenizer_evo, region

    model_evo, tokenizer_evo, region = tracer.time(
        "setup", "Load EvoDiff model + region sequence", None,
        accession, "", _do_setup,
        details_fn=lambda r: {"region_length": len(r[2]), "model": "OA_DM_38M"},
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
                f"{len(fp_idx)} residues", f"shape {own_fp.shape}",
                {"n_fingerprint_residues": len(fp_idx)},
            )
        elif fingerprint_embeddings and accession in fingerprint_embeddings:
            own_fp = fingerprint_embeddings[accession]
            tracer.add("fingerprint", "Reuse reference fingerprint (aggressive)",
                       None, accession, f"shape {own_fp.shape}", {})

    results = _run_generation_loop(
        tracer, model_evo, tokenizer_evo, region, chunks,
        own_fp, fingerprint_embeddings, embed_fn, n_variants, progress_cb,
        parent_accession=accession,
    )

    return {
        "source": "curated",
        "accession": accession,
        "name": cfg["name"],
        "mode": mode,
        "region_length": L,
        "n_catalytic": len(cat_rel),
        "chunks": chunks,
        "variants": results,
        "trace": tracer.steps,
    }


# ======================================================================
# PATH B — Custom sequence
# ======================================================================
def generate_variants_custom(
    sequence: str,
    region_start: int = 0,
    region_end: Optional[int] = None,
    catalytic_residues: Optional[list] = None,
    auto_catalytic: bool = True,
    mode: str = "preserved",
    n_variants: int = 3,
    fingerprint_embeddings: Optional[dict] = None,
    embed_fn: Optional[Callable] = None,
    safe_pool: Optional[np.ndarray] = None,
    progress_cb: Optional[Callable[[str, float], None]] = None,
) -> dict:
    """Generate variants of an arbitrary user-supplied protein."""
    seq = (sequence or "").strip().upper().replace(" ", "").replace("\n", "")
    if len(seq) < 50:
        raise ValueError("Sequence must be at least 50 residues")
    if len(seq) > 500:
        seq = seq[:500]

    tracer = _Tracer()

    def _do_setup():
        model_evo, tokenizer_evo = _load_evo()
        return model_evo, tokenizer_evo

    model_evo, tokenizer_evo = tracer.time(
        "setup", "Load EvoDiff model", None, f"{len(seq)} residues", "",
        _do_setup,
        details_fn=lambda r: {"model": "OA_DM_38M", "input_length": len(seq)},
    )

    lo = max(0, region_start)
    hi = len(seq) if region_end is None else min(len(seq), region_end)
    region = seq[lo:hi]
    L = len(region)

    tracer.add(
        "setup", "Configure region and redesign chunks", 0.0, None,
        f"region=[{lo}:{hi}] ({L} residues)",
        "", {"region_start": lo, "region_end": hi, "region_length": L},
    )

    # Determine catalytic residues
    cat_rel = []
    auto_used = False
    if catalytic_residues:
        cat_rel = [c - lo for c in catalytic_residues if lo <= c < hi]
    elif auto_catalytic and embed_fn is not None and safe_pool is not None:
        region_emb = tracer.time(
            "fingerprint", "Embed region for auto-catalytic detection", None,
            f"{L} residues", "", lambda: embed_fn(region),
            details_fn=lambda e: {"embedding_shape": list(e.shape)},
        )
        cat_rel = _auto_catalytic_positions(region_emb, safe_pool, n_residues=42)
        auto_used = True
        tracer.add(
            "fingerprint", "Auto-detect catalytic-like positions", 0.0, None,
            f"top 42 least-safe residues", f"{len(cat_rel)} positions",
            {"n_residues": len(cat_rel), "method": "least_similar_to_safe"},
        )

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
        raise ValueError("No redesign region left after protecting catalytic windows")

    tracer.add(
        "setup", f"Redesign chunks ({mode})", 0.0, None,
        f"{len(chunks)} chunks" if mode == "preserved" else "1 chunk (aggressive)",
        f"{chunks}",
        {"chunks": chunks, "total_redesign_residues": sum(b - a for a, b in chunks)},
    )

    # Build own fingerprint
    own_fp = None
    if embed_fn is not None:
        if 'region_emb' not in dir():
            region_emb = tracer.time(
                "fingerprint", "Embed region for fingerprint extraction", None,
                f"{L} residues", "", lambda: embed_fn(region),
                details_fn=lambda e: {"embedding_shape": list(e.shape)},
            )
        if fp_idx:
            own_fp = region_emb[fp_idx]
            tracer.add(
                "fingerprint", "Extract fingerprint", 0.0, None,
                f"{len(fp_idx)} residues", f"shape {own_fp.shape}",
                {"n_fingerprint_residues": len(fp_idx)},
            )

    results = _run_generation_loop(
        tracer, model_evo, tokenizer_evo, region, chunks,
        own_fp, fingerprint_embeddings, embed_fn, n_variants, progress_cb,
        parent_accession="custom",
    )

    return {
        "source": "custom",
        "accession": "custom",
        "name": f"Custom sequence ({len(seq)} aa)",
        "mode": mode,
        "region_length": L,
        "n_catalytic": len(cat_rel),
        "auto_catalytic": auto_used,
        "chunks": chunks,
        "variants": results,
        "trace": tracer.steps,
    }


# ======================================================================
# Shared generation loop
# ======================================================================
def _run_generation_loop(
    tracer, model_evo, tokenizer_evo, region, chunks,
    own_fp, fingerprint_embeddings, embed_fn, n_variants, progress_cb,
    parent_accession="custom",
):
    results = []
    for i in range(n_variants):
        if progress_cb:
            progress_cb(f"variant {i+1}/{n_variants}", i / max(n_variants, 1))

        seed = stable_seed(parent_accession, i)
        torch.manual_seed(seed)
        np.random.seed(seed)

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

            current = tracer.time(
                "evodiff", f"EvoDiff inpaint — chunk {ci+1}/{len(chunks)}", i,
                f"redesign [{a}:{b}]", "", _do_inpaint,
                details_fn=lambda s: {"chunk_start": a, "chunk_end": b,
                                       "chunk_length": b - a, "output_length": len(s)},
            )
        variant = current

        def _do_score():
            return _score_variant(variant, region, own_fp,
                                   fingerprint_embeddings, embed_fn)
        identity, bl, fp_own, cross_ranking = tracer.time(
            "score", "Identity + BLAST + fingerprint + cross-match", i,
            "variant vs parent + all fingerprints", "", _do_score,
            details_fn=lambda r: {
                "identity_pct": round(r[0] * 100, 2),
                "blast": round(r[1], 4),
                "fingerprint": round(r[2], 4) if r[2] is not None else None,
                "n_cross_matches": len(r[3]),
            },
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
            "cross_ranking": cross_ranking,
        })

        tracer.add(
            "verdict", f"Variant v{i} verdict", 0.0, i,
            f"BLAST={bl:.3f}, FP={fp_own:.3f}" if fp_own is not None else f"BLAST={bl:.3f}",
            f"{'EVADED' if bl < 0.30 else 'caught'} by BLAST / "
            f"{'CAUGHT' if (fp_own or 0) >= 0.85 else 'clear'} by fingerprint",
            {},
        )

    if progress_cb:
        progress_cb("done", 1.0)
    return results
