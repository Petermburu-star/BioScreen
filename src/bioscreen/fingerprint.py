"""
BioScreen — Ordered catalytic fingerprint.

A compact representation of the residues that determine a toxin's chemistry.
Fingerprints are extracted by finding the residues least similar to a
safe-protein background. Queries are scored by both per-residue similarity
and sequence-order consistency.

Two constraints eliminate the false positives that naive top-k similarity
produces on unrelated proteins:

    score = top20_similarity × order_consistency
"""
from __future__ import annotations

import pickle
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import torch


# ----------------------------------------------------------------------
# Standalone scoring
# ----------------------------------------------------------------------

def extract_fingerprint(
    toxin_embeddings: np.ndarray,
    safe_pool: np.ndarray,
    n_residues: int = 42,
) -> np.ndarray:
    """
    Return the indices of the n_residues that differ most from the safe pool.

    toxin_embeddings : (L_toxin, D)
    safe_pool        : (N_safe_residues, D)
    """
    t = toxin_embeddings / (np.linalg.norm(toxin_embeddings, axis=1, keepdims=True) + 1e-8)
    s = safe_pool / (np.linalg.norm(safe_pool, axis=1, keepdims=True) + 1e-8)
    max_sim = (t @ s.T).max(axis=1)
    idx = np.argsort(max_sim)[:n_residues]
    return np.sort(idx)


def score_fingerprint(
    query_embeddings: np.ndarray,
    fingerprint_embeddings: np.ndarray,
    top_k: int = 20,
) -> Dict[str, float]:
    """
    Score a query against a fingerprint.

    query_embeddings       : (L_query, D)
    fingerprint_embeddings : (N_fp, D)

    Returns:
        combined         — top20_similarity × order_consistency
        top_k_similarity — mean of the top-k best per-fingerprint-residue matches
        order_consistency — fraction of consecutive fingerprint residues whose
                            best matches appear in increasing query position
    """
    q = query_embeddings / (np.linalg.norm(query_embeddings, axis=1, keepdims=True) + 1e-8)
    f = fingerprint_embeddings / (np.linalg.norm(fingerprint_embeddings, axis=1, keepdims=True) + 1e-8)

    sims = f @ q.T                          # (N_fp, L_query)
    best_pos = sims.argmax(axis=1)          # (N_fp,)
    best_val = sims.max(axis=1)             # (N_fp,)

    if len(best_pos) < 2:
        order = 0.0
    else:
        order = float(np.mean(best_pos[1:] > best_pos[:-1]))

    k = min(top_k, len(best_val))
    top_k_sim = float(np.sort(best_val)[-k:].mean())

    return {
        "combined": top_k_sim * order,
        "top_k_similarity": top_k_sim,
        "order_consistency": order,
    }


# ----------------------------------------------------------------------
# FingerprintEngine
# ----------------------------------------------------------------------

class FingerprintEngine:
    """
    Screening engine that uses ordered catalytic fingerprints.

    Loads a fingerprint pickle produced by scripts/build_fingerprints.py.

    Usage:
        engine = FingerprintEngine(model, tokenizer, fingerprint_path)
        result = engine.screen(sequence)
    """

    THRESHOLD = 0.85
    TOP_K = 20

    def __init__(
        self,
        model,
        tokenizer,
        fingerprint_path: Path,
        threshold: float = THRESHOLD,
        top_k: int = TOP_K,
    ):
        self.model = model
        self.tokenizer = tokenizer
        self.threshold = threshold
        self.top_k = top_k

        data = pickle.loads(Path(fingerprint_path).read_bytes())
        self.fingerprints: Dict[str, np.ndarray] = data["fingerprints"]
        self.names: Dict[str, str] = data["names"]
        self.n_residues: int = data.get("n_residues", 42)

    # ----------------------------------------------------------------
    def embed(self, sequence: str) -> np.ndarray:
        inputs = self.tokenizer([sequence.upper().strip()], return_tensors="pt", padding=True)
        with torch.inference_mode():
            out = self.model(**inputs, output_hidden_states=True)
        return out.hidden_states[-1][0, 1:-1].cpu().numpy()

    # ----------------------------------------------------------------
    def screen(self, sequence: str) -> Dict:
        """Score a sequence against every fingerprint. Returns the best match."""
        query = self.embed(sequence)

        scores = {}
        for acc, fp in self.fingerprints.items():
            scores[acc] = score_fingerprint(query, fp, top_k=self.top_k)

        best_acc = max(scores, key=lambda a: scores[a]["combined"])
        best = scores[best_acc]

        # Sort others for display
        ranked = sorted(
            [(acc, s["combined"]) for acc, s in scores.items()],
            key=lambda kv: -kv[1],
        )

        return {
            "fingerprint_score": round(best["combined"], 4),
            "fingerprint_top_k": round(best["top_k_similarity"], 4),
            "fingerprint_order": round(best["order_consistency"], 4),
            "fingerprint_best": self.names[best_acc],
            "fingerprint_best_accession": best_acc,
            "fingerprint_flagged": bool(best["combined"] >= self.threshold),
            "fingerprint_threshold": self.threshold,
            "fingerprint_ranking": [
                {"accession": acc, "name": self.names[acc], "score": round(s, 4)}
                for acc, s in ranked[:5]
            ],
        }
