"""
BioScreen — Functional screening with per-residue local alignment.

Fixes the mean-pooling problem by preserving residue-level information.

hazard_score = local_alignment_to_toxins − local_alignment_to_safe_refs

local_alignment(query, reference) = mean of top-K best per-residue matches
where K = 50 (a "motif-sized" window)
"""
import pickle
import numpy as np
from pathlib import Path
import torch


def _local_alignment(query_emb, ref_emb, top_k=50):
    """
    Compute the local alignment score between query and reference.

    For each query residue, find its best match in the reference.
    Return the mean of the top-k best matches.
    """
    # Normalize
    q_norms = np.linalg.norm(query_emb, axis=1, keepdims=True)
    r_norms = np.linalg.norm(ref_emb, axis=1, keepdims=True)
    q = query_emb / (q_norms + 1e-8)
    r = ref_emb / (r_norms + 1e-8)

    # Similarity matrix: (L_query, L_ref)
    sim = q @ r.T

    # For each query residue, take its best match in the reference
    best_per_query_residue = sim.max(axis=1)

    # Top-k mean (motif-sized window)
    k = min(top_k, len(best_per_query_residue))
    top_k_values = np.sort(best_per_query_residue)[-k:]

    return float(top_k_values.mean())


class BioScreenV2:
    """
    Alignment-based screening using per-residue ESM-C embeddings.
    """
    def __init__(self, model, tokenizer, per_residue_path,
                 hazard_threshold=0.02, top_k=50,
                 toxin_similarity_floor=0.90):
        self.model = model
        self.tokenizer = tokenizer
        self.top_k = top_k

        data = pickle.loads(Path(per_residue_path).read_bytes())
        self.toxin_accessions = data["toxins"]["accessions"]
        self.toxin_names = data["toxins"]["names"]
        self.toxin_embeddings = data["toxins"]["embeddings"]  # dict: acc -> (L, D)

        self.safe_accessions = data["safes"]["accessions"]
        self.safe_names = data["safes"]["names"]
        self.safe_embeddings = data["safes"]["embeddings"]

        self.hazard_threshold = hazard_threshold
        self.toxin_similarity_floor = toxin_similarity_floor

    def embed_per_residue(self, sequence):
        inputs = self.tokenizer([sequence.upper().strip()],
                                 return_tensors="pt", padding=True)
        with torch.inference_mode():
            output = self.model(**inputs, output_hidden_states=True)
        return output.hidden_states[-1][0, 1:-1].cpu().numpy()

    def screen(self, sequence):
        query_emb = self.embed_per_residue(sequence)

        # Best local alignment to toxins
        toxin_scores = []
        for acc in self.toxin_accessions:
            ref = self.toxin_embeddings[acc]
            score = _local_alignment(query_emb, ref, top_k=self.top_k)
            toxin_scores.append(score)
        toxin_scores = np.array(toxin_scores)
        best_toxin_idx = int(np.argmax(toxin_scores))
        best_toxin_score = float(toxin_scores[best_toxin_idx])

        # Best local alignment to safe references
        safe_scores = []
        for acc in self.safe_accessions:
            ref = self.safe_embeddings[acc]
            score = _local_alignment(query_emb, ref, top_k=self.top_k)
            safe_scores.append(score)
        safe_scores = np.array(safe_scores)
        best_safe_idx = int(np.argmax(safe_scores))
        best_safe_score = float(safe_scores[best_safe_idx])

        hazard_delta = best_toxin_score - best_safe_score
        flagged = (hazard_delta >= self.hazard_threshold
                   and best_toxin_score >= self.toxin_similarity_floor)

        return {
            "best_toxin_score": round(best_toxin_score, 4),
            "best_toxin_name": self.toxin_names[best_toxin_idx],
            "best_safe_score": round(best_safe_score, 4),
            "best_safe_name": self.safe_names[best_safe_idx],
            "hazard_delta": round(hazard_delta, 4),
            "flagged": flagged,
        }
