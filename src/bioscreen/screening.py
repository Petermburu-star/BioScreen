"""
BioScreen — Functional screening with contrastive scoring.

Mean-pooled embeddings don't distinguish toxins from safe proteins.
Fix: hazard_score = max_sim_to_toxins − max_sim_to_safe_reference

A protein is hazardous only if it looks like a toxin AND NOT like a safe
protein.
"""
import pickle
import numpy as np
from pathlib import Path
import torch


class BioScreen:
    def __init__(self, model, tokenizer,
                 toxin_embeddings_path,
                 safe_embeddings_path,
                 hazard_threshold=0.05,
                 toxin_similarity_floor=0.90):
        self.model = model
        self.tokenizer = tokenizer

        toxin_data = pickle.loads(Path(toxin_embeddings_path).read_bytes())
        self.toxin_embeddings = toxin_data["embeddings"]
        self.toxin_names = toxin_data["names"]
        norms = np.linalg.norm(self.toxin_embeddings, axis=1, keepdims=True)
        self.toxin_embeddings_norm = self.toxin_embeddings / (norms + 1e-8)

        safe_data = pickle.loads(Path(safe_embeddings_path).read_bytes())
        self.safe_embeddings = safe_data["embeddings"]
        self.safe_names = safe_data["names"]
        norms = np.linalg.norm(self.safe_embeddings, axis=1, keepdims=True)
        self.safe_embeddings_norm = self.safe_embeddings / (norms + 1e-8)

        self.hazard_threshold = hazard_threshold
        self.toxin_similarity_floor = toxin_similarity_floor

    def embed(self, sequence):
        inputs = self.tokenizer([sequence.upper().strip()],
                                 return_tensors="pt", padding=True)
        with torch.inference_mode():
            output = self.model(**inputs, output_hidden_states=True)
        emb = output.hidden_states[-1]
        return emb[0, 1:-1].mean(dim=0).cpu().numpy()

    def screen(self, sequence):
        query_emb = self.embed(sequence)
        query_norm = query_emb / (np.linalg.norm(query_emb) + 1e-8)

        toxin_sims = self.toxin_embeddings_norm @ query_norm
        best_toxin_idx = int(np.argmax(toxin_sims))
        best_toxin_sim = float(toxin_sims[best_toxin_idx])

        safe_sims = self.safe_embeddings_norm @ query_norm
        best_safe_idx = int(np.argmax(safe_sims))
        best_safe_sim = float(safe_sims[best_safe_idx])

        hazard_delta = best_toxin_sim - best_safe_sim

        flagged = (hazard_delta >= self.hazard_threshold
                   and best_toxin_sim >= self.toxin_similarity_floor)

        return {
            "best_toxin_sim": round(best_toxin_sim, 4),
            "best_toxin_name": self.toxin_names[best_toxin_idx],
            "best_safe_sim": round(best_safe_sim, 4),
            "best_safe_name": self.safe_names[best_safe_idx],
            "hazard_delta": round(hazard_delta, 4),
            "flagged": flagged,
        }
