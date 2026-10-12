"""
BioScreen — context-free fingerprint.

A second detection layer that is independent of protein language models.
Instead of comparing per-residue embeddings (which are context-dependent),
this layer compares:

    1. Residue identity at catalytic positions
    2. BLOSUM62 substitution scores at catalytic positions
    3. AAindex physicochemical profiles over the ±window fingerprint
    4. Residue identity over the whole fingerprint window

All four signals are computed per-residue and are independent of the
surrounding sequence. They cannot shift when a scaffold is redesigned,
which is the failure mode of the ESM-based fingerprint.

See docs/HANDOFF.md session 4 for the P11439 case that motivated this:
ESM scored 0.6748 on a variant that preserved 50/50 fingerprint residues
verbatim. This layer scores the same variant at 0.9773.

No PyTorch, no ESM, no network access. Pure Python + numpy + Bio.BLOSUM62.
"""
from __future__ import annotations

from typing import Iterable, Optional

import numpy as np
from Bio.Align import substitution_matrices


# --------------------------------------------------------------------------
# AAindex physicochemical scales
# --------------------------------------------------------------------------
# Four well-characterised scales, taken from the AAindex1 database
# (Kawashima et al., 2008). Values are verified against the published
# papers, not parsed dynamically — a small, self-contained subset.

HYDROPATHY = {  # KYTJ820101 — Kyte & Doolittle (1982)
    "A": 1.8,  "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8,  "K": -3.9, "M": 1.9,  "F": 2.8,  "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

POLARITY = {  # GRAR740102 — Grantham (1974)
    "A": 8.1,  "R": 10.5, "N": 11.6, "D": 13.0, "C": 5.5,
    "Q": 10.5, "E": 12.3, "G": 9.0,  "H": 10.4, "I": 5.2,
    "L": 4.9,  "K": 11.3, "M": 5.7,  "F": 5.2,  "P": 8.0,
    "S": 9.2,  "T": 8.6,  "W": 5.4,  "Y": 6.2,  "V": 5.9,
}

ISOELECTRIC = {  # ZIMJ680104 — Zimmerman et al. (1968)
    "A": 6.00,  "R": 10.76, "N": 5.41, "D": 2.77, "C": 5.07,
    "Q": 5.65,  "E": 3.22,  "G": 5.97, "H": 7.59, "I": 6.02,
    "L": 5.98,  "K": 9.74,  "M": 5.74, "F": 5.48, "P": 6.30,
    "S": 5.68,  "T": 5.60,  "W": 5.89, "Y": 5.66, "V": 5.96,
}

HELIX = {  # CHOP780101 — Chou & Fasman (1978)
    "A": 1.42, "R": 0.98, "N": 0.67, "D": 1.01, "C": 0.70,
    "Q": 1.11, "E": 1.51, "G": 0.57, "H": 1.00, "I": 1.08,
    "L": 1.21, "K": 1.16, "M": 1.45, "F": 1.13, "P": 0.57,
    "S": 0.77, "T": 0.83, "W": 1.08, "Y": 0.69, "V": 1.06,
}

DEFAULT_SCALES = {
    "hydropathy":  HYDROPATHY,
    "polarity":    POLARITY,
    "isoelectric": ISOELECTRIC,
    "helix":       HELIX,
}

# Component weights for the *inner* score (after the identity gate).
# Calibrated against a 25-protein safe panel: 0/13 fingerprints exceed
# the 0.60 threshold under this form (vs 2/13 under the previous additive
# form). Identity itself is used as a multiplicative gate, not a weight.
# Threshold 0.60 and these weights are provisional — recalibration against
# a 60-protein panel is accelerator-phase work.
DEFAULT_WEIGHTS = {
    "blosum":      0.30,
    "aaindex":     0.20,
    "fp_identity": 0.50,
}


# --------------------------------------------------------------------------
# Public class
# --------------------------------------------------------------------------


def identity_tolerance(native_aa: str, query_aa: str, blosum) -> float:
    """Smooth identity: 1.0 for exact match, partial for conservative.

    Y->F (BLOSUM ~3)  -> ~0.83
    E->D (BLOSUM ~2)  -> ~0.67
    R->K (BLOSUM ~2)  -> ~0.67
    Y->S (BLOSUM ~-1) -> 0.0
    Y->A (BLOSUM ~-2) -> 0.0

    Conservative substitutions get partial credit; non-conservative
    substitutions get zero. Consistent with the motif gate's
    CONSERVATIVE_SUBSTITUTIONS table.
    """
    if native_aa == query_aa:
        return 1.0
    try:
        b = blosum[native_aa][query_aa]
    except (KeyError, IndexError):
        return 0.0
    # Map BLOSUM score [-2, +4] to [0, 1]
    return max(0.0, min(1.0, (b + 2.0) / 6.0))


class ContextFreeFingerprint:
    """Context-free fingerprint for one toxin's functional region.

    Parameters
    ----------
    region_sequence : str
        The toxin's functional region (e.g. P11439 430-638).
    catalytic_positions : list[int]
        0-based indices into region_sequence where catalytic residues sit.

        IMPORTANT — choose these deliberately. UniProt's "Binding site"
        annotation includes many second-shell and cofactor contacts that
        tolerate substitution (e.g. P11439 has 12 NAD+ contacts but only
        4 literature-curated catalytic residues). Using the wider set
        asks the wrong question: does the variant preserve positions
        EvoDiff was never told to protect. Use the catalytic set from
        TOXIN_CONFIG or an equivalent literature source.
    window : int, default 6
        Half-width of the fingerprint window. Matches the ESM fingerprint
        definition (±6 residues around each catalytic position).
    scales : dict[str, dict[str, float]], optional
        Physicochemical scales. Defaults to the four AAindex scales above.
    weights : dict[str, float], optional
        Component weights for the combined score.

    Attributes
    ----------
    fingerprint_indices : list[int]
        Sorted, deduplicated 0-based indices of the fingerprint window.
    catalytic_residues : list[str]
        Single-letter codes of the catalytic residues in the native region.
    """

    def __init__(
        self,
        region_sequence: str,
        catalytic_positions: Iterable[int],
        window: int = 6,
        scales: Optional[dict] = None,
        weights: Optional[dict] = None,
    ) -> None:
        self.region = region_sequence.upper()
        self.cat_positions = sorted(set(int(p) for p in catalytic_positions))
        self.window = int(window)
        self.scales = scales if scales is not None else DEFAULT_SCALES
        self.weights = weights if weights is not None else DEFAULT_WEIGHTS

        L = len(self.region)
        fp = set()
        for p in self.cat_positions:
            for off in range(-self.window, self.window + 1):
                q = p + off
                if 0 <= q < L:
                    fp.add(q)
        self.fingerprint_indices = sorted(fp)
        self.catalytic_residues = [self.region[p] for p in self.cat_positions]

        self._blosum = substitution_matrices.load("BLOSUM62")

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def score(self, query_region: str) -> dict:
        """Score a query region against the fingerprint.

        Parameters
        ----------
        query_region : str
            Same length as the native region (no alignment step).

        Returns
        -------
        dict
            Component scores and the combined score.
        """
        q = query_region.upper()
        if len(q) != len(self.region):
            raise ValueError(
                f"query length {len(q)} != region length {len(self.region)}"
            )

        identity      = self._identity(q)
        blosum_norm   = self._blosum_norm(q)
        aaindex_corr  = self._aaindex_corr(q)
        fp_identity   = self._fp_identity(q)

        # Gated form. Identity at the catalytic positions is a *precondition*,
        # not a score. A chance match of two catalytic residues in an unrelated
        # protein should not add to the combined score — it should gate the
        # secondary signals. This form was calibrated on 25 safe proteins:
        # additive  -> safe_max 0.61 (2/13 fingerprints fail at threshold 0.60)
        # gated     -> safe_max 0.41 (0/13 fail)
        # On P11439 true positives: native 0.9659, aggressive 0.0000.
        inner = (
            self.weights["fp_identity"] * fp_identity
            + self.weights["blosum"]    * blosum_norm
            + self.weights["aaindex"]   * aaindex_corr
        )
        combined = identity * inner

        return {
            "identity":    round(float(identity), 4),
            "blosum":      round(float(blosum_norm), 4),
            "aaindex":     round(float(aaindex_corr), 4),
            "fp_identity": round(float(fp_identity), 4),
            "combined":    round(float(combined), 4),
        }

    # ------------------------------------------------------------------
    # Component scores
    # ------------------------------------------------------------------
    def _identity(self, q: str) -> float:
        vals = [
            identity_tolerance(self.catalytic_residues[i], q[p], self._blosum)
            for i, p in enumerate(self.cat_positions)
        ]
        return float(np.mean(vals))

    def _blosum_norm(self, q: str) -> float:
        # BLOSUM62 identity scores are typically 4-8; conservative
        # substitutions 1-3; radical substitutions negative. We map
        # [-3, +8] -> [0, 1], capped.
        scores = [
            self._blosum[self.catalytic_residues[i]][q[p]]
            for i, p in enumerate(self.cat_positions)
        ]
        mean = float(np.mean(scores))
        return max(0.0, min(1.0, (mean + 3.0) / 11.0))

    def _aaindex_corr(self, q: str) -> float:
        corrs = []
        for scale in self.scales.values():
            n = np.array([scale[self.region[i]] for i in self.fingerprint_indices])
            v = np.array([scale[q[i]]          for i in self.fingerprint_indices])
            if n.std() > 1e-9 and v.std() > 1e-9:
                corrs.append(float(np.corrcoef(n, v)[0, 1]))
        if not corrs:
            return 0.0
        return max(0.0, float(np.mean(corrs)))

    def _fp_identity(self, q: str) -> float:
        hits = sum(1 for i in self.fingerprint_indices if q[i] == self.region[i])
        return hits / len(self.fingerprint_indices)



# --------------------------------------------------------------------------
# Motif patterns per mechanism class
# --------------------------------------------------------------------------
# A catalytic motif is a set of residue *types* at fixed relative
# positions. If the query doesn't reproduce them, it cannot be a
# functional variant — reject regardless of embedding similarity.
# Conservative substitutions are permitted (Y<->F, E<->D, R<->K, H<->K/R).

MOTIFS = {
    "N-glycosidase":     {"required": ["Y", "Y", "E", "R"]},
    "Zinc protease":     {"required": ["H", "E", "H"]},
    "rRNA depurination": {"required": ["E", "R"]},
    "Adenylate cyclase": {"required": ["H"]},
    "ADP-ribosylation":  {"required": None},   # variable active-site chemistry
}

CONSERVATIVE_SUBSTITUTIONS = {
    "Y": {"F", "W"}, "F": {"Y", "W"},
    "E": {"D"},      "D": {"E"},
    "R": {"K"},      "K": {"R"},
    "H": {"K", "R"},
    "N": {"Q"},      "Q": {"N"},
    "S": {"T"},      "T": {"S"},
    "A": {"G", "S"}, "G": {"A", "S"},
    "V": {"I", "L"}, "I": {"V", "L"}, "L": {"I", "V"},
}

# Threshold for the per-site max verdict.
PER_SITE_THRESHOLD = 0.60


# --------------------------------------------------------------------------
# Public helper — motif check
# --------------------------------------------------------------------------
def motif_passes(query_region: str, native_region: str,
                 cat_positions: list[int], mechanism: str
                 ) -> tuple[bool, str]:
    """Check whether the query preserves the mechanism's catalytic motif.

    Returns (passes, reason). ``reason`` is a short diagnostic string.
    """
    required = MOTIFS.get(mechanism, {}).get("required", "unknown")
    if required is None:
        return True, "no-motif-class"
    if required == "unknown":
        return True, f"unknown-mechanism:{mechanism}"

    if len(query_region) != len(native_region):
        return False, f"length-mismatch({len(query_region)}!={len(native_region)})"

    for idx in cat_positions:
        n_aa = native_region[idx]
        q_aa = query_region[idx]
        if q_aa == n_aa:
            continue
        if q_aa in CONSERVATIVE_SUBSTITUTIONS.get(n_aa, set()):
            continue
        return False, f"site-{idx}:{n_aa}->{q_aa}"
    return True, "motif-preserved"


# --------------------------------------------------------------------------
# Patch the ContextFreeFingerprint class
# --------------------------------------------------------------------------
def _per_site_max(self, query_region: str) -> float:
    """Score each catalytic site's local window independently.

    Returns the max across sites. Unlike the combined score, this does
    NOT require every site to be preserved — a variant that keeps any
    one functional site fires.
    """
    q = query_region.upper()
    if len(q) != len(self.region):
        return 0.0

    scores = []
    for cat_idx in self.cat_positions:
        lo = max(0, cat_idx - self.window)
        hi = min(len(self.region), cat_idx + self.window + 1)
        idx = list(range(lo, hi))

        native_res = self.region[cat_idx]
        query_res = q[cat_idx]
        identity = identity_tolerance(native_res, query_res, self._blosum)

        blosum = self._blosum[native_res][query_res]
        bl_norm = max(0.0, min(1.0, (blosum + 3.0) / 11.0))

        corrs = []
        for scale in self.scales.values():
            n = np.array([scale[self.region[i]] for i in idx])
            v = np.array([scale[q[i]] for i in idx])
            if n.std() > 1e-9 and v.std() > 1e-9:
                corrs.append(float(np.corrcoef(n, v)[0, 1]))
        aai = max(0.0, float(np.mean(corrs))) if corrs else 0.0

        fp_ident = sum(1 for i in idx if q[i] == self.region[i]) / len(idx)
        scores.append(identity * (0.50 * fp_ident + 0.30 * bl_norm + 0.20 * aai))

    return max(scores) if scores else 0.0


def _motif_check(self, query_region: str, mechanism: str
                 ) -> tuple[bool, str]:
    """Convenience wrapper — calls motif_passes with self's region."""
    return motif_passes(query_region, self.region, self.cat_positions, mechanism)


def _verdict(self, query_region: str, mechanism: str) -> dict:
    """Combined verdict: motif gate AND per-site max.

    Returns a dict with the individual signals and the final decision.
    The verdict is:
        fires = motif_pass AND per_site_max >= PER_SITE_THRESHOLD
    """
    motif_ok, reason = self._motif_check(query_region, mechanism)
    psm = self._per_site_max(query_region) if motif_ok else 0.0
    fires = bool(motif_ok and psm >= PER_SITE_THRESHOLD)
    return {
        "fires": fires,
        "motif_pass": bool(motif_ok),
        "motif_reason": reason,
        "per_site_max": round(float(psm), 4),
        "combined_score": round(float(self.score(query_region)["combined"]), 4),
        "threshold": PER_SITE_THRESHOLD,
    }


# --------------------------------------------------------------------------
# Attach to the class
# --------------------------------------------------------------------------
ContextFreeFingerprint._per_site_max = _per_site_max
ContextFreeFingerprint._motif_check = _motif_check
ContextFreeFingerprint.per_site_max = _per_site_max
ContextFreeFingerprint.motif_check = _motif_check
ContextFreeFingerprint.verdict = _verdict

__all__ = [
    "ContextFreeFingerprint",
    "DEFAULT_SCALES",
    "DEFAULT_WEIGHTS",
    "MOTIFS",
    "CONSERVATIVE_SUBSTITUTIONS",
    "PER_SITE_THRESHOLD",
    "motif_passes",
    "identity_tolerance",
]
