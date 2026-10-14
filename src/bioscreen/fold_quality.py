"""Fold-confidence gate for geometry-based screening.

A geometric spec — catalytic triad, Zn coordination, 3Di token — is only
meaningful if the input structure is a real fold. ESMFold will produce a
PDB for any sequence, including random coil. pLDDT is the per-residue
confidence that the predicted backbone is correct.

Rule: before any geometry check, verify mean pLDDT >= 50 (0-100 scale).
Below 50 is effectively disordered (Jumper et al. 2021).

SCALE NOTE
----------
The transformers ESMFold implementation writes pLDDT in the 0-1 range,
while the original ESMFold paper and most literature report 0-100. This
module auto-detects the scale from the values themselves (max <= 1.5
implies 0-1 scale) and normalizes to 0-100 before comparison. Threshold
values are always in 0-100.
"""

import numpy as np
from Bio.PDB import PDBParser

MEAN_PLDDT_MIN = 50.0            # 0-100 scale
ACTIVE_SITE_PLDDT_MIN = 50.0     # 0-100 scale

_parser = PDBParser(QUIET=True)


def _to_100(vals):
    """Normalize pLDDT values to the 0-100 scale, in place-safe form."""
    if len(vals) == 0:
        return vals
    return vals * 100.0 if vals.max() <= 1.5 else vals


def plddt_from_pdb(pdb_path, as_100=True):
    """Extract per-residue pLDDT from an ESMFold PDB (CA B-factor column).

    Returns values on the 0-100 scale by default. Pass as_100=False to
    get raw B-factor values as stored.
    """
    st = _parser.get_structure("x", str(pdb_path))
    vals = []
    for residue in st[0].get_residues():
        if "CA" in residue:
            vals.append(residue["CA"].get_bfactor())
    arr = np.array(vals, dtype=float) if vals else np.array([])
    return _to_100(arr) if as_100 else arr


def fold_is_interpretable(pdb_path, min_mean=MEAN_PLDDT_MIN,
                          min_activesite=ACTIVE_SITE_PLDDT_MIN,
                          active_site_slice=None):
    """Return (ok: bool, reason: str, stats: dict).

    All pLDDT comparisons are on the normalized 0-100 scale, regardless
    of what scale the PDB stored.
    """
    plddt = plddt_from_pdb(pdb_path)   # already 0-100
    if len(plddt) == 0:
        return False, "no residues with CA in PDB", {}

    stats = {
        "n_residues": int(len(plddt)),
        "mean_plddt": float(plddt.mean()),
        "median_plddt": float(np.median(plddt)),
        "max_plddt": float(plddt.max()),
        "frac_ge70": float((plddt >= 70).mean()),
    }

    if plddt.mean() < min_mean:
        return False, f"mean pLDDT {plddt.mean():.1f} < {min_mean}", stats

    if active_site_slice is not None:
        s, e = active_site_slice
        as_mean = float(plddt[s:e].mean())
        stats["active_site_plddt"] = as_mean
        if as_mean < min_activesite:
            return False, (f"active-site pLDDT {as_mean:.1f} "
                           f"< {min_activesite}"), stats

    return True, "structure interpretable", stats


def geometry_check_guarded(pdb_path, geometry_fn, *args, **kwargs):
    """Run geometry_fn only if the fold passes the gate.

    If the fold is not interpretable, returns evaluable=False and
    passes=None — NOT passes=False. We decline to evaluate rather than
    declaring a negative.
    """
    ok, reason, stats = fold_is_interpretable(pdb_path)
    if not ok:
        return {
            "evaluable": False,
            "gate_reason": reason,
            "gate_stats": stats,
            "passes": None,
            "fired": False,
        }
    result = geometry_fn(pdb_path, *args, **kwargs)
    result["evaluable"] = True
    result["gate_stats"] = stats
    return result
