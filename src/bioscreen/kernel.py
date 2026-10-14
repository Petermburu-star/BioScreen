"""
BioScreen product kernel v1.0.3. Fully offline.
Detection:   top combined score vs threshold (DB meta, 0.60).
Attribution: candidates need catalytic identity >= 0.6 AND window identity
             >= 0.5. Ties in window identity (within 0.02 in combined score)
             are reported as a tied family LIST — e.g. Botulinum A / Botulinum E
             share an identical HExxH catalytic window, so no window-level
             detector can separate them; we say so instead of guessing.
"""
import json, re, time
from pathlib import Path
try:
    from .contextfree import ContextFreeFingerprint
except ImportError:
    from bioscreen.contextfree import ContextFreeFingerprint

_REPO = Path(__file__).resolve().parents[2]
_DB_PATHS = [_REPO / "data" / "product_db" / "bioscreen_db.json",
             Path(__file__).resolve().parent / "bioscreen_db.json"]
_DB = None

def db_path():
    for p in _DB_PATHS:
        if p.exists():
            return p
    raise FileNotFoundError("bioscreen_db.json not found - run the freeze cell")

def load_db(path=None, force=False):
    global _DB
    if _DB is not None and not force and path is None:
        return _DB
    raw = json.loads((Path(path) if path else db_path()).read_text(encoding="utf-8"))
    entries = [{**e, "cf": ContextFreeFingerprint(e["region_sequence"],
                                                  e["catalytic_region_rel"])}
               for e in raw["entries"]]
    _DB = {"meta": {k: v for k, v in raw.items() if k != "entries"},
           "entries": entries}
    return _DB

def screen_sequence(seq, db=None):
    """Raw sequence in -> biosecurity verdict out."""
    if db is None:
        db = load_db()
    t0 = time.perf_counter()
    s = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", str(seq).upper())
    if len(s) < 50:
        raise ValueError(f"sequence too short after cleaning: {len(s)} aa (min 50)")
    hits = []
    for e in db["entries"]:
        try:
            sc = e["cf"].score(s)
        except Exception:
            continue
        hits.append((float(sc["combined"]), float(sc["identity"]),
                     float(sc["fp_identity"]), e))
    if not hits:
        raise RuntimeError("no fingerprint produced a score")
    hits.sort(key=lambda h: -h[0])
    thr = float(db["meta"]["threshold"])
    top_sc = hits[0][0]
    cands = [h for h in hits if h[1] >= 0.6 and h[2] >= 0.5]
    attrib = None
    if cands:
        max_w = max(h[2] for h in cands)
        tied = [h for h in cands if max_w - h[2] < 1e-9]
        best_s = max(h[0] for h in tied)
        tied = [h for h in tied if best_s - h[0] <= 0.02]
        tied.sort(key=lambda h: h[3]["name"])
        _, ident, fpid, e = tied[0]
        attrib = {"accession": e["accession"], "name": e["name"],
                  "names": [h[3]["name"] for h in tied],
                  "class": e["class"], "mechanism": e["mechanism"],
                  "verification": e["verification"],
                  "catalytic_identity": round(ident, 2),
                  "window_identity": round(fpid, 2),
                  "score": round(tied[0][0], 4)}
    return {"status": "FLAGGED" if top_sc >= thr else "CLEAR",
            "top_score": round(top_sc, 4), "threshold": thr,
            "best_match": attrib,
            "top5": [{"accession": h[3]["accession"], "name": h[3]["name"],
                       "class": h[3]["class"], "score": round(h[0], 4)}
                      for h in hits[:5]],
            "n_scanned": len(db["entries"]), "query_len": len(s),
            "seconds": round(time.perf_counter() - t0, 2),
            "db_version": db["meta"]["db_version"]}
