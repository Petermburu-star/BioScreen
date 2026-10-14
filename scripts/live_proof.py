"""
BioScreen LIVE PROOF v2 — reviewer-driven end-to-end demonstration.
Reviewer picks toxin (+ seed for global RNGs). Variants generated LIVE,
scored LIVE; controls must pass or the protocol invalidates itself.

Usage: python scripts/live_proof.py --accession P02879 --seed 31415 --n 1
"""
import argparse, hashlib, json, random, sys, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

def sha(s): return hashlib.sha256(s.encode()).hexdigest()[:16]
def pid(a, b):
    n = min(len(a), len(b))
    return sum(x == y for x, y in zip(a[:n], b[:n])) / n

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--accession", default="P02879")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--n", type=int, default=1)
    args = ap.parse_args()
    t_start = time.time()
    print("=" * 74)
    print(f"BIOSCREEN LIVE PROOF v2 | toxin {args.accession} | seed {args.seed}")
    print("  variants generated live; controls enforce themselves")
    print("=" * 74)

    db = json.loads((ROOT / "data/product_db/bioscreen_db.json").read_text())
    e = next(x for x in db["entries"] if x["accession"] == args.accession)
    region = e["region_sequence"]
    cat_rel = e["catalytic_region_rel"]
    print(f"\n[1] TARGET: {e['name']} ({e['accession']}) class={e['class']}")
    print(f"    catalytic (1-based precursor): {e['catalytic_1based']}")
    print(f"    residues at catalytic positions: {e['residues']}")
    print(f"    verification: {e['verification']}")

    # HARD GATE — residues at the positions we are about to protect must
    # match the DB's verified residue string. No off-by-one survives this.
    at = "".join(region[c] for c in cat_rel)
    assert at == e["residues"], f"frame mismatch: '{at}' != '{e['residues']}'"
    print(f"    frame gate: residues at region-rel {cat_rel} == '{at}'  OK")

    print(f"\n[2] LIVE GENERATION (EvoDiff; seed {args.seed} applied to global RNGs)")
    try:
        import torch
        torch.manual_seed(args.seed)
    except Exception:
        pass
    random.seed(args.seed)
    from bioscreen.generate import generate_variants_custom
    outs = {}
    for mode in ("preserved", "aggressive"):
        t0 = time.time()
        outs[mode] = generate_variants_custom(
            sequence=region,              # region as the sequence
            region_start=0, region_end=None,
            catalytic_residues=cat_rel,   # 0-based region-relative (line-49 convention)
            auto_catalytic=False,
            mode=mode, n_variants=args.n,
            protection_window=8)          # protects +-8, superset of CF +-6
        print(f"    {mode}: {args.n} variant(s) in {time.time()-t0:.0f}s")

    def variants_of(out):
        if isinstance(out, dict):
            for k in ("variants", "results", "data"):
                if isinstance(out.get(k), list):
                    return out[k]
            return []
        return out if isinstance(out, list) else []
    seqs = {m: [v["sequence"] if isinstance(v, dict) else str(v)
                for v in variants_of(o)] for m, o in outs.items()}

    stored = {x["region_sequence"] for x in db["entries"]}
    print("\n[3] FRESHNESS (new sequences, exist in no stored dataset)")
    for m in seqs:
        for i, s in enumerate(seqs[m]):
            print(f"    {m}#{i}: sha={sha(s)} len={len(s)} "
                  f"ident-to-native={pid(s, region)*100:.1f}% in-db={s in stored}")

    from bioscreen.kernel import screen_sequence
    from bioscreen.contextfree import ContextFreeFingerprint as CFF
    fp = CFF(region, cat_rel)
    expect = {"preserved": "FLAGGED", "aggressive": "CLEAR"}
    rows, ok = [], True
    print("\n[4] SCORING (CF kernel, full decomposition)")
    for m in seqs:
        for i, s in enumerate(seqs[m]):
            r = screen_sequence(s)
            dec = fp.score(s)
            good = r["status"] == expect[m]
            ok &= good
            rows.append(dict(tag=f"{m}#{i}", sha=sha(s), status=r["status"],
                             expect=expect[m], score=r["top_score"], dec=dec))
            print(f"    {m}#{i}: {r['status']:<8} CF={r['top_score']:.4f} "
                  f"[gate={dec['identity']:.2f} blosum={dec['blosum']:.2f} "
                  f"aaindex={dec['aaindex']:.2f} fp={dec['fp_identity']:.2f}] "
                  f"-> {'PASS' if good else 'DEMO-INVALID'}")

    print("\n  CONTROL ASSERTION: " +
          ("PROTOCOL VALID" if ok else "DEMO INVALIDATED — investigate"))
    print(f"  total wall time: {time.time()-t_start:.0f}s")
    outp = ROOT / "results" / f"live_proof_{args.accession}_{args.seed}.json"
    outp.write_text(json.dumps(dict(accession=args.accession, seed=args.seed,
                                    generator="bioscreen.generate.generate_variants_custom",
                                    rows=rows, controls_valid=ok), indent=2),
                    encoding="utf-8")
    print(f"  saved: {outp.relative_to(ROOT)}")

if __name__ == "__main__":
    main()
