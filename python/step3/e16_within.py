"""E16: hypothesis H — corner pairs are within-I balanced for EVERY box cut.
For fail cases: for corner pairs (a_{i*},b_{j*}) and (a_{i*+1},b_{j*+1}),
compute p_I(x,y) for every cut I in the box containing both; check all in [1/3,2/3].
Also record p_P and the crossing masses to see the full balance ledger.
"""
import sys
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e11_corner.py").resolve(), encoding="utf-8").read().split("def main():")[0])

third = Fraction(1, 3); twothird = Fraction(2, 3)


def within_ledger(less, n, z, part, an, gx, gy):
    """for pair (gx,gy): ledger over cuts; return dict with within-range flag,
    crossing masses, pP."""
    A, B = part
    Q = [v for v in range(n) if v != z]
    xl, yl = Q.index(gx), Q.index(gy)
    k = A.index(xl) + 1 if xl in A else None
    l = B.index(yl) + 1 if yl in B else None
    mu = an['mu']
    cross = anti = Fraction(0)
    wvals = []
    for (i, j), w in mu.items():
        x_in = k <= i; y_in = l <= j
        if x_in and not y_in: cross += w
        elif y_in and not x_in: anti += w
        elif x_in and y_in:
            xf, yf, tot = inter_count(less, [Q[v] for v in A[:i]], [Q[v] for v in B[:j]], gx, gy)
            wvals.append((w, Fraction(xf, tot)))
        else:
            xf, yf, tot = inter_count(less, [Q[v] for v in A[i:]], [Q[v] for v in B[j:]], gx, gy)
            wvals.append((w, Fraction(xf, tot)))
    inrange = all(third <= v <= twothird for (_, v) in wvals)
    pP = pval(less, n, gx, gy)
    ledger = cross + sum(w * v for (w, v) in wvals)
    assert ledger == pP, ("ledger mismatch", ledger, pP)
    return dict(inrange=inrange, cross=cross, anti=anti, pP=pP,
                wmin=min((v for _, v in wvals), default=None),
                wmax=max((v for _, v in wvals), default=None))


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    agg = dict(pairs=0, H_ok=0, H_bad=0, ledger_ok=0, ledger_bad=0)
    badex = []
    for n in range(4, nmax + 1):
        npairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        for mask in range(1 << len(npairs)):
            rel = [npairs[k] for k in range(len(npairs)) if (mask >> k) & 1]
            L = closure(n, rel)
            if not is_poset(L, n): continue
            if antichain_width(L, n) < 2: continue
            for z in range(n):
                Qz = [v for v in range(n) if v != z]
                Ql, _ = induced(L, n, Qz)
                if antichain_width(Ql, n - 1) > 2: continue
                r = analyze(L, n, z)
                if r['zpair']: continue
                for part in chain_partitions(Ql, n - 1):
                    cc = corner_check(L, n, z, part)
                    if cc is None: continue
                    an = cc['an']; istar = cc['istar']; jstar = cc['jstar']
                    A, B = part; s, t = an['s'], an['t']
                    Q = [v for v in range(n) if v != z]
                    for (ka, lb) in [(istar, jstar), (istar + 1, jstar + 1)]:
                        if ka == 0 or lb == 0 or ka > s or lb > t: continue
                        gx, gy = Q[A[ka - 1]], Q[B[lb - 1]]
                        if L[gx][gy] or L[gy][gx]: continue
                        agg['pairs'] += 1
                        led = within_ledger(L, n, z, part, an, gx, gy)
                        if led['inrange']: agg['H_ok'] += 1
                        else:
                            agg['H_bad'] += 1
                            if len(badex) < 8:
                                badex.append((n, z, rel, (ka, lb), led, an, part))
                        # ledger check: does cross/anti + within-range imply pP in [1/3,2/3]?
                        lo = led['wmin'] * (1 - led['cross'] - led['anti']) if led['wmin'] is not None else Fraction(0)
                        hi = led['wmax'] * (1 - led['cross'] - led['anti']) + led['cross'] if led['wmax'] is not None else led['cross']
                        if led['wmin'] is None:
                            lo, hi = led['cross'], led['cross']
                        ok = (third <= led['pP'] <= twothird)
                        implies = (lo >= third) and (hi <= twothird)
                        if ok == implies or (ok and not implies):
                            agg['ledger_ok'] += 1
                        else:
                            agg['ledger_bad'] += 1
    print(agg)
    for (n, z, rel, cell, led, an, part) in badex:
        print("H-BAD n=%d z=%d cell=%s pairP=%s cross=%s anti=%s within=[%s,%s]" % (
            n, z, cell, led['pP'], led['cross'], led['anti'], led['wmin'], led['wmax']))
        print("   rel=%s part=%s" % (rel, part))
        print("   mu=%s" % {k: str(v) for k, v in sorted(an['mu'].items())})


if __name__ == '__main__':
    main()
