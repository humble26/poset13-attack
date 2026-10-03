"""E11: verify the sharp corner-pair claim.

theta(a_k) = 0 (a_k in D), M_A(k) (window), 1 (a_k in U); monotone along chain.
Fail => theta jumps the interval [1/3,2/3]: clean cuts i* (last theta<1/3), j*.
CLAIM: one of the corner pairs (a_{i*+di}, b_{j*+dj}), di,dj in {0,1}, that is
INCOMPARABLE, is BALANCED. Record which cells win per cut-config.
"""
import itertools, sys
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main()")[0])

third = Fraction(1, 3)
twothird = Fraction(2, 3)


def corner_check(less, n, z, part):
    A, B = part
    an = anatomy(less, n, z, part)
    if an is None: return None
    s, t = an['s'], an['t']
    Q = [v for v in range(n) if v != z]
    MAv, MBv = an['MAv'], an['MBv']
    # theta over full chains (1-based index; index s+1 = past-the-end theta=1)
    def thA(k):
        if k == 0: return Fraction(0)
        g = Q[A[k-1]]
        if less[g][z]: return Fraction(0)
        if less[z][g]: return Fraction(1)
        return MAv[k]
    def thB(l):
        if l == 0: return Fraction(0)
        g = Q[B[l-1]]
        if less[g][z]: return Fraction(0)
        if less[z][g]: return Fraction(1)
        return MBv[l]
    istar = max([k for k in range(0, s+1) if thA(k) < third], default=None)
    jstar = max([l for l in range(0, t+1) if thB(l) < third], default=None)
    assert istar is not None and jstar is not None
    # verify jump structure: theta > 2/3 above the cut
    for k in range(istar+1, s+1):
        assert thA(k) > twothird, ("thetaA not jumping", k, thA(k))
    for l in range(jstar+1, t+1):
        assert thB(l) > twothird, ("thetaB not jumping", l, thB(l))
    cells = {}
    rescue = None
    for di in (0, 1):
        for dj in (0, 1):
            ka, lb = istar + di, jstar + dj
            if ka == 0 or lb == 0 or ka > s or lb > t: continue
            gx, gy = Q[A[ka-1]], Q[B[lb-1]]
            if less[gx][gy] or less[gy][gx]:
                cells[(di, dj)] = ('comp', None)
                continue
            p = pval(less, n, gx, gy)
            bal = third <= p <= twothird
            cells[(di, dj)] = ('inc', p)
            if bal and rescue is None:
                rescue = (di, dj, gx, gy, p)
    return dict(an=an, istar=istar, jstar=jstar, cells=cells, rescue=rescue,
                sideA='alllo' if an['iHi'] is None else ('allhi' if an['iLo'] is None else 'jump'),
                sideB='alllo' if an['jHi'] is None else ('allhi' if an['jLo'] is None else 'jump'))


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    agg = {}
    bad = []
    ncase = 0
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
                    ncase += 1
                    key = (cc['sideA'], cc['sideB'])
                    a = agg.setdefault(key, dict(n=0, rescued=0, cell={}, mixedfirst=0, nocell=0, onlycomp=0))
                    a['n'] += 1
                    if cc['rescue'] is None:
                        incs = [kk for kk, vv in cc['cells'].items() if vv[0] == 'inc']
                        if not incs:
                            a['onlycomp'] += 1
                            if len(bad) < 10: bad.append((n, z, rel, part, cc))
                        else:
                            a['nocell'] += 1
                            if len(bad) < 10: bad.append((n, z, rel, part, cc))
                    else:
                        a['rescued'] += 1
                        di, dj = cc['rescue'][0], cc['rescue'][1]
                        # first winning cell in scan order; record same-side vs mixed
                        a['cell'][(di, dj)] = a['cell'].get((di, dj), 0) + 1
    print("corner claim over %d (case,partition)s:" % ncase)
    for key in sorted(agg):
        a = agg[key]
        print("  cfg %s: n=%d rescued=%d no-inc-cell=%d all-comp=%d cells=%s" % (
            key, a['n'], a['rescued'], a['nocell'], a['onlycomp'],
            {kk: vv for kk, vv in sorted(a['cell'].items())}))
    for (n, z, rel, part, cc) in bad:
        print("BAD n=%d z=%d rel=%s A=%s B=%s istar=%d jstar=%d cfg=%s" % (
            n, z, rel, part[0], part[1], cc['istar'], cc['jstar'], (cc['sideA'], cc['sideB'])))
        print("   cells:", cc['cells'])
        print("   MAv:", {k: str(v) for k, v in sorted(cc['an']['MAv'].items())})
        print("   MBv:", {k: str(v) for k, v in sorted(cc['an']['MBv'].items())})
        print("   mu:", {k: str(v) for k, v in sorted(cc['an']['mu'].items())})


if __name__ == '__main__':
    main()
