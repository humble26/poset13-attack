"""E24: what does the gap-weighting do to the p-matrix?

For fail cases with multi-point box (z non-extreme):
 1. p_P(k,l) vs p_Q(k,l) per cell: shift magnitude and sign vs the theta-cuts
 2. central-symmetry deficit: d(k,l) = p_P(k,l) + p_P(s+1-k, t+1-l) - 1:
    range and sign pattern
 3. which cells are balanced in P but not Q (and vice versa); where do the
    S1 survivors sit (deep in D/U? near window?)
"""
import random
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e11_corner.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e22b_s3.py").resolve(), encoding="utf-8").read().split("random.seed(314)")[0].replace("def pmatrix", "def pmatrixOLD"))
exec(open((Path(__file__).parent / "e22b_s3.py").resolve(), encoding="utf-8").read().split("def pmatrix(less")[1].split("random.seed")[0].join(["def pmatrix(less", ""]))

third = Fraction(1, 3); twothird = Fraction(2, 3)
random.seed(2024)


def pQ_cell(Ql, m, Q, gx, gy, extsQ):
    return Fraction(sum(1 for e in extsQ if e.index(Q.index(gx)) < e.index(Q.index(gy))), len(extsQ))


def main():
    n = 6
    npairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    stats = dict(cases=0, dmin=Fraction(0), dmax=Fraction(0),
                 shift_pos=0, shift_neg=0, shift_zero=0)
    dvals = []
    surv_deep = surv_win = 0
    examples = []
    tried = 0
    while tried < 20000 and stats['cases'] < 400:
        tried += 1
        mask = random.randrange(1 << len(npairs))
        rel = [npairs[k] for k in range(len(npairs)) if (mask >> k) & 1]
        L = closure(n, rel)
        if not is_poset(L, n): continue
        if antichain_width(L, n) < 2: continue
        z = random.randrange(n)
        Qz = [v for v in range(n) if v != z]
        Ql, _ = induced(L, n, Qz)
        if antichain_width(Ql, n - 1) > 2: continue
        r = analyze(L, n, z)
        if r['zpair']: continue
        for part in chain_partitions(Ql, n - 1):
            cc = corner_check(L, n, z, part)
            if cc is None: continue
            an = cc['an']
            if len(an['mu']) < 2: continue   # multi-cut box only
            stats['cases'] += 1
            A, B = part
            s, t = an['s'], an['t']
            Q = [v for v in range(n) if v != z]
            Ql2, _ = induced(L, n, Q)
            extsQ = []
            def rec(placed, seq):
                if len(seq) == n - 1: extsQ.append(list(seq)); return
                for v in range(n - 1):
                    if (placed >> v) & 1: continue
                    if all((placed >> u) & 1 for u in range(n - 1) if Ql2[u][v]):
                        rec(placed | 1 << v, seq + [v])
            rec(0, [])
            pm = pmatrix(L, n, z, part, an)
            istar = cc['istar']; jstar = cc['jstar']
            for k in range(1, s + 1):
                for l in range(1, t + 1):
                    rel_, pP = pm[(k, l)]
                    if rel_ != 'inc': continue
                    gx, gy = Q[A[k - 1]], Q[B[l - 1]]
                    pQv = Fraction(sum(1 for e in extsQ
                                       if e.index(Q.index(gx)) < e.index(Q.index(gy))), len(extsQ))
                    # central-symmetry deficit
                    k2, l2 = s + 1 - k, t + 1 - l
                    if (k2, l2) in pm and pm[(k2, l2)][0] == 'inc':
                        pP2 = pm[(k2, l2)][1]
                        d = pP + pP2 - 1
                        dvals.append(d)
                        if d < stats['dmin']: stats['dmin'] = d
                        if d > stats['dmax']: stats['dmax'] = d
                    # shift sign relative to cuts
                    sh = pP - pQv
                    if sh > 0: stats['shift_pos'] += 1
                    elif sh < 0: stats['shift_neg'] += 1
                    else: stats['shift_zero'] += 1
            if len(examples) < 3:
                examples.append((rel, z, part, istar, jstar,
                                 {kk: str(v) for kk, v in sorted(an['mu'].items())}))
    import collections
    hist = collections.Counter(dvals)
    print(stats)
    print("d (central-symmetry deficit) histogram:", dict(sorted(hist.items())[:12]))
    for (rel_, z_, part_, ist_, jst_, mu_) in examples:
        print("example rel=%s z=%d A=%s B=%s i*=%s j*=%s mu=%s" % (rel_, z_, part_[0], part_[1], ist_, jst_, mu_))


if __name__ == '__main__':
    main()
