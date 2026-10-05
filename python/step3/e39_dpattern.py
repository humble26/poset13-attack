"""E39: anatomy of the balanced D-pair (for proving Lemma 9.D).
For fail + multi-cut + DxD-incomparable-exists cases:
 - how many balanced D-pairs per case; their (k,l) positions;
 - is the balanced pair's within-D value p_D also balanced?
 - is its within-Imax (max-mu cut) value also balanced?
 - is the balanced pair always on the D-matrix staircase boundary?
"""
import random, collections
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e11_corner.py").resolve(), encoding="utf-8").read().split("def main():")[0])

third = Fraction(1, 3); twothird = Fraction(2, 3)
random.seed(616)


def p_within(L, n, Ilist, x, y):
    Il, _ = induced(L, n, Ilist)
    extsI = []
    pos = {v: k for k, v in enumerate(Ilist)}
    def rec(placed, seq):
        if len(seq) == len(Ilist): extsI.append(list(seq)); return
        for v in Ilist:
            if (placed >> pos[v]) & 1: continue
            if all((placed >> pos[u]) & 1 for u in Ilist if L[u][v]):
                rec(placed | 1 << pos[v], seq + [v])
    rec(0, [])
    if not extsI: return None
    return Fraction(sum(1 for e in extsI if e.index(x) < e.index(y)), len(extsI))


def main():
    agg = dict(cases=0, balcnt=collections.Counter(),
               dD_bal=0, dD_tot=0, imax_bal=0, imax_tot=0,
               pos_hist=collections.Counter())
    fails = []
    for n in (5, 6):
        npairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        if n == 5:
            it = [(mask, z) for mask in range(1 << len(npairs)) for z in range(n)]
        else:
            it = [(random.randrange(1 << len(npairs)), random.randrange(n)) for _ in range(20000)]
        for mask, z in it:
            rel = [npairs[k] for k in range(len(npairs)) if (mask >> k) & 1]
            L = closure(n, rel)
            if not is_poset(L, n): continue
            if antichain_width(L, n) < 2: continue
            Qz = [v for v in range(n) if v != z]
            Ql, _ = induced(L, n, Qz)
            if antichain_width(Ql, n - 1) > 2: continue
            r = analyze(L, n, z)
            if r['zpair']: continue
            for part in chain_partitions(Ql, n - 1):
                cc = corner_check(L, n, z, part)
                if cc is None: continue
                an = cc['an']
                iD, jD = an['iD'], an['jD']
                if iD < 1 or jD < 1: continue
                if len(an['mu']) < 2: continue
                Q = [v for v in range(n) if v != z]
                A, B = part
                # DxD pairs with p-values
                pairs = []
                for k in range(1, iD + 1):
                    for l in range(1, jD + 1):
                        x, y = Q[A[k - 1]], Q[B[l - 1]]
                        if L[x][y] or L[y][x]: continue
                        p = pval(L, n, x, y)
                        pairs.append((k, l, x, y, p))
                if not pairs: continue
                agg['cases'] += 1
                bals = [(k, l, x, y, p) for (k, l, x, y, p) in pairs if third <= p <= twothird]
                agg['balcnt'][len(bals)] += 1
                if not bals:
                    if len(fails) < 3: fails.append((n, z, part, [(k, l, str(p)) for k, l, _, _, p in pairs]))
                    continue
                # pattern: positions of balanced pairs
                for (k, l, _, _, _) in bals:
                    agg['pos_hist'][(k, l)] += 1
                # check a balanced pair's within-D and within-Imax values
                k, l, x, y, p = bals[0]
                pD = p_within(L, n, sorted(set(A[:iD]) | set(B[:jD])), x, y)
                if third <= pD <= twothird: agg['dD_bal'] += 1
                agg['dD_tot'] += 1
                # max-mu cut
                Imax = max(an['mu'], key=an['mu'].get)
                Ilist = sorted(set(A[:Imax[0]]) | set(B[:Imax[1]]))
                pIm = p_within(L, n, Ilist, x, y)
                if pIm is not None and third <= pIm <= twothird: agg['imax_bal'] += 1
                if pIm is not None: agg['imax_tot'] += 1
    print(agg['cases'], "cases | balanced-pairs-per-case hist:", dict(sorted(agg['balcnt'].items())))
    print("balanced pair's within-D value balanced:", agg['dD_bal'], "/", agg['dD_tot'])
    print("balanced pair's within-Imax value balanced:", agg['imax_bal'], "/", agg['imax_tot'])
    print("balanced positions (top 8):", agg['pos_hist'].most_common(8))
    for f in fails: print("  NONE:", f)


if __name__ == '__main__':
    main()
