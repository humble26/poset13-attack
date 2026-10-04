"""E30: weighted vertex mass WV(i*,j*) in fail cases.
WV = sum over M in ext(Q) whose prefix of length i*+j* equals I(i*,j*) of gap(M)/e(P).
Verify the four edge masses alpha1, alpha2, beta1, beta2 are each <= WV
(generic corner, all four cells incomparable), and record WV's distribution.
"""
import random
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e11_corner.py").resolve(), encoding="utf-8").read().split("def main():")[0])

third = Fraction(1, 3); twothird = Fraction(2, 3)
random.seed(1234)


def gap_of(M, Dpos, Upos, m):
    # D, U as Q-label sets; t = max position of D (1-based), nu = min position of U
    t = 0
    for v in Dpos:
        p = M.index(v) + 1
        if p > t: t = p
    nu = m + 1
    for v in Upos:
        p = M.index(v) + 1
        if p < nu: nu = p
    return (nu - t) if nu > t else 0


def main():
    ok_dom = bad_dom = 0
    wv_hist = {}
    n_case = 0
    for n in (5, 6):
        npairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        if n == 5:
            it = [(mask, random.randrange(n)) for mask in range(1 << len(npairs))
                  for _ in (0,) for z in [None] for mask2 in [0]]
            it = [(mask, z) for mask in range(1 << len(npairs)) for z in range(n)]
        else:
            it = []
            while len(it) < 30000:
                it.append((random.randrange(1 << len(npairs)), random.randrange(n)))
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
                ist, jst = cc['istar'], cc['jstar']
                s_, t_ = an['s'], an['t']
                if ist == 0 or jst == 0 or ist + 1 > s_ or jst + 1 > t_: continue
                Q = [v for v in range(n) if v != z]
                A, B = part
                gx = {}
                okgen = True
                for (nm, kk, ll) in (('m1', ist, jst + 1), ('m2', ist + 1, jst),
                                     ('c1', ist, jst), ('c2', ist + 1, jst + 1)):
                    x, y = Q[A[kk - 1]], Q[B[ll - 1]]
                    if L[x][y] or L[y][x]: okgen = False; break
                    gx[nm] = (x, y)
                if not okgen: continue
                n_case += 1
                m1 = pval(L, n, *gx['m1']); m2 = pval(L, n, *gx['m2'])
                c1 = pval(L, n, *gx['c1']); c2 = pval(L, n, *gx['c2'])
                a1 = c1 - m2; a2 = m1 - c2; b1 = m1 - c1; b2 = c2 - m2
                # WV: enumerate ext(Q)
                Ql2, _ = induced(L, n, Qz)
                m = n - 1
                extsQ = []
                def rec(placed, seq):
                    if len(seq) == m: extsQ.append(list(seq)); return
                    for v in range(m):
                        if (placed >> v) & 1: continue
                        if all((placed >> u) & 1 for u in range(m) if Ql2[u][v]):
                            rec(placed | 1 << v, seq + [v])
                rec(0, [])
                Dq = {v for v in range(m) if L[Q[v]][z]}
                Uq = {v for v in range(m) if L[z][Q[v]]}
                eP = count_ext(L, n)
                wv = Fraction(0)
                Iset = set(A[:ist]) | set(B[:jst])
                for M in extsQ:
                    if set(M[:ist + jst]) == Iset:
                        wv += gap_of(M, Dq, Uq, m)
                wv = Fraction(wv, eP)
                edges = {'a1': a1, 'a2': a2, 'b1': b1, 'b2': b2}
                if all(v_ <= wv for v_ in edges.values()): ok_dom += 1
                else:
                    bad_dom += 1
                    if bad_dom <= 3:
                        print("DOM-FAIL", n, z, {k: str(v_) for k, v_ in edges.items()}, "WV=", str(wv))
                key = wv
                wv_hist[key] = wv_hist.get(key, 0) + 1
    print("generic corner cases:", n_case, "| edge<=WV: ok=%d bad=%d" % (ok_dom, bad_dom))
    big = sum(v for k, v in wv_hist.items() if k > third)
    tot = sum(wv_hist.values())
    print("WV > 1/3 in %d / %d cases" % (big, tot))
    print("WV top values:", sorted(((k, v) for k, v in wv_hist.items()), reverse=True)[:6])


if __name__ == '__main__':
    main()
