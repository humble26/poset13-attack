"""E12: sharp claim C1 test.
In a fail case with theta-cuts i*, j*:
  - (0,0) cell (a_{i*}, b_{j*}) both-lo: balanced when incomparable?
  - (1,1) cell (a_{i*+1}, b_{j*+1}) both-hi: balanced when incomparable?
  - (0,1) mixed lo-hi: p > 2/3?   (1,0) mixed hi-lo: p < 1/3?
"""
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main()")[0])

third = Fraction(1, 3); twothird = Fraction(2, 3)


def corner_check(less, n, z, part):
    A, B = part
    an = anatomy(less, n, z, part)
    if an is None: return None
    s, t = an['s'], an['t']
    Q = [v for v in range(n) if v != z]
    MAv, MBv = an['MAv'], an['MBv']
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
    cells = {}
    for di in (0, 1):
        for dj in (0, 1):
            ka, lb = istar + di, jstar + dj
            if ka == 0 or lb == 0 or ka > s or lb > t: continue
            gx, gy = Q[A[ka-1]], Q[B[lb-1]]
            if less[gx][gy] or less[gy][gx]:
                cells[(di, dj)] = ('comp', None)
                continue
            p = pval(less, n, gx, gy)
            cells[(di, dj)] = ('inc', p)
    return dict(an=an, istar=istar, jstar=jstar, cells=cells)
stat = dict(c00_inc=0, c00_bal=0, c11_inc=0, c11_bal=0,
            c01_inc=0, c01_hi=0, c10_inc=0, c10_lo=0)
bad = []
nmax = 5
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
                cells = cc['cells']
                for key, incst, balst in [((0, 0), 'c00_inc', 'c00_bal'),
                                          ((1, 1), 'c11_inc', 'c11_bal'),
                                          ((0, 1), 'c01_inc', 'c01_hi'),
                                          ((1, 0), 'c10_inc', 'c10_lo')]:
                    c = cells.get(key)
                    if c is None or c[0] != 'inc': continue
                    p = c[1]
                    stat[incst] += 1
                    if key in ((0, 0), (1, 1)):
                        if third <= p <= twothird: stat[balst] += 1
                        else: bad.append((n, z, rel, part, key, str(p)))
                    elif key == (0, 1):
                        if p > twothird: stat[balst] += 1
                        else: bad.append((n, z, rel, part, key, str(p)))
                    else:
                        if p < third: stat[balst] += 1
                        else: bad.append((n, z, rel, part, key, str(p)))
print(stat)
print("violations:", len(bad))
for b in bad[:8]:
    print(b)
