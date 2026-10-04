"""E32: (mu, WV) joint distribution at the corner in fail cases (n=6 sampled).
For generic-corner fail cases: measure
  mu = mu(i*,j*)  (corner cut mass, the balance-delivering quantity)
  WV = weighted vertex mass (what the configurations force > 1/3)
  the four edges, the configuration label under the hypothesis-sign pattern,
and check: edge <= WV, mu <= WV, and the WV/mu ratio distribution.
"""
import random, time
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e11_corner.py").resolve(), encoding="utf-8").read().split("def main():")[0])

third = Fraction(1, 3); twothird = Fraction(2, 3)
random.seed(4321)
t0 = time.time()


def gap_of(M, Dq, Uq, m):
    t = 0
    for v in Dq:
        p = M.index(v) + 1
        if p > t: t = p
    nu = m + 1
    for v in Uq:
        p = M.index(v) + 1
        if p < nu: nu = p
    return (nu - t) if nu > t else 0


n = 6
npairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
stats = dict(cases=0, dom_ok=0, dom_bad=0, mu_le_wv=0, mu_gt_wv=0)
wvs = []; mus = []; ratios = []
cfg = {}
tried = 0
while tried < 40000 and time.time() - t0 < 480 and stats['cases'] < 400:
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
        ist, jst = cc['istar'], cc['jstar']
        s_, t_ = an['s'], an['t']
        if ist == 0 or jst == 0 or ist + 1 > s_ or jst + 1 > t_: continue
        Q = [v for v in range(n) if v != z]
        A, B = part
        cells = {}
        okgen = True
        for (nm, kk, ll) in (('m1', ist, jst + 1), ('m2', ist + 1, jst),
                             ('c1', ist, jst), ('c2', ist + 1, jst + 1)):
            x, y = Q[A[kk - 1]], Q[B[ll - 1]]
            if L[x][y] or L[y][x]: okgen = False; break
            cells[nm] = pval(L, n, x, y)
        if not okgen: continue
        m1, m2, c1, c2 = cells['m1'], cells['m2'], cells['c1'], cells['c2']
        # configuration under hypothesis-sign pattern (c1,c2 each: hi/lo/bal)
        def lab(v_):
            if v_ > twothird: return 'hi'
            if v_ < third: return 'lo'
            return 'bal'
        key = (lab(c1), lab(c2))
        cfg[key] = cfg.get(key, 0) + 1
        # WV
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
        Iset = set(A[:ist]) | set(B[:jst])
        wv = sum(gap_of(M, Dq, Uq, m) for M in extsQ if set(M[:ist + jst]) == Iset)
        wv = Fraction(wv, eP)
        mu = an['mu'].get((ist, jst), Fraction(0))
        a1 = c1 - m2; a2 = m1 - c2; b1 = m1 - c1; b2 = c2 - m2
        edges = (a1, a2, b1, b2)
        stats['cases'] += 1
        if all(e_ <= wv for e_ in edges): stats['dom_ok'] += 1
        else: stats['dom_bad'] += 1
        if mu <= wv: stats['mu_le_wv'] += 1
        else: stats['mu_gt_wv'] += 1
        wvs.append(wv); mus.append(mu)
        if wv > 0: ratios.append(mu / wv)

print(stats)
import collections
print("configs (c1lab,c2lab):", dict(sorted(cfg.items())))
print("WV>1/3:", sum(1 for w in wvs if w > third), "/", len(wvs))
print("mu>1/3:", sum(1 for x in mus if x > third), "/", len(mus))
if ratios:
    rmax = max(ratios); rmin = min(ratios)
    avg = sum(ratios) / len(ratios)
    print("mu/WV ratio: min=%s max=%s avg=%.3f" % (rmin, rmax, float(avg)))
print("WV top:", sorted(set(wvs), reverse=True)[:5])
print("mu top:", sorted(set(mus), reverse=True)[:5])
