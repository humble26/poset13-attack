"""E33: comprehensive verification of GCB (generic corner balance).
Claim: in fail cases, if the 2x2 corner cells m1,m2,c1,c2 all exist and are
mutually incomparable (generic), then BOTH c1=(a_{i*},b_{j*}) and
c2=(a_{i*+1},b_{j*+1}) are balanced.
Test: n=4,5 exhaustive + n=6,7 sampled. Also record WV/mu constants.
"""
import random
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e11_corner.py").resolve(), encoding="utf-8").read().split("def main():")[0])

third = Fraction(1, 3); twothird = Fraction(2, 3)
random.seed(2024)


def gcb_check(L, n, z, part):
    cc = corner_check(L, n, z, part)
    if cc is None: return None
    an = cc['an']
    ist, jst = cc['istar'], cc['jstar']
    s_, t_ = an['s'], an['t']
    if ist == 0 or jst == 0 or ist + 1 > s_ or jst + 1 > t_: return None
    Q = [v for v in range(n) if v != z]
    A, B = part
    cells = {}
    for (nm, kk, ll) in (('m1', ist, jst + 1), ('m2', ist + 1, jst),
                         ('c1', ist, jst), ('c2', ist + 1, jst + 1)):
        x, y = Q[A[kk - 1]], Q[B[ll - 1]]
        if L[x][y] or L[y][x]: return None  # not generic
        cells[nm] = pval(L, n, x, y)
    bal1 = third <= cells['c1'] <= twothird
    bal2 = third <= cells['c2'] <= twothird
    return (bal1, bal2, cells)


agg = dict(n5=dict(gen=0, both=0, only1=0, only2=0, none=0),
           n6=dict(gen=0, both=0, only1=0, only2=0, none=0),
           n7=dict(gen=0, both=0, only1=0, only2=0, none=0))
wv_mu = []
for n in (4, 5):
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
                res = gcb_check(L, n, z, part)
                if res is None: continue
                b1, b2, cells = res
                key = 'n%d' % n
                agg[key]['gen'] += 1
                if b1 and b2: agg[key]['both'] += 1
                elif b1: agg[key]['only1'] += 1
                elif b2: agg[key]['only2'] += 1
                else: agg[key]['none'] += 1
for n, nsamp in ((6, 20000), (7, 8000)):
    npairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    tried = 0
    while tried < nsamp:
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
            res = gcb_check(L, n, z, part)
            if res is None: continue
            b1, b2, cells = res
            key = 'n%d' % n
            agg[key]['gen'] += 1
            if b1 and b2: agg[key]['both'] += 1
            elif b1: agg[key]['only1'] += 1
            elif b2: agg[key]['only2'] += 1
            else: agg[key]['none'] += 1
for key in ('n4', 'n5', 'n6', 'n7'):
    a = agg.get(key)
    if a: print(key, "generic:", a['gen'], "| both-bal:", a['both'],
                "only-c1:", a['only1'], "only-c2:", a['only2'], "none:", a['none'])
