"""E28: under fail, classify the 2x2 corner configurations (A)-(D) and the
forced separation inequalities. Generic case: all four cells
m1=(i*,j*+1), m2=(i*+1,j*), c1=(i*,j*), c2=(i*+1,j*+1) exist and are incomparable.
Forced under the contradiction hypothesis (all Q-pairs unbalanced):
 (A) c1>2/3 & c2<1/3 => alpha1=p(c1)-p(m2)>1/3 and alpha2=p(m1)-p(c2)>1/3
 (B) c1<1/3 & c2>2/3 => beta1=p(m1)-p(c1)>1/3 and beta2=p(c2)-p(m2)>1/3
 (C) c1>2/3 & c2>2/3 => alpha1>1/3 and beta2>1/3
 (D) c1<1/3 & c2<1/3 => beta1>1/3 and alpha2>1/3
Also: z-separations zeta_A=colmass(i*), zeta_B=rowmass(j*) (>1/3 by Lemma G).
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
random.seed(808)

import collections
cfg_hist = collections.Counter()
viol = []
n_case = 0
for n in (5, 6):
    npairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    masks = range(1 << len(npairs)) if n == 5 else None
    if n == 5:
        mask_iter = list(range(1 << len(npairs)))
        zs = lambda L: range(n)
    else:
        mask_iter = (random.randrange(1 << len(npairs)) for _ in range(20000))
        zs = lambda L: [random.randrange(n)]
    for mask in mask_iter:
        rel = [npairs[k] for k in range(len(npairs)) if (mask >> k) & 1]
        L = closure(n, rel)
        if not is_poset(L, n): continue
        if antichain_width(L, n) < 2: continue
        for z in zs(L):
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
                if ist + 1 > s_ or jst + 1 > t_ or ist == 0 or jst == 0: continue
                Q = [v for v in range(n) if v != z]
                A, B = part
                cells = {}
                okgen = True
                for (nm, kk, ll) in (('m1', ist, jst + 1), ('m2', ist + 1, jst),
                                     ('c1', ist, jst), ('c2', ist + 1, jst + 1)):
                    gx, gy = Q[A[kk - 1]], Q[B[ll - 1]]
                    if L[gx][gy] or L[gy][gx]: okgen = False; break
                    cells[nm] = pval(L, n, gx, gy)
                if not okgen: continue
                n_case += 1
                m1, m2, c1, c2 = cells['m1'], cells['m2'], cells['c1'], cells['c2']
                # sanity: pinned
                if not (m1 > twothird and m2 < third):
                    viol.append(('pin', n, z, str(m1), str(m2))); continue
                if c1 > twothird and c2 < third: cfg = 'A'
                elif c1 < third and c2 > twothird: cfg = 'B'
                elif c1 > twothird: cfg = 'C'
                elif c2 > twothird: cfg = 'D'
                else: cfg = '??'; viol.append(('??', n, z, str(c1), str(c2))); continue
                cfg_hist[cfg] += 1
                # forced inequalities
                a1 = c1 - m2; a2 = m1 - c2
                b1 = m1 - c1; b2 = c2 - m2
                chk = {
                    'A': a1 > third and a2 > third,
                    'B': b1 > third and b2 > third,
                    'C': a1 > third and b2 > third,
                    'D': b1 > third and a2 > third,
                }[cfg]
                if not chk: viol.append(('ineq', n, z, cfg, str(a1), str(a2), str(b1), str(b2)))
                # zeta heaviness
                ca = an['colmass'].get(ist, Fraction(0))
                rb = an['rowmass'].get(jst, Fraction(0))
                if ca <= third or rb <= third:
                    viol.append(('heavy', n, z, str(ca), str(rb)))
print("generic corner cases:", n_case)
print("configurations:", dict(cfg_hist))
print("violations:", len(viol))
for v in viol[:6]: print("  ", v)
