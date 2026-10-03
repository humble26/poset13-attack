"""E14: n=6 sampled verification of the empirical claims.
Claims on fail cases (z-pairs all unbalanced), per (case,partition):
  A) mixed corner cells strictly oriented: p(a_{i*},b_{j*+1}) > 2/3, p(a_{i*+1},b_{j*}) < 1/3
  B) same-side corner cells balanced when incomparable (C1) -- count violations
  C) conservation: some Q-pair balanced (should be 100%)
  D) rescue pair category: same-side (lo,lo)/(hi,hi) only?
"""
import random, sys
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e11_corner.py").resolve(), encoding="utf-8").read().split("def main():")[0])

third = Fraction(1, 3); twothird = Fraction(2, 3)
random.seed(2026)

NSAMP = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
n = 6
npairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
agg = dict(cases=0, fail=0, parts=0, conserv_bad=0,
           A01=0, A01_ok=0, A10=0, A10_ok=0,
           c00_inc=0, c00_bad=0, c11_inc=0, c11_bad=0,
           rescue_side=dict(lo=0, hi=0, mixed=0, none=0))
badexamples = []
tried = 0
while tried < NSAMP:
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
    agg['cases'] += 1
    if r['zpair']: continue
    agg['fail'] += 1
    for part in chain_partitions(Ql, n - 1):
        cc = corner_check(L, n, z, part)
        if cc is None: continue
        agg['parts'] += 1
        an = cc['an']
        # C) conservation
        if not an['bal_pairs']:
            agg['conserv_bad'] += 1
        cells = cc['cells']
        # A) mixed cells
        c01 = cells.get((0, 1)); c10 = cells.get((1, 0))
        if c01 and c01[0] == 'inc':
            agg['A01'] += 1
            if c01[1] > twothird: agg['A01_ok'] += 1
            elif len(badexamples) < 8: badexamples.append(('A01', n, z, rel, part, str(c01[1])))
        if c10 and c10[0] == 'inc':
            agg['A10'] += 1
            if c10[1] < third: agg['A10_ok'] += 1
            elif len(badexamples) < 8: badexamples.append(('A10', n, z, rel, part, str(c10[1])))
        # B) same-side corners
        c00 = cells.get((0, 0)); c11 = cells.get((1, 1))
        if c00 and c00[0] == 'inc':
            agg['c00_inc'] += 1
            if not (third <= c00[1] <= twothird):
                agg['c00_bad'] += 1
                if len(badexamples) < 8: badexamples.append(('c00', n, z, rel, part, str(c00[1])))
        if c11 and c11[0] == 'inc':
            agg['c11_inc'] += 1
            if not (third <= c11[1] <= twothird):
                agg['c11_bad'] += 1
                if len(badexamples) < 8: badexamples.append(('c11', n, z, rel, part, str(c11[1])))
        # D) balanced-pair sides
        sides = set()
        for (ba, bb, pp, cxA, cyB, lx, hx, ly, hy, k, l) in an['bal_pairs']:
            fx = 'lo' if (cxA == 'D' or (cxA == 'W' and an['MAv'][k] < third)) else 'hi'
            fy = 'lo' if (cyB == 'D' or (cyB == 'W' and an['MBv'][l] < third)) else 'hi'
            sides.add(fx + fy)
        if not sides: agg['rescue_side']['none'] += 1
        elif sides <= {'lolo'}: agg['rescue_side']['lo'] += 1
        elif sides <= {'hihi'}: agg['rescue_side']['hi'] += 1
        else: agg['rescue_side']['mixed'] += 1

print(agg)
print("bad examples:", len(badexamples))
for b in badexamples:
    print(b)
