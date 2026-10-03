"""E15: test the SMOOTHING conjecture.
Adding an element w with w > x and w > y (above both) moves p(x,y) toward 1/2:
  |p_{P+w}(x,y) - 1/2| <= |p_P(x,y) - 1/2|
Dually for below-both. Also test the stronger exact-monotonicity of distance.
"""
import random
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])

random.seed(7)
half = Fraction(1, 2)
ok_ab = bad_ab = ok_bl = bad_bl = 0
strict = 0
badex = []
for t in range(4000):
    n = random.choice([4, 5, 6])
    rel = [(i, j) for i in range(n) for j in range(i + 1, n) if random.random() < 0.35]
    L = closure(n, rel)
    if not is_poset(L, n): continue
    # find incomparable pair x,y
    incs = [(x, y) for x in range(n) for y in range(x + 1, n)
            if not (L[x][y] or L[y][x])]
    if not incs: continue
    x, y = random.choice(incs)
    p0 = pval(L, n, x, y)
    d0 = abs(p0 - half)
    if d0 == 0: continue
    mode = random.choice(['above', 'below'])
    # build P + w
    n1 = n + 1
    rel1 = [(i, j) for i in range(n) for j in range(n) if L[i][j]]
    if mode == 'above':
        rel1 += [(x, n), (y, n)]
    else:
        rel1 += [(n, x), (n, y)]
    L1 = closure(n1, rel1)
    if not is_poset(L1, n1): continue
    p1 = pval(L1, n1, x, y)
    d1 = abs(p1 - half)
    if mode == 'above':
        if d1 <= d0: ok_ab += 1; strict += (d1 < d0)
        else:
            bad_ab += 1
            if len(badex) < 6: badex.append(('above', n, x, y, str(p0), str(p1)))
    else:
        if d1 <= d0: ok_bl += 1; strict += (d1 < d0)
        else:
            bad_bl += 1
            if len(badex) < 6: badex.append(('below', n, x, y, str(p0), str(p1)))
print("above-both: ok=%d bad=%d | below-both: ok=%d bad=%d | strict=%d" %
      (ok_ab, bad_ab, ok_bl, bad_bl, strict))
for b in badex: print(b)
