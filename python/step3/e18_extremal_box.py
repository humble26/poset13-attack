"""E18: structure of extremal (delta = min) cases at n=6.
Check: do all delta==1/3 target-class cases have a single-point box (z extreme)?
"""
import itertools, sys
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])

third = Fraction(1, 3)


def delta_of(less, n):
    best = Fraction(0); arg = None
    for x, y in itertools.combinations(range(n), 2):
        if less[x][y] or less[y][x]: continue
        p = pval(less, n, x, y)
        m = min(p, 1 - p)
        if m > best: best = m; arg = (x, y, p)
    return best, arg


n = 6
npairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
all_cases = []
for mask in range(1 << len(npairs)):
    rel = [npairs[k] for k in range(len(npairs)) if (mask >> k) & 1]
    L = closure(n, rel)
    if not is_poset(L, n): continue
    if antichain_width(L, n) < 2: continue
    for z in range(n):
        Qz = [v for v in range(n) if v != z]
        Ql, _ = induced(L, n, Qz)
        if antichain_width(Ql, n - 1) > 2: continue
        d, arg = delta_of(L, n)
        all_cases.append((d, tuple(rel), z))
mn = min(c[0] for c in all_cases)
hard = [c for c in all_cases if c[0] == mn]
print("n=6: cases=%d min delta=%s hard=%d" % (len(all_cases), mn, len(hard)))
boxsizes = {}
single = 0
for (d, rel, z) in hard:
    L = closure(n, list(rel))
    Qz = [v for v in range(n) if v != z]
    Ql, _ = induced(L, n, Qz)
    parts = chain_partitions(Ql, n - 1)
    sizes = set()
    for part in parts:
        an = anatomy(L, n, z, part)
        if an is None: continue
        sizes.add(len(an['mu']))
    key = min(sizes) if sizes else -1
    boxsizes[key] = boxsizes.get(key, 0) + 1
    if key == 1: single += 1
print("hard cases by min-box-size over partitions:", boxsizes)
print("single-cut (box=1):", single, "/", len(hard))
# also: how many hard cases have width(P) == 2 (i.e. P itself width-2)?
w2 = 0
for (d, rel, z) in hard:
    L = closure(n, list(rel))
    if antichain_width(L, n) == 2: w2 += 1
print("hard cases with width(P)==2:", w2, "/", len(hard))
