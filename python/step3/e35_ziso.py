"""E35: verify the z-isolated reduction (p_P = p_Q on Q-pairs when D=U=empty),
and dissect the n=7 extremal family (deep-D balance = mu-average of corner races).
"""
import random
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e11_corner.py").resolve(), encoding="utf-8").read().split("def main():")[0])

third = Fraction(1, 3)
random.seed(99)

# Part 1: z-isolated => p_P = p_Q on all Q-pairs
ok = bad = 0
for t_ in range(600):
    n = random.choice([4, 5, 6])
    rel = [(i, j) for i in range(n) for j in range(i + 1, n) if random.random() < 0.35]
    L = closure(n, rel)
    if not is_poset(L, n): continue
    z = n - 1
    if any(L[v][z] or L[z][v] for v in range(n - 1)): continue  # z isolated?
    Qz = [v for v in range(n) if v != z]
    Ql, _ = induced(L, n, Qz)
    eQ = count_ext(Ql, n - 1)
    extsQ = []
    def rec(placed, seq):
        if len(seq) == n - 1: extsQ.append(list(seq)); return
        for v in Qz:
            if (placed >> Qz.index(v)) & 1: continue
            if all((placed >> Qz.index(u)) & 1 for u in Qz if L[u][v]):
                rec(placed | 1 << Qz.index(v), seq + [v])
    rec(0, [])
    incs = [(x, y) for x in Qz for y in Qz if x < y and not (L[x][y] or L[y][x])]
    if not incs: continue
    x, y = random.choice(incs)
    pP = pval(L, n, x, y)
    pQ = Fraction(sum(1 for e in extsQ if e.index(Qz.index(x)) < e.index(Qz.index(y))), eQ)
    if pP == pQ: ok += 1
    else: bad += 1; print("ZISO FAIL", n, x, y, pP, pQ)
print("z-isolated p_P=p_Q: ok=%d bad=%d" % (ok, bad))

# Part 2: the n=7 extremal anatomy — deep-D pair as mu-average of corner races
n = 7
rel = [(0, 1), (0, 3), (0, 4), (0, 6), (1, 3), (2, 4), (2, 6), (3, 6), (4, 5), (4, 6)]
z = 4
L = closure(n, rel)
Q = [v for v in range(n) if v != z]
# chains A=(0,1,3,5), B=(2,6); D={0,2}, U={5,6}; window A={1,3}
a1, b1 = 0, 2
mu = {(1, 1): Fraction(8, 25), (2, 1): Fraction(9, 25), (3, 1): Fraction(8, 25)}
A = [0, 1, 3, 5]; B = [2, 6]
tot = Fraction(0)
for (i, j), w in mu.items():
    Iset = set(A[:i]) | set(B[:j])
    Ilist = sorted(Iset)
    Il, _ = induced(L, n, Ilist)
    pI = Fraction(sum(1 for _ in range(1)), 1) if False else None
    extsI = []
    def rec2(placed, seq):
        if len(seq) == len(Ilist): extsI.append(list(seq)); return
        for v in Ilist:
            if (placed >> v) & 1: continue
            if all((placed >> u) & 1 for u in Ilist if L[u][v]):
                rec2(placed | 1 << v, seq + [v])
    rec2(0, [])
    pI = Fraction(sum(1 for e in extsI if e.index(a1) < e.index(b1)), len(extsI))
    tot += w * pI
    print("cut (%d,%d): I=%s p_I(a1,b1)=%s w=%s" % (i, j, Ilist, pI, w))
pP = pval(L, n, a1, b1)
print("mu-average =", tot, " vs p_P(a1,b1) =", pP, " match:", tot == pP)
print("balanced:", third <= pP <= 1 - third)
