import random
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])

random.seed(23)
shown = 0
for t in range(600):
    n = random.choice([4, 5, 6])
    rel = [(i, j) for i in range(n) for j in range(i + 1, n) if random.random() < 0.4]
    L = closure(n, rel)
    if not is_poset(L, n): continue
    z = n - 1
    Q = list(range(n - 1))
    Ql, idx = induced(L, n, Q)
    parts = chain_partitions(Ql, n - 1)
    if not parts: continue
    incs = [(x, y) for x in Q for y in Q if x < y and not (L[x][y] or L[y][x])]
    if not incs: continue
    x, y = incs[0]
    Dg = [v for v in Q if L[v][z]]; Ug = [v for v in Q if L[z][v]]
    extsP = []
    def recP(placed, seq):
        if len(seq) == n: extsP.append(list(seq)); return
        for v in range(n):
            if (placed >> v) & 1: continue
            if all((placed >> u) & 1 for u in range(n) if L[u][v]):
                recP(placed | 1 << v, seq + [v])
    recP(0, [])
    extsQ = [[v for v in e if v != z] for e in extsP]
    m = n - 1
    tot = 0
    for M in extsQ:
        posd = [M.index(v) + 1 for v in Dg]
        posu = [M.index(v) + 1 for v in Ug]
        tt = max(posd) if posd else 0
        nu = min(posu) if posu else m + 1
        tot += (nu - tt) if nu > tt else 0
    eP = count_ext(L, n)
    eQ = count_ext(Ql, n - 1)
    print("n=%d rel=%s z=%d x=%d y=%d D=%s U=%s" % (n, rel, z, x, y, Dg, Ug))
    print("  e(P)=%d  e(Q)=%d  #extQ=%d  sum_gap=%d" % (eP, eQ, len(extsQ), tot))
    print("  extsQ:", extsQ)
    shown += 1
    if shown >= 3: break
