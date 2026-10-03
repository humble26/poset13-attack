import random
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])

random.seed(11)
found = 0
for t in range(800):
    n = random.choice([4, 5, 6])
    rel = [(i, j) for i in range(n) for j in range(i + 1, n) if random.random() < 0.4]
    L = closure(n, rel)
    if not is_poset(L, n): continue
    m = n - 1
    Q = [v for v in range(n) if v != n - 1]
    Ql, _ = induced(L, n, Q)
    parts = chain_partitions(Ql, m)
    if not parts: continue
    A, B = parts[0]
    if not A or not B: continue
    x = Q[random.choice(A)]; y = Q[random.choice(B)]
    if L[x][y] or L[y][x]: continue
    exts = []
    def rec(placed, seq):
        if len(seq) == m: exts.append(list(seq)); return
        for v in Q:
            if (placed >> Q.index(v)) & 1: continue
            if all((placed >> Q.index(u)) & 1 for u in Q if L[u][v]):
                rec(placed | 1 << Q.index(v), seq + [v])
    rec(0, [])
    bf = sum(1 for e in exts if e.index(x) < e.index(y))
    xf, yf, tot = inter_count(L, [Q[v] for v in A], [Q[v] for v in B], x, y)
    if xf != bf or tot != len(exts):
        print("FAIL n=%d rel=%s x=%d y=%d  DP(xf,tot)=(%d,%d) BF=(%d,%d)" % (n, rel, x, y, xf, tot, bf, len(exts)))
        print("  Q=%s A=%s B=%s" % (Q, A, B))
        print("  exts:", exts)
        print("  L matrix (Q only):")
        for i in Q:
            print("   ", [1 if L[i][j] else 0 for j in Q])
        found += 1
        if found >= 2: break
