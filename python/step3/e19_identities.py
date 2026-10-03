"""E19: verify the two identities the session-3 analysis rests on.
(1) Delta identity: for a_k < a_{k+1} on a chain and any b_l:
    p(a_k,b_l) - p(a_{k+1},b_l) = P_L(a_k < b_l < a_{k+1}), L uniform ext(P).
(2) M-frame reweighting (Lemma C specialization):
    p_P(x,y) = sum_M gap(M)*[x<_M y] / sum_M gap(M), M in ext(Q),
    gap(M) = (#legal z-slots) = minU(M) - t(M) if all-D-before-all-U else 0,
    t = last position of D(z), minU = first position of U(z) (1-based, m+1 if none).
"""
import random
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])

random.seed(23)
ok1 = bad1 = ok2 = bad2 = 0
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
    A, B = parts[0]
    Ag = [Q[v] for v in A]; Bg = [Q[v] for v in B]
    Dg = [v for v in Q if L[v][z]]; Ug = [v for v in Q if L[z][v]]
    extsP = []
    def recP(placed, seq):
        if len(seq) == n: extsP.append(list(seq)); return
        for v in range(n):
            if (placed >> v) & 1: continue
            if all((placed >> u) & 1 for u in range(n) if L[u][v]):
                recP(placed | 1 << v, seq + [v])
    recP(0, [])
    extsQ = []
    def recQ(placed, seq):
        if len(seq) == n - 1: extsQ.append(list(seq)); return
        for v in Q:
            if (placed >> idx[v]) & 1: continue
            if all((placed >> idx[u]) & 1 for u in Q if L[u][v]):
                recQ(placed | 1 << idx[v], seq + [v])
    recQ(0, [])
    # (1) delta identity on adjacent chain pair of A
    if len(Ag) >= 2 and Bg:
        k1 = random.randrange(len(Ag) - 1)
        ak, ak1 = Ag[k1], Ag[k1 + 1]
        bl = random.choice(Bg)
        if not (L[ak][bl] or L[bl][ak]) and not (L[ak1][bl] or L[bl][ak1]):
            p1 = pval(L, n, ak, bl); p2 = pval(L, n, ak1, bl)
            between = sum(1 for e in extsP
                          if e.index(ak) < e.index(bl) < e.index(ak1))
            lhs = p1 - p2
            rhs = Fraction(between, len(extsP))
            if lhs == rhs: ok1 += 1
            else: bad1 += 1; print("DELTA FAIL", n, ak, ak1, bl, lhs, rhs)
    # (2) M-frame reweighting for an incomparable pair in Q
    incs = [(x, y) for x in Q for y in Q
            if x < y and not (L[x][y] or L[y][x])]
    if not incs: continue
    x, y = random.choice(incs)
    m = n - 1
    def gap(M):
        posd = [M.index(v) + 1 for v in Dg if v in M]
        posu = [M.index(v) + 1 for v in Ug if v in M]
        t = max(posd) if posd else 0
        nu = min(posu) if posu else m + 1
        return nu - t if nu > t else 0
    num = den = 0
    for M in extsQ:
        g = gap(M)
        den += g
        if M.index(x) < M.index(y): num += g
    lhs = pval(L, n, x, y)
    rhs = Fraction(num, den) if den else None
    if rhs is not None and lhs == rhs: ok2 += 1
    else: bad2 += 1; print("MFRAME FAIL", n, x, y, lhs, rhs, den)
print("delta identity: ok=%d bad=%d | m-frame reweighting: ok=%d bad=%d" % (ok1, bad1, ok2, bad2))
