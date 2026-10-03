"""E23 (clean): free-corridor facts.
 (a) corner values: p(1,t) = 1 - 1/C(s+t,s), p(s,1) = 1/C(s+t,s)
 (b) center cell (1,1) balanced for s=t=2
 (c) Vandermonde level budget: sum_k w^{(k)}_h = s/(t+1) for every h,
     w^{(k)}_h = P(exactly h B-elements before a_k)
"""
from fractions import Fraction
from math import comb

ok_a = bad_a = ok_b = bad_b = ok_c = bad_c = 0
for s in range(2, 6):
    for t in range(2, 6):
        m = s + t
        less = [[False] * m for _ in range(m)]
        for i in range(s):
            for j in range(i + 1, s): less[i][j] = True
        for i in range(s, m):
            for j in range(s, m):
                if i < j: less[i][j] = True
        A = list(range(s)); B = list(range(s, m))
        exts = []
        def rec(placed, seq):
            if len(seq) == m: exts.append(list(seq)); return
            for v in range(m):
                if (placed >> v) & 1: continue
                if all((placed >> u) & 1 for u in range(m) if less[u][v]):
                    rec(placed | 1 << v, seq + [v])
        rec(0, [])
        eQ = len(exts)
        assert eQ == comb(m, s)
        # (a)
        p1t = Fraction(sum(1 for e in exts if e.index(A[0]) < e.index(B[t - 1])), eQ)
        ps1 = Fraction(sum(1 for e in exts if e.index(A[s - 1]) < e.index(B[0])), eQ)
        if p1t == 1 - Fraction(1, comb(m, s)) and ps1 == Fraction(1, comb(m, s)):
            ok_a += 1
        else:
            bad_a += 1; print("corner FAIL", s, t, p1t, ps1)
        # (b)
        if s == 2 and t == 2:
            p11 = Fraction(sum(1 for e in exts if e.index(A[0]) < e.index(B[0])), eQ)
            if Fraction(1, 3) <= p11 <= Fraction(2, 3): ok_b += 1
            else: bad_b += 1; print("center FAIL", s, t, p11)
        # (c)
        for h in range(t + 1):
            tot = sum(1 for e in exts
                      for k in range(1, s + 1)
                      if sum(1 for b in B if e.index(b) < e.index(A[k - 1])) == h)
            if tot == Fraction(s, t + 1) * eQ: ok_c += 1
            else: bad_c += 1; print("vandermonde FAIL", s, t, h, Fraction(tot, eQ))
print("corner ok/bad: %d/%d | center ok/bad: %d/%d | vandermonde ok/bad: %d/%d"
      % (ok_a, bad_a, ok_b, bad_b, ok_c, bad_c))
