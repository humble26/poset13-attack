"""E29: verify the edge/vertex mass identities on free corridors.
alpha1 = P(b_{j*} between a_{i*}, a_{i*+1}) = up-edge mass at (i*,j*-1)->(i*,j*)
       = C(i*+j*-1, i*) * C(s-i*+t-j*, s-i*) / C(s+t, s)
V(i*,j*) = vertex mass at (i*,j*) = C(i*+j*, i*) * C(s-i*+t-j*, s-i*) / C(s+t,s)
alpha1 <= V (using an edge implies passing its head vertex).
Also beta1 = left-edge mass into (i*,j*) and beta1 <= V.
"""
from fractions import Fraction
from math import comb

ok = bad = 0
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
        for ist in range(1, s):
            for jst in range(1, t):
                a1 = Fraction(sum(1 for e in exts
                                  if e.index(A[ist-1]) < e.index(B[jst-1]) < e.index(A[ist])), eQ)
                a2 = Fraction(sum(1 for e in exts
                                  if e.index(A[ist-1]) < e.index(B[jst]) < e.index(A[ist])), eQ)
                b1 = Fraction(sum(1 for e in exts
                                  if e.index(B[jst-1]) < e.index(A[ist-1]) < e.index(B[jst])), eQ)
                b2 = Fraction(sum(1 for e in exts
                                  if e.index(B[jst-1]) < e.index(A[ist]) < e.index(B[jst])), eQ)
                up1 = Fraction(comb(ist + jst - 1, ist) * comb(s - ist + t - jst, s - ist), comb(m, s))
                up2b = Fraction(comb(ist + jst, ist) * comb(s - ist + t - jst - 1, s - ist), comb(m, s))
                # beta1: E-step into column ist at height jst: edge (ist-1,jst)->(ist,jst) [in-left]
                left1 = Fraction(comb(ist - 1 + jst, ist - 1) * comb(s - ist + t - jst, s - ist), comb(m, s))
                # beta2: E-step into column ist+1 at height jst: edge (ist,jst)->(ist+1,jst) [out-right]
                left2 = Fraction(comb(ist + jst, ist) * comb(s - ist - 1 + t - jst, s - ist - 1), comb(m, s))
                V = Fraction(comb(ist + jst, ist) * comb(s - ist + t - jst, s - ist), comb(m, s))
                conds = [a1 == up1, a2 == up2b, b1 == left1, b2 == left2,
                         a1 <= V, a2 <= V, b1 <= V, b2 <= V]
                if all(conds): ok += 1
                else:
                    bad += 1
                    if bad <= 3: print("FAIL", s, t, ist, jst, a1, up1, a2, up2b, b1, left1, b2, left2, V)
print("edge/vertex identities: ok=%d bad=%d" % (ok, bad))
