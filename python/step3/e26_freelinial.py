"""E26: verify the NEW 3-case proof of 1/3-2/3 for the free corridor.
The proof CONSTRUCTS the balanced cell; verify each construction step:
 case M (t/2 <= s <= 2t): p(1,1) = s/(s+t) in [1/3,2/3]
 case S > 2T (s > 2t): p(1,1) > 2/3, p(s,1) < 1/3,
     every column step = C(s-k+t-1,t-1)/C(s+t,s) < 1/3, and the column LANDS
     (some p(k,1) in [1/3,2/3])
 case t > 2s: dual (row).
Also verify the closed forms p(k,1) = C(s-k+t,t)/C(s+t,s) and step identity.
"""
from fractions import Fraction
from math import comb

ok = bad = 0
fails = []
for s in range(1, 9):
    for t in range(1, 9):
        m = s + t
        # exact p(k,1) via closed form and via brute force enumeration
        col = [Fraction(comb(s - k + t, t), comb(m, s)) for k in range(1, s + 1)]
        # brute force
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
        colBF = [Fraction(sum(1 for e in exts if e.index(A[k]) < e.index(B[0])), eQ)
                 for k in range(s)]
        if col != colBF:
            bad += 1; fails.append(("closed-form", s, t)); continue
        # steps
        steps = [col[k] - col[k + 1] for k in range(s - 1)]
        case = None
        if Fraction(t, 2) <= s <= 2 * t:
            case = 'M'
            good = Fraction(1, 3) <= col[0] <= Fraction(2, 3)
        elif s > 2 * t:
            case = 'S'
            maxstep = max(steps) if steps else Fraction(0)
            good = col[0] > Fraction(2, 3) and col[-1] < Fraction(1, 3) \
                and maxstep < Fraction(1, 3) \
                and any(Fraction(1, 3) <= c <= Fraction(2, 3) for c in col)
            # verify the step bound identity
            if steps:
                bound_ok = maxstep <= Fraction(s * t, (s + t) * (s + t - 1))
            else:
                bound_ok = True
            good = good and bound_ok
        else:
            case = 'T'
            # t > 2s: row 1 INCREASES from p(1,1)=s/(s+t) < 1/3 to p(1,t) > 2/3
            row = [Fraction(sum(1 for e in exts if e.index(A[0]) < e.index(B[l])), eQ)
                   for l in range(t)]
            rowsteps = [row[l + 1] - row[l] for l in range(t - 1)]
            maxstep = max(rowsteps) if rowsteps else Fraction(0)
            good = row[0] < Fraction(1, 3) and row[-1] > Fraction(2, 3) \
                and maxstep < Fraction(1, 3) \
                and any(Fraction(1, 3) <= c <= Fraction(2, 3) for c in row)
        # ultimate check: some balanced cell exists in the row/col used by the proof
        if case == 'T':
            bal = any(Fraction(1, 3) <= p_ <= Fraction(2, 3) for p_ in row)
        else:
            bal = any(Fraction(1, 3) <= p_ <= Fraction(2, 3) for p_ in col)
        if good and bal: ok += 1
        else:
            bad += 1; fails.append((case, s, t, good, bal))
print("E26: constructed-proof verification ok=%d bad=%d" % (ok, bad))
for f in fails[:6]: print("  FAIL:", f)
