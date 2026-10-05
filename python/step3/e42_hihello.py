"""E42: (m1,m2) fail-pattern frequencies over chain-box fail cases (n<=8).
Patterns: (lo,lo), (lo,hi), (hi,hi). Also verify the pair census:
in (hi,hi) chain-box cases, is EVERY pair really unbalanced (theorem would fail)?
Uses exhaustive_n.py-style enumeration with independent implementation cross-check.
"""
import sys, itertools, collections
from fractions import Fraction
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from my_poset import my_close, my_is_poset, my_induced, my_count_ext

third = Fraction(1, 3)
twothird = Fraction(2, 3)


def width_le2(Q, mm):
    for a, b, c in itertools.combinations(range(mm), 3):
        if not (Q[a][b] or Q[b][a] or Q[a][c] or Q[c][a] or Q[b][c] or Q[c][b]): return False
    return True


def pval_my(Q, mm, x, y):
    L2 = [row[:] for row in Q]
    L2[x][y] = True
    for k in range(mm):
        for i in range(mm):
            if L2[i][k]:
                for j in range(mm):
                    if L2[k][j]: L2[i][j] = True
    return Fraction(my_count_ext(L2, mm), my_count_ext(Q, mm))


def main():
    hist = collections.Counter()
    hihi_cases = []
    n_pair_check = 0
    for n in range(4, 9):
        pairs = list(itertools.combinations(range(n), 2))
        np_ = len(pairs)
        for mask in range(1 << np_):
            rel = [pairs[k] for k in range(np_) if (mask >> k) & 1]
            L = my_close(n, rel)
            if not my_is_poset(L, n): continue
            # need some z with width(Q)<=2 and chain-box: D,U chains, W={w1,w2}, w1<w2
            for z in range(n):
                Qz = [v for v in range(n) if v != z]
                mm = n - 1
                Ql = my_induced(L, Qz)   # relabel in Qz order
                idx = {v: k for k, v in enumerate(Qz)}
                if not width_le2(Ql, mm): continue
                # W = global elements incomparable to z (use full poset L)
                Wz = [v for v in Qz if not (L[v][z] or L[z][v])]
                if len(Wz) != 2: continue
                w1, w2 = Wz
                if Ql[w2][w1]:
                    w1, w2 = w2, w1
                if not Ql[w1][w2]: continue  # need w1 < w2
                # D,U chains? D = {v: v<z}, U = {v: z<v}
                Dz = [v for v in Qz if L[v][z]]
                Uz = [v for v in Qz if L[z][v]]
                def is_chain(els):
                    return all(Ql[a][b] or Ql[b][a] for a, b in itertools.combinations(els, 2)) or len(els) <= 1
                if not (is_chain(Dz) and is_chain(Uz)): continue
                # m1, m2 via p-values
                zz = z
                m1 = pval_my(L, n, zz, Q[w1])
                m2 = pval_my(L, n, zz, Q[w2])
                def lab(v):
                    if v < third: return 'lo'
                    if v > twothird: return 'hi'
                    return 'bal'
                p1, p2 = lab(m1), lab(m2)
                hist[(p1, p2)] += 1
                if p1 == 'hi' and p2 == 'hi':
                    hihi_cases.append((n, z, rel, str(m1), str(m2)))
                n_pair_check += 1
    print("chain-box fail-relevant instances:", n_pair_check)
    print("(m1,m2) patterns:", dict(hist))
    print("hi-hi cases:", len(hihi_cases))
    for c_ in hihi_cases[:3]: print("   HIHI:", c_[0], "z=", c_[1], "m=", c_[3], c_[4])


if __name__ == '__main__':
    main()
