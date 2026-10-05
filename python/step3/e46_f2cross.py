"""E46: cross-tabulate the F2 sweep-failure instances (all d*W1 pairs > 2/3)
by fail pattern (m1, m2 labels) and parameter structure. Also check whether
the (w_k, u) U-side pairs rescue them (balanced)."""
import random, collections
from fractions import Fraction
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from my_poset import my_close, my_is_poset, my_induced, my_count_ext

third = Fraction(1, 3)
twothird = Fraction(2, 3)


def build(r, a1, a2, sp, b1, b2):
    n = r + 2 + sp + 1
    z = n - 1
    w1, w2 = r, r + 1
    u = list(range(r + 2, z))
    rel = [(i, i + 1) for i in range(r - 1)]
    rel.append((w1, w2))
    rel += [(u[k], u[k + 1]) for k in range(sp - 1)]
    for j in range(a1 - 1): rel.append((j, w1))
    for j in range(a2 - 1): rel.append((j, w2))
    for j in range(b1, sp): rel.append((w1, u[j]))
    for j in range(b2, sp): rel.append((w2, u[j]))
    for d in range(r): rel.append((d, z))
    for uu in u: rel.append((z, uu))
    return my_close(n, rel), n, w1, w2, u, z


def pval(L, n, x, y):
    L2 = [row[:] for row in L]
    L2[x][y] = True
    for k in range(n):
        for i in range(n):
            if L2[i][k]:
                for j in range(n):
                    if L2[k][j]: L2[i][j] = True
    return Fraction(my_count_ext(L2, n), my_count_ext(L, n))


def main():
    cross = collections.Counter()
    rescue = collections.Counter()
    n_f2 = 0
    examples = []
    for r in range(1, 7):
        for a1 in range(1, r + 2):
            for a2 in range(a1, r + 2):
                for sp in range(0, 5):
                    for b1 in range(0, sp + 1):
                        for b2 in range(b1, sp + 1):
                            L, n, w1, w2, u, z = build(r, a1, a2, sp, b1, b2)
                            if not my_is_poset(L, n): continue
                            m1 = pval(L, n, z, w1)
                            m2 = pval(L, n, z, w2)
                            # sweep all > 2/3?
                            allhi = True
                            for i in range(a1, r + 1):
                                d = i - 1
                                if not (third < pval(L, n, d, w1) <= 1): allhi = False; break
                            if not allhi: continue
                            # this is an F2 instance; classify fail pattern
                            def lab(v):
                                if v < third: return 'lo'
                                if v > twothird: return 'hi'
                                return 'bal'
                            key = (lab(m1), lab(m2))
                            n_f2 += 1
                            cross[key] += 1
                            # check (w_k, u) rescue: w1/u and w2/u incomparable pairs
                            ok = False
                            for w in (w1, w2):
                                for uu in u:
                                    if L[w][uu] or L[uu][w]: continue
                                    p = pval(L, n, w, uu)
                                    if third <= p <= twothird: ok = True; break
                                if ok: break
                            rescue['yes' if ok else 'no'] += 1
                            if not ok and len(examples) < 3:
                                examples.append((r, a1, a2, sp, b1, b2, str(m1), str(m2)))
    print("F2 instances:", n_f2)
    print("fail patterns:", dict(cross))
    print("(w,u) rescue:", dict(rescue))
    for e in examples: print("  NO-RESCUE:", e)


if __name__ == '__main__':
    main()
