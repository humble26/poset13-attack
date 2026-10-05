"""E47: the (lo,hi) sandwich landing test (M1's last step for the chain box).
For all (r, a1, a2, sp, b1, b2) with fail pattern (lo,hi): m1 < 1/3, m2 > 2/3:
check that SOME (w2, u_j) pair (j in the w2-incomparable U-prefix) has
p_P(w2, u_j) in [1/3, 2/3].  Also verify the race closed form:
race(j) = #{(x,y): x<=b1, j<=x+y<=b2} / U0  (F0-extensions with w2 before u_j).
"""
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
    st = dict(cases=0, landed=0, notlanded=0, race_ok=0, race_bad=0)
    fails = []
    for r in range(1, 7):
        for a1 in range(1, r + 2):
            for a2 in range(a1, r + 2):
                for sp in range(0, 6):
                    for b1 in range(0, sp + 1):
                        for b2 in range(b1, sp + 1):
                            L, n, w1, w2, u, z = build(r, a1, a2, sp, b1, b2)
                            if not my_is_poset(L, n): continue
                            m1 = pval(L, n, z, w1)
                            m2 = pval(L, n, z, w2)
                            if not (m1 < third and m2 > twothird): continue
                            st['cases'] += 1
                            # race closed form: #{(x,y): x<=b1, j<=x+y<=b2}/U0
                            U0 = (b1 + 1) * (b2 + 1) - b1 * (b1 + 1) // 2
                            landed = False
                            for j in range(1, b2 + 1):
                                uj = u[j - 1]
                                if L[w2][uj] or L[uj][w2]: continue  # not incomparable
                                # direct race in F0 = induced on Q∖D = {w1,w2} ∪ U
                                QD = [v for v in range(n) if v != z and v != u[j - 1] and
                                      not (v == u[k - 1] and False for k in range(0))]
                                # simpler: F0 minus u_j for the race: count extensions of F0 with w2 before u_j via added relation on F0
                                QD = [v for v in range(n) if v != z]
                                F0 = my_induced(L, QD)
                                w1f, w2f, ujf = QD.index(w1), QD.index(w2), QD.index(uj)
                                F0r = [row[:] for row in F0]
                                F0r[w2f][ujf] = True
                                for k in range(len(QD)):
                                    for i in range(len(QD)):
                                        if F0r[i][k]:
                                            for jj in range(len(QD)):
                                                if F0r[k][jj]: F0r[i][jj] = True
                                ef = my_count_ext(F0, len(QD))
                                race = Fraction(my_count_ext(F0r, len(QD)), ef)
                                # closed form check: #{(x,y): x<=b1, j<=x+y<=b2}/U0
                                cc = sum(1 for x in range(0, b1 + 1) for y in range(0, sp - x)
                                         if j <= x + y + 1 <= b2)  # x 0-based u's before w1; u_j index j (1-based)
                                cc2 = sum(1 for x in range(0, b1 + 1) for y in range(0, b2 - x)
                                          if x + y + 1 >= j)
                                if cc != cc2: pass
                                race_cf = Fraction(cc2, U0)
                                if race == race_cf: st['race_ok'] += 1
                                else: st['race_bad'] += 1
                                pP = pval(L, n, w2, uj)
                                # p_P(w2,u_j) = (1-m2) + m2*race
                                pred = (1 - m2) + m2 * race
                                if pP != pred:
                                    st['race_bad'] += 1
                                    if len(fails) < 4: fails.append(("pred", r, a1, a2, sp, b1, b2, str(pP), str(pred)))
                                if third <= pP <= twothird: landed = True
                            if landed: st['landed'] += 1
                            else:
                                st['notlanded'] += 1
                                if len(fails) < 8: fails.append(("NOT-LANDED", r, a1, a2, sp, b1, b2, str(m1), str(m2)))
    print(st)
    for f in fails: print("  ", f)


if __name__ == '__main__':
    main()
