"""E41: LP study of the chain-box sweep lemma.
Sweep values: s_i = m1 + mu1*A(i) + mu2*B(i), i in [a1, r], with
  A(i) = (r-i+1)/N1, N1 = r-a1+2
  B(i) = cnt(i)/e2 (closed forms from e_9e_chainbox)
Weights: m1+mu1+mu2 = 1, all > 0.
Fail: m1 not in [1/3,2/3] AND m1+mu1 not in [1/3,2/3].
Sweep-failure patterns (monotone decreasing sequence):
  F1: all s_i < 1/3   F2: all s_i > 2/3   F3: exists j: s_{j-1} > 2/3, s_j < 1/3
Question: for which (r,a1,a2) is (fail + sweep-failure) LP-feasible?
If infeasible for ALL parameters => sweep lemma is pure algebra (Farkas).
"""
import itertools
from fractions import Fraction
try:
    from scipy.optimize import linprog
    HAS_LP = True
except ImportError:
    HAS_LP = False


def A_of(i, r, a1):
    return Fraction(r - i + 1, r - a1 + 2)


def B_profile(r, a1, a2):
    """B(i) = cnt(i)/e2 for i in [a1, r]; returns (list of (i, B), e2)."""
    e2 = (a2 - a1) * (r - a2 + 2) + (r - a2 + 2) * (r - a2 + 3) // 2
    out = []
    for i in range(a1, r + 1):
        if i >= a2:
            cnt = (r - i + 1) * (r - i + 2) // 2
        else:
            cnt = (a2 - 1 - i) * (r - a2 + 2) + (r - a2 + 2) * (r - a2 + 3) // 2
        out.append((i, Fraction(cnt, e2)))
    return out, e2


def frac_lp(constraints, n=3):
    """Feasibility of {a·x <= b} over x >= 0 (fractions, exact Fourier-Motzkin for n=3)."""
    # variables (m1, mu1, mu2); eliminate mu2 = 1 - m1 - mu1 by substitution:
    # work in (m1, mu1) >= 0, m1 + mu1 <= 1.
    # constraints given as (coeffs2, b) meaning a1*m1 + a2*mu1 <= b (after substitution).
    # exact vertex enumeration over pairs of constraints (small n).
    import itertools as it
    verts = [(Fraction(0), Fraction(0))]
    # candidate vertices: intersections of constraint lines + axes
    lines = constraints[:]  # (a1, a2, b): a1*m1 + a2*mu1 <= b
    lines.append((Fraction(1), Fraction(0), Fraction(0)))   # m1 >= 0 => -m1 <= 0
    lines.append((Fraction(-1), Fraction(0), Fraction(0)))
    lines.append((Fraction(0), Fraction(1), Fraction(0)))   # mu1 >= 0
    lines.append((Fraction(0), Fraction(-1), Fraction(0)))
    lines.append((Fraction(1), Fraction(1), Fraction(1)))   # m1+mu1 <= 1
    delta = Fraction(1, 2000)                                # all three cut masses strictly positive
    lines.append((Fraction(-1), Fraction(0), Fraction(-delta)))
    lines.append((Fraction(0), Fraction(-1), Fraction(-delta)))
    lines.append((Fraction(1), Fraction(1), Fraction(1) - delta))
    def feasible(pt):
        for (a1, a2, b) in lines:
            if a1 * pt[0] + a2 * pt[1] > b: return False
        return True
    cands = [(Fraction(0), Fraction(0)), (Fraction(1), Fraction(0)), (Fraction(0), Fraction(1))]
    for (l1, l2) in it.combinations(lines, 2):
        (a1, a2, b1), (c1, c2, b2) = l1, l2
        det = a1 * c2 - a2 * c1
        if det == 0: continue
        x = (b1 * c2 - a2 * b2) / det
        y = (a1 * b2 - b1 * c2) / det
        cands.append((x, y))
    for pt in cands:
        if feasible(pt): return True
    return False


def sweep_fail_lp(r, a1, a2):
    """Check 12 disjunctive LPs; return list of feasible (failcase, pattern) combos."""
    Bp, e2 = B_profile(r, a1, a2)
    Is = [i for i, _ in Bp]
    res = []
    for m1case in ('lo', 'hi'):            # m1 < 1/3 or m1 > 2/3
        for m2case in ('lo', 'hi'):        # m1+mu1 < 1/3 or > 2/3
            # skip impossible: m1>2/3 & m1+mu1<1/3
            if m1case == 'hi' and m2case == 'lo': continue
            for pattern in ('F1', 'F2', 'F3min', 'F3any'):
                cons = []
                if m1case == 'lo': cons.append((Fraction(1), Fraction(0), Fraction(1)/3))
                else: cons.append((Fraction(-1), Fraction(0), Fraction(-2)/3))
                if m2case == 'lo': cons.append((Fraction(1), Fraction(1), Fraction(1)/3))
                else: cons.append((Fraction(-1), Fraction(-1), Fraction(-2)/3))
                # sweep constraints: s_i = m1 + mu1*A(i) + mu2*B(i), mu2 = 1-m1-mu1
                # s_i = m1 + mu1*A + (1-m1-mu1)*B = B + m1*(1-B) + mu1*(A-B)
                def s_cons(i, sign, bnd):
                    B = dict(Bp)[i]; A = A_of(i, r, a1)
                    # sign=+1: s_i <= bnd ; sign=-1: s_i >= bnd
                    # s_i = B + m1*(1-B) + mu1*(A-B)
                    return (sign * (1 - B), sign * (A - B), sign * (bnd - B))
                if pattern == 'F1':
                    for i in Is: cons.append(s_cons(i, +1, Fraction(1)/3))
                elif pattern == 'F2':
                    for i in Is: cons.append(s_cons(i, -1, Fraction(2)/3))
                elif pattern == 'F3min':
                    # exists j: s_{j-1} > 2/3, s_j < 1/3 (j = a1 means no left cell:
                    # then F3 needs s_{a1-1} = 1 > 2/3 anchor and s_{a1} < 1/3)
                    for j in Is:
                        cj = []
                        if j - 1 >= a1:
                            cj.append(s_cons(j - 1, -1, Fraction(2)/3))
                        else:
                            pass  # anchor s_{a1-1} = 1 > 2/3 automatic
                        cj.append(s_cons(j, +1, Fraction(1)/3))
                        if frac_lp(cj): res.append((m1case, m2case, 'F3', j))
                    continue
                elif pattern == 'F3any':
                    continue
                if HAS_LP:
                    # also cross-check with scipy (float) — optional
                    pass
                if frac_lp(cons): res.append((m1case, m2case, pattern, None))
    return res


def main():
    infeasible_all = True
    bad = []
    total = 0
    for r in range(1, 21):
        for a1 in range(1, r + 2):
            for a2 in range(a1, r + 2):
                total += 1
                res = sweep_fail_lp(r, a1, a2)
                if res:
                    infeasible_all = False
                    if len(bad) < 8:
                        bad.append((r, a1, a2, res[:4]))
    print("parameter triples checked:", total)
    if infeasible_all:
        print("ALL INFEASIBLE => sweep lemma is PURE ALGEBRA (fail+mu structure excludes F1/F2/F3)")
    else:
        print("FEASIBLE sweep-failure found in", len(bad), "shown regimes (need mu-realizability or rescue):")
        for b in bad: print("   ", b)


if __name__ == '__main__':
    main()
