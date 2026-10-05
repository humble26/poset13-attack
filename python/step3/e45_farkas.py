"""E45: extract exact Farkas certificates for the chain-box sweep lemma.
For each (r, a1, a2) and each fail/pattern combo: the constraint set (in variables
m1, mu1; mu2 = 1 - m1 - mu1) is:
  (fail)  m1 </> 1/3 ;  m1+mu1 </> 1/3
  (sweep) m1 + mu1*A(i) + mu2*B(i) </> 1/3, 2/3   for each sweep i
  (pos)   m1, mu1, mu2 >= 0
Infeasibility certificate: nonneg lambda over constraints with
sum(lambda*a) = 0, sum(lambda*b) = 0, sum(lambda*c) < 0.
Enumerate support <= 3 for exact certificates; report symbolic structure.
"""
import itertools
from fractions import Fraction


def A_of(i, r, a1):
    return Fraction(r - i + 1, r - a1 + 2)


def B_profile(r, a1, a2):
    e2 = (a2 - a1) * (r - a2 + 2) + (r - a2 + 2) * (r - a2 + 3) // 2
    out = []
    for i in range(a1, r + 1):
        if i >= a2:
            cnt = (r - i + 1) * (r - i + 2) // 2
        else:
            cnt = (a2 - 1 - i) * (r - a2 + 2) + (r - a2 + 2) * (r - a2 + 3) // 2
        out.append((i, Fraction(cnt, e2)))
    return out, e2


def s_coeff(i, r, a1, a2):
    """s_i = c0 + c1*m1 + c2*mu1 (mu2 = 1-m1-mu1). Returns (c0, c1, c2).
    The constraint 's_i <= t' is c1*m1 + c2*mu1 <= t - c0."""
    B = dict(B_profile(r, a1, a2)[0])[i]
    A = A_of(i, r, a1)
    # s_i = B + m1*(1-B) + mu1*(A-B)
    return (B, 1 - B, A - B)


def build_constraints(r, a1, a2, m1case, m2case, pattern, j=None):
    """returns list of (a1, a2, c, name): a1*m1 + a2*mu1 <= c"""
    cons = []
    if m1case == 'lo': cons.append((Fraction(1), Fraction(0), Fraction(1, 3), 'm1<1/3'))
    else: cons.append((Fraction(-1), Fraction(0), Fraction(-2, 3), 'm1>2/3'))
    if m2case == 'lo': cons.append((Fraction(1), Fraction(1), Fraction(1, 3), 'm2<1/3'))
    else: cons.append((Fraction(-1), Fraction(-1), Fraction(-2, 3), 'm2>2/3'))
    Bp, e2 = B_profile(r, a1, a2)
    for (i, B) in Bp:
        (c0, c1, c2) = s_coeff(i, r, a1, a2)
        if pattern == 'F1':
            cons.append((c1, c2, Fraction(1, 3) - c0, f's{i}<1/3'))
        elif pattern == 'F2':
            cons.append((-c1, -c2, c0 - Fraction(2, 3), f's{i}>2/3'))
        elif pattern == 'F3':
            if i < j:
                cons.append((-c1, -c2, c0 - Fraction(2, 3), f's{i}>2/3'))
            elif i == j:
                cons.append((c1, c2, Fraction(1, 3) - c0, f's{i}<1/3'))
            # i > j: unconstrained in F3's jump definition (between: any)
            # hmm — actually F3's jump: s_{j-1} > 2/3, s_j < 1/3 — the rows
            # between j-1 and j don't exist (adjacent indices) — but the sweep
            # index steps by 1: the jump is between consecutive i's ✓
    # positivity: m1 >= 0, mu1 >= 0, mu2 >= 0
    cons.append((Fraction(-1), Fraction(0), Fraction(0), 'm1>=0'))
    cons.append((Fraction(0), Fraction(-1), Fraction(0), 'mu1>=0'))
    cons.append((Fraction(1), Fraction(1), Fraction(1), 'mu2>=0'))
    return cons


def farkas_cert(cons):
    """find nonneg lambda over cons: sum(lam*a)=0, sum(lam*b)=0, sum(lam*c)<0.
    Try supports of size 2 and 3. Returns (support, lambdas) or None."""
    n = len(cons)
    # size 2: needs a-vectors opposite
    for i, j in itertools.combinations(range(n), 2):
        (a1, a2, c, n1_), (b1, b2, cc, n2_) = cons[i], cons[j]
        # lam_i*(a1,a2) + lam_j*(b1,b2) = 0
        if (a1, a2) == (0, 0):
            if c < 0: return [(i, Fraction(1)), (j, Fraction(0))]
            continue
        if (b1, b2) == (0, 0): continue
        # lam_i/lam_j = -(b1,b2)/(a1,a2): need proportional opposite
        if a1 * b2 - a2 * b1 != 0: continue
        if a1 == 0 or b1 == 0:
            if a1 == 0 and b1 == 0: continue
        # ratio: lam_i*a1 + lam_j*b1 = 0 and a2,lam: lam_i*a2 = -lam_j*b2
        if a1 * b2 == a2 * b1 and a1 * b2 <= 0 and (a1, a2) != (0, 0):
            # opposite direction
            if a1 != 0:
                lam_i = Fraction(abs(b1), abs(a1)) if a1 * b1 <= 0 else None
            elif a2 != 0:
                lam_i = Fraction(abs(b2), abs(a2)) if a2 * b2 <= 0 else None
            if lam_i is None: continue
            lam_j = Fraction(1)
            if a1 != 0: lam_i = Fraction(-b1 * lam_j, a1)
            else: lam_i = Fraction(-b2 * lam_j, a2)
            if lam_i < 0: lam_i, lam_j = -lam_i, -lam_j
            if lam_i < 0 or lam_j < 0: continue
            tot_c = lam_i * c + lam_j * cc
            if tot_c < 0: return [(i, lam_i), (j, lam_j)]
    # size 3: enumerate triples with rational solve
    for tri in itertools.combinations(range(n), 3):
        (v1, v2, v3) = [cons[k][:2] for k in tri]
        (c1, c2, c3) = [cons[k][2] for k in tri]
        # solve lam1*v1 + lam2*v2 + lam3*v3 = 0 with lam >= 0 (1-dim solution space)
        # cross product of v1, v2 gives the normal direction of their span:
        det = v1[0] * v2[1] - v1[1] * v2[0]
        if det == 0: continue
        # v3 must be in span(v1, v2): always true in 2D
        # solve lam1*v1 + lam2*v2 = -lam3*v3: set lam3 = 1:
        # lam1*v1 + lam2*v2 = -v3
        # 2x2 solve:
        (a1, a2) = v1; (b1, b2) = v2; (t1, t2) = (-v3[0], -v3[1])
        det = a1 * b2 - a2 * b1
        if det == 0: continue
        lam1 = Fraction(t1 * b2 - a2 * t2, det)
        lam2 = Fraction(a1 * t2 - t1 * b1, det)
        lam3 = Fraction(1)
        if lam1 < 0 or lam2 < 0: continue
        tot_c = lam1 * c1 + lam2 * c2 + lam3 * c3
        if tot_c < 0: return [(tri[0], lam1), (tri[1], lam2), (tri[2], lam3)]
    return None


def main():
    n_param = 0
    cert_count = 0
    no_cert = []
    for r in range(1, 9):
        for a1 in range(1, r + 2):
            for a2 in range(a1, r + 2):
                n_param += 1
                Bp, e2 = B_profile(r, a1, a2)
                if not Bp: continue  # empty sweep: vacuous
                for m1case in ('lo', 'hi'):
                    for m2case in ('lo', 'hi'):
                        if m1case == 'hi' and m2case == 'lo': continue
                        for pattern in ('F1', 'F2', 'F3all'):
                            if pattern == 'F3all':
                                # F3: jump at each possible j (separately)
                                for j in range(a1, r + 1):
                                    cons = build_constraints(r, a1, a2, m1case, m2case, 'F3', j)
                                    cert = farkas_cert(cons)
                                    if cert is None:
                                        print(f"NO-CERT r={r} a1={a1} a2={a2} {m1case}/{m2case} F3 j={j}")
                                        no_cert.append((r, a1, a2, m1case, m2case, 'F3', j))
                                    else: cert_count += 1
                            else:
                                cons = build_constraints(r, a1, a2, m1case, m2case, pattern)
                                cert = farkas_cert(cons)
                                if cert is None:
                                    print(f"NO-CERT r={r} a1={a1} a2={a2} {m1case}/{m2case} {pattern}")
                                    no_cert.append((r, a1, a2, m1case, m2case, pattern))
                                else: cert_count += 1
    print(f"params: {n_param}, certificates found: {cert_count}, missing: {len(no_cert)}")
    for nc in no_cert[:6]: print("   MISSING:", nc)


if __name__ == '__main__':
    main()
