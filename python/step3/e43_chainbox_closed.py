"""E43 (clean): verify the chain-box closed forms
  e(Q∖D) = beta2-beta1+1,  e(Q∖I1) = beta2+1,  e(Q∖I2) = 1,
  m1 = W/e(P) <= 1/2,
by building the full poset P from (r, a1, a2, sp, b1, b2) and counting induced sub-posets.
"""
from fractions import Fraction
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from my_poset import my_close, my_is_poset, my_induced, my_count_ext


def build(r, a1, a2, sp, b1, b2):
    """returns (L, n, labels). Labels: d_i = i-1 (i=1..r), w1=r, w2=r+1, u_j = r+1+j (j=1..sp)."""
    n = r + 2 + sp
    w1, w2 = r, r + 1
    u = list(range(r + 2, n))
    rel = [(i, i + 1) for i in range(r - 1)]          # d_i < d_{i+1}
    rel.append((w1, w2))                               # w1 < w2
    rel += [(u[k], u[k + 1]) for k in range(sp - 1)]   # U chain
    for j in range(a1 - 1): rel.append((j, w1))        # d_j < w1 (j < a1-1... 0-based: j < a1-1)
    for j in range(a2 - 1): rel.append((j, w2))
    for j in range(b2, sp): rel.append((w1, u[j]))     # w1 < u_j (0-based j >= b2)
    for j in range(b2, sp): rel.append((w2, u[j]))
    L = my_close(n, rel)
    return L, n, w1, w2, u


def main():
    ok = bad = 0
    fails = []
    m1_ok = m1_bad = 0
    for r in range(1, 7):
        for a1 in range(1, r + 2):
            for a2 in range(a1, r + 2):
                for sp in range(0, 5):
                    for b1 in range(0, sp + 1):
                        for b2 in range(b1, sp + 1):
                            L, n, w1, w2, u = build(r, a1, a2, sp, b1, b2)
                            if not my_is_poset(L, n): continue
                            Dset = list(range(r))
                            I1 = Dset + [w1]
                            I2 = I1 + [w2]
                            Wf = b2 - b1 + 1
                            def e_of(subset):
                                sub = sorted(subset)
                                Ls = my_induced(L, sub)
                                return my_count_ext(Ls, len(sub))
                            eQD = e_of([v for v in range(n) if v not in set(Dset)])
                            eQI1 = e_of([v for v in range(n) if v not in set(I1)])
                            eQI2 = e_of([v for v in range(n) if v not in set(I2)])
                            checks = [
                                ("e(QD)", eQD, Wf),
                                ("e(QI1)", eQI1, b2 + 1),
                                ("e(QI2)", eQI2, 1),
                            ]
                            for (nm, got, want) in checks:
                                if got == want: ok += 1
                                else:
                                    bad += 1
                                    if len(fails) < 5: fails.append((nm, r, a1, a2, sp, b1, b2, got, want))
                            # m1 = e(QD)*1 / e(P); check <= 1/2
                            eP = my_count_ext(L, n)
                            m1 = Fraction(eQD, eP)
                            if m1 <= Fraction(1, 2): m1_ok += 1
                            else:
                                m1_bad += 1
                                if len(fails) < 8: fails.append(("m1>1/2", r, a1, a2, sp, b1, b2, str(m1)))
    print("closed-form checks: ok=%d bad=%d | m1<=1/2: ok=%d bad=%d" % (ok, bad, m1_ok, m1_bad))
    for f in fails: print("  FAIL:", f)


if __name__ == '__main__':
    main()
