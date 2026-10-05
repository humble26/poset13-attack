"""E44: full parameter-space sweep test for the chain-box case (M1).
Enumerate (r, a1, a2, sp, b1, b2): build P, compute exact sweep values
s_i = p_P(d_i, w1) for i in [a1, r], classify: landed (some s_i in [1,3,2/3])
vs failure patterns F1 (all < 1/3) / F2 (all > 2/3) / F3 (jump over).
Also cross-check the corrected closed forms U0, e2, cnt(i) per instance.
"""
from fractions import Fraction
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent))
from my_poset import my_close, my_is_poset, my_induced, my_count_ext

third = Fraction(1, 3)
twothird = Fraction(2, 3)


def build(r, a1, a2, sp, b1, b2):
    # labels: d_i = i-1 (i=1..r), w1=r, w2=r+1, u_j = r+1+j (j=1..sp), z = n-1
    n = r + 2 + sp + 1
    z = n - 1
    w1, w2 = r, r + 1
    u = list(range(r + 2, z))
    rel = [(i, i + 1) for i in range(r - 1)]
    rel.append((w1, w2))
    rel += [(u[k], u[k + 1]) for k in range(sp - 1)]
    for j in range(a1 - 1): rel.append((j, w1))
    for j in range(a2 - 1): rel.append((j, w2))
    for j in range(b1, sp): rel.append((w1, u[j]))   # w1 < u_j (0-based j >= b1)
    for j in range(b2, sp): rel.append((w2, u[j]))
    for d in range(r): rel.append((d, z))            # d < z
    for uu in u: rel.append((z, uu))                 # z < u
    L = my_close(n, rel)
    return L, n, w1, w2, u, z


def main():
    st = dict(sweep=0, landed=0, F1=0, F2=0, F3=0)
    fails = []
    closed_ok = closed_bad = 0
    for r in range(1, 7):
        for a1 in range(1, r + 2):
            for a2 in range(a1, r + 2):
                for sp in range(0, 5):
                    for b1 in range(0, sp + 1):
                        for b2 in range(b1, sp + 1):
                            L, n, w1, w2, u, z = build(r, a1, a2, sp, b1, b2)
                            if not my_is_poset(L, n): continue
                            Dset = list(range(r))
                            I1 = Dset + [w1]
                            I2 = I1 + [w2]
                            Uset = u
                            # corrected closed forms
                            U0 = (b1 + 1) * (b2 + 1) - b1 * (b1 + 1) // 2
                            N1 = r - a1 + 2
                            N2 = r - a2 + 2
                            e2 = (a2 - a1) * N2 + N2 * (N2 + 1) // 2
                            eP = U0 + N1 * (b2 + 1) + e2
                            # exact e(P) check
                            eP_direct = my_count_ext(L, n)
                            if eP != eP_direct:
                                closed_bad += 1
                                if closed_bad <= 3: fails.append(("eP", r, a1, a2, sp, b1, b2, eP, eP_direct))
                            else: closed_ok += 1
                            # sweep values: p_P(d_i, w1) via count with added relation
                            vals = []
                            for i in range(a1, r + 1):
                                d = i - 1
                                L2 = my_close(n, [(a, b) for a in range(n) for b in range(n) if L[a][b]] + [(d, w1)])
                                s_i = Fraction(my_count_ext(L2, n), eP_direct)
                                # closed form check
                                cnt = None
                                if i >= a2: cnt = (r - i + 1) * (r - i + 2) // 2
                                else: cnt = (a2 - 1 - i) * N2 + N2 * (N2 + 1) // 2
                                s_closed = Fraction(U0 + (b2 + 1) * (r - i + 1) + cnt, eP)
                                if s_i != s_closed:
                                    closed_bad += 1
                                    if closed_bad <= 6: fails.append(("s_i", r, a1, a2, sp, b1, b2, i, str(s_i), str(s_closed)))
                                vals.append(s_i)
                            if not vals: continue
                            st['sweep'] += 1
                            if any(third <= v <= twothird for v in vals): st['landed'] += 1
                            elif all(v < third for v in vals): st['F1'] += 1
                            elif all(v > twothird for v in vals): st['F2'] += 1
                            else:
                                st['F3'] += 1
                                if len(fails) < 8:
                                    fails.append(("F3", r, a1, a2, sp, b1, b2, [str(v) for v in vals]))
    print(st)
    print("closed-form checks: ok=%d bad=%d" % (closed_ok, closed_bad))
    for f in fails: print("  ", f)


if __name__ == '__main__':
    main()
