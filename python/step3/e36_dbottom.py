"""E36: the D-bottom pair claim. In fail cases with iD>=1, jD>=1 and a1||b1:
is p_P(a1,b1) in [1/3,2/3]? Also collect race profiles R(I)=P_I(first=a1)
to test monotonicity along the box, and the capping structure.
"""
import random
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e11_corner.py").resolve(), encoding="utf-8").read().split("def main():")[0])

third = Fraction(1, 3); twothird = Fraction(2, 3)
random.seed(31415)


def races(L, n, Q, A, B, box, iD, jD):
    """R(I(i,j)) = P_I(first = a1) = e(I \\ a1)/e(I) for box cuts with i,j >= 1."""
    out = {}
    for (i, j) in box:
        if i < 1 or j < 1: continue
        Ilist = A[:i] + B[:j]
        pos = {v: k for k, v in enumerate(Ilist)}
        extsI = []
        def rec(placed, seq):
            if len(seq) == len(Ilist): extsI.append(list(seq)); return
            for v in Ilist:
                if (placed >> pos[v]) & 1: continue
                if all((placed >> pos[u]) & 1 for u in Ilist if L[u][v]):
                    rec(placed | 1 << pos[v], seq + [v])
        rec(0, [])
        a1 = Q[A[0]]
        R = Fraction(sum(1 for e in extsI if e[0] == a1), len(extsI))
        out[(i, j)] = R
    return out


def main():
    ok = bad = 0
    mono_ok = mono_bad = 0
    fails = []
    n_case = 0
    for n in (5, 6):
        npairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        if n == 5:
            it = [(mask, z) for mask in range(1 << len(npairs)) for z in range(n)]
        else:
            it = [(random.randrange(1 << len(npairs)), random.randrange(n)) for _ in range(30000)]
        for mask, z in it:
            rel = [npairs[k] for k in range(len(npairs)) if (mask >> k) & 1]
            L = closure(n, rel)
            if not is_poset(L, n): continue
            if antichain_width(L, n) < 2: continue
            Qz = [v for v in range(n) if v != z]
            Ql, _ = induced(L, n, Qz)
            if antichain_width(Ql, n - 1) > 2: continue
            r = analyze(L, n, z)
            if r['zpair']: continue
            for part in chain_partitions(Ql, n - 1):
                cc = corner_check(L, n, z, part)
                if cc is None: continue
                an = cc['an']
                iD, jD = an['iD'], an['jD']
                if iD < 1 or jD < 1: continue
                Q = [v for v in range(n) if v != z]
                A, B = part
                a1, b1 = Q[A[0]], Q[B[0]]
                if L[a1][b1] or L[b1][a1]: continue
                n_case += 1
                p = pval(L, n, a1, b1)
                if third <= p <= twothird: ok += 1
                else:
                    bad += 1
                    if len(fails) < 4:
                        fails.append((n, z, part, str(p), {k: str(v) for k, v in sorted(an['mu'].items())}))
                # race monotonicity along i (fixed j) and reverse along j
                box = list(an['mu'].keys())
                R = races(L, n, Q, A, B, box, iD, jD)
                ks = sorted(R)
                for (i, j) in ks:
                    if (i + 1, j) in R:
                        if R[(i, j)] <= R[(i + 1, j)]: mono_ok += 1
                        else: mono_bad += 1
                    if (i, j + 1) in R:
                        if R[(i, j)] >= R[(i, j + 1)]: mono_ok += 1
                        else: mono_bad += 1
    print("D-bottom pair cases:", n_case, "| balanced:", ok, "unbalanced:", bad)
    for f in fails: print("  UNBAL:", f)
    print("race monotone (R inc in i, dec in j): ok=%d bad=%d" % (mono_ok, mono_bad))


if __name__ == '__main__':
    main()
