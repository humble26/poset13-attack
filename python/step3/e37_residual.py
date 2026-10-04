"""E37: in residual fail cases (multi-cut, D and U both chain-like),
classify where the balanced Q-pairs live: (D,W), (W,W), (W,U), (D,U), etc.
"""
import random, collections
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e11_corner.py").resolve(), encoding="utf-8").read().split("def main():")[0])

third = Fraction(1, 3); twothird = Fraction(2, 3)
random.seed(9999)


def chainlike(L, n, Q, A, B, lo, hi):
    """cross pairs between A[loA..hiA-1] and B[loB..hiB-1] all comparable?
    Here: D-chainlike = all a in A[:iD] x b in B[:jD] comparable; U similarly."""
    return True  # placeholder, real check inline


def main():
    hist = collections.Counter()
    n_res = 0
    nocross = 0
    for n in (5, 6):
        npairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        if n == 5:
            it = [(mask, z) for mask in range(1 << len(npairs)) for z in range(n)]
        else:
            it = [(random.randrange(1 << len(npairs)), random.randrange(n)) for _ in range(25000)]
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
                iD, jD, iU, jU = an['iD'], an['jD'], an['iU'], an['jU']
                if len(an['mu']) < 2: continue
                s_, t_ = an['s'], an['t']
                Q = [v for v in range(n) if v != z]
                A, B = part
                # DxD incomparable exists?
                dinc = any(not (L[Q[A[k-1]]][Q[B[l-1]]] or L[Q[B[l-1]]][Q[A[k-1]]])
                           for k in range(1, iD + 1) for l in range(1, jD + 1)) if (iD >= 1 and jD >= 1) else False
                uinc = any(not (L[Q[A[k-1]]][Q[B[l-1]]] or L[Q[B[l-1]]][Q[A[k-1]]])
                           for k in range(iU + 1, s_ + 1) for l in range(jU + 1, t_ + 1)) if (s_ - iU >= 1 and t_ - jU >= 1) else False
                if dinc or uinc: continue  # covered by 9.8
                n_res += 1
                # find balanced pairs and classify
                cats = set()
                for k in range(1, s_ + 1):
                    for l in range(1, t_ + 1):
                        x, y = Q[A[k - 1]], Q[B[l - 1]]
                        if L[x][y] or L[y][x]: continue
                        p = pval(L, n, x, y)
                        if not (third <= p <= twothird): continue
                        ca = 'D' if k <= iD else ('U' if k > iU else 'W')
                        cb = 'D' if l <= jD else ('U' if l > jU else 'W')
                        cats.add(ca + cb)
                if not cats:
                    nocross += 1
                    continue
                for c_ in cats: hist[c_] += 1
    print("residual (chain-like D,U) instances:", n_res, "| no balanced pair:", nocross)
    print("balanced-pair categories:", dict(hist))


if __name__ == '__main__':
    main()
