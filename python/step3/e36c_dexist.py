"""E36c: the D-route EXISTENCE claim: fail + multi-cut box + DxD has an
incomparable pair => SOME DxD pair is balanced. Dual for UxU. Combined coverage.
"""
import random
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e11_corner.py").resolve(), encoding="utf-8").read().split("def main():")[0])

third = Fraction(1, 3); twothird = Fraction(2, 3)
random.seed(2718)


def main():
    stats = dict(cases=0, d_avail=0, d_ok=0, d_bad=0, u_avail=0, u_ok=0, u_bad=0,
                 both_missing=0, covered=0, uncovered=0)
    fails = []
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
                iD, jD = an['iD'], an['jD']
                iU, jU = an['iU'], an['jU']
                if len(an['mu']) < 2: continue  # multi-cut (single-cut closed)
                s_, t_ = an['s'], an['t']
                Q = [v for v in range(n) if v != z]
                A, B = part
                stats['cases'] += 1
                # DxD pairs
                dpairs = []
                if iD >= 1 and jD >= 1:
                    for k in range(1, iD + 1):
                        for l in range(1, jD + 1):
                            x, y = Q[A[k - 1]], Q[B[l - 1]]
                            if not (L[x][y] or L[y][x]):
                                dpairs.append(pval(L, n, x, y))
                upairs = []
                if s_ - iU >= 1 and t_ - jU >= 1:
                    for k in range(iU + 1, s_ + 1):
                        for l in range(jU + 1, t_ + 1):
                            x, y = Q[A[k - 1]], Q[B[l - 1]]
                            if not (L[x][y] or L[y][x]):
                                upairs.append(pval(L, n, x, y))
                d_avail = len(dpairs) > 0
                u_avail = len(upairs) > 0
                if d_avail:
                    stats['d_avail'] += 1
                    if any(third <= p_ <= twothird for p_ in dpairs): stats['d_ok'] += 1
                    else:
                        stats['d_bad'] += 1
                        if len(fails) < 4:
                            fails.append(('D', n, z, part, [str(p_) for p_ in dpairs]))
                if u_avail:
                    stats['u_avail'] += 1
                    if any(third <= p_ <= twothird for p_ in upairs): stats['u_ok'] += 1
                    else:
                        stats['u_bad'] += 1
                        if len(fails) < 8:
                            fails.append(('U', n, z, part, [str(p_) for p_ in upairs]))
                if not d_avail and not u_avail:
                    stats['both_missing'] += 1
                if (d_avail and any(third <= p_ <= twothird for p_ in dpairs)) or \
                   (u_avail and any(third <= p_ <= twothird for p_ in upairs)):
                    stats['covered'] += 1
                else:
                    stats['uncovered'] += 1
    print(stats)
    for f in fails: print("  FAIL:", f)


if __name__ == '__main__':
    main()
