"""测 p_P 相对 p_D 的位移结构（D×D 带内，fail 案例）：
 S1  p_P(x,y) >= p_D(x,y) 恒成立？        （单向位移猜想）
 S2  p_P - p_D <= 1/3 恒成立？             （位移上界）
 S3  min/max 位移幅度分布；存活对的 p_D 位置分布
"""
import sys, collections
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e11_corner.py").resolve(), encoding="utf-8").read().split("def main():")[0])

third = Fraction(1, 3); twothird = Fraction(2, 3)


def p_within(less, n, Ilist, x, y):
    Il, idx = induced(less, n, Ilist); m = len(Ilist); eI = count_ext(Il, m)
    L2 = closure(m, [(i, j) for i in range(m) for j in range(m) if Il[i][j]] + [(idx[x], idx[y])])
    return Fraction(count_ext(L2, m), eI)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    st = collections.Counter()
    maxshift = Fraction(-9); minshift = Fraction(9)
    for mask in range(1 << (n * (n - 1) // 2)):
        npairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        rel = [npairs[k] for k in range(len(npairs)) if (mask >> k) & 1]
        L = closure(n, rel)
        if not is_poset(L, n) or antichain_width(L, n) < 2: continue
        for z in range(n):
            Qz = [v for v in range(n) if v != z]; Ql, _ = induced(L, n, Qz)
            if antichain_width(Ql, n - 1) > 2: continue
            r = analyze(L, n, z)
            if r['zpair']: continue
            for part in chain_partitions(Ql, n - 1):
                cc = corner_check(L, n, z, part)
                if cc is None: continue
                an = cc['an']; iD, jD = an['iD'], an['jD']
                if len(an['mu']) < 2 or iD < 1 or jD < 1: continue
                A, B = part; Q = [v for v in range(n) if v != z]
                band = [(k, l) for k in range(1, iD+1) for l in range(1, jD+1)
                        if not (L[Q[A[k-1]]][Q[B[l-1]]] or L[Q[B[l-1]]][Q[A[k-1]]])]
                if not band: continue
                Imin = [Q[a] for a in A[:iD]] + [Q[b] for b in B[:jD]]
                for (k, l) in band:
                    x, y = Q[A[k-1]], Q[B[l-1]]
                    pd = p_within(L, n, Imin, x, y)
                    pp = pval(L, n, x, y)
                    d = pp - pd
                    st['pairs'] += 1
                    if d >= 0: st['S1_ge'] += 1
                    else: st['S1_lt'] += 1
                    if d > third: st['S2_gt'] += 1
                    if d > maxshift: maxshift = d
                    if d < minshift: minshift = d
    print(f"n={n} D×D 带内对数 {st['pairs']}")
    print(f"  S1  p_P >= p_D : {st['S1_ge']}    p_P < p_D : {st['S1_lt']}")
    print(f"  S2  p_P-p_D > 1/3 : {st['S2_gt']}")
    print(f"  位移范围: [{minshift}, {maxshift}]")


if __name__ == '__main__':
    main()
