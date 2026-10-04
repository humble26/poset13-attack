"""路线 T：D 本身 width-2（若 D×D 有不可比对），Linial 保证 D 内存在平衡对
（p_D ∈ [1/3,2/3]，p_D = 盒最小点 I=D 处的 p_I）。
测试：这些"Linial 见证对"的 p_P 是否也落在 [1/3,2/3]？
若恒真 => 9.D 归约为 Linial on D + 一条"平均稳定性"引理。

同时测几个变体：
 T1 存在 p_D ∈ [1/3,2/3] 且 p_P ∈ [1/3,2/3] 的对
 T2 全部 p_D ∈ [1/3,2/3] 的对里 p_P 的范围
 T3 p_D 最接近 1/2 的对的 p_P 是否在区间
"""
import sys, collections
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e11_corner.py").resolve(), encoding="utf-8").read().split("def main():")[0])

third = Fraction(1, 3); twothird = Fraction(2, 3)


def p_within(less, n, Ilist, x, y):
    Il, idx = induced(less, n, Ilist)
    m = len(Ilist); eI = count_ext(Il, m)
    L2 = closure(m, [(i, j) for i in range(m) for j in range(m) if Il[i][j]] + [(idx[x], idx[y])])
    return Fraction(count_ext(L2, m), eI)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    st = collections.Counter()
    offs = []
    for mask in range(1 << (n * (n - 1) // 2)):
        npairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        rel = [npairs[k] for k in range(len(npairs)) if (mask >> k) & 1]
        L = closure(n, rel)
        if not is_poset(L, n) or antichain_width(L, n) < 2: continue
        for z in range(n):
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
                if len(an['mu']) < 2 or iD < 1 or jD < 1: continue
                A, B = part; Q = [v for v in range(n) if v != z]
                band = [(k, l) for k in range(1, iD+1) for l in range(1, jD+1)
                        if not (L[Q[A[k-1]]][Q[B[l-1]]] or L[Q[B[l-1]]][Q[A[k-1]]])]
                if not band: continue
                st['cases'] += 1
                Imin = [Q[a] for a in A[:iD]] + [Q[b] for b in B[:jD]]
                rows = []
                for (k, l) in band:
                    x, y = Q[A[k-1]], Q[B[l-1]]
                    pd = p_within(L, n, Imin, x, y)
                    pp = pval(L, n, x, y)
                    rows.append((k, l, pd, pp))
                # D 内有 Linial 见证？
                wit = [(k, l, pd, pp) for (k, l, pd, pp) in rows if third <= pd <= twothird]
                if not wit:
                    st['no_D_witness'] += 1
                    continue
                st['has_D_witness'] += 1
                if any(third <= pp <= twothird for (_, _, _, pp) in wit):
                    st['T1_ok'] += 1
                else:
                    st['T1_BAD'] += 1
                    if len(offs) < 6:
                        offs.append((n, z, rel, part, [(k, l, str(pd), str(pp)) for (k, l, pd, pp) in rows]))
                # T3: p_D 最接近 1/2 的对
                best = min(rows, key=lambda t: abs(t[2] - Fraction(1, 2)))
                if third <= best[3] <= twothird: st['T3_ok'] += 1
                else: st['T3_bad'] += 1
    print(f"n={n} 案例 {st['cases']}")
    print(f"  有 D-Linial 见证: {st['has_D_witness']}  无: {st['no_D_witness']}")
    print(f"  T1 见证对的 p_P 也平衡: {st['T1_ok']}   不平衡: {st['T1_BAD']}")
    print(f"  T3 (p_D 最接近1/2 的对) p_P 平衡: {st['T3_ok']}  不平衡: {st['T3_bad']}")
    for o in offs: print("  OFF:", o)


if __name__ == '__main__':
    main()
