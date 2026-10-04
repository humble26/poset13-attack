"""测试"跳跃界"路线：
若 D 块阶梯带内每步行/列跳跃 <= 1/3，则单调序列只能以 <=1/3 的步长越过，
不能从 <=1/3 跳到 >=2/3 而不落进 [1/3,2/3]（除非恰为 1/3 -> 2/3，两端点都在区间内）。
=> 若带内存在某行取值 <=1/3 且同行另有 >=2/3，则 9.D 得证。

本脚本测：
  (A) D 块带内最大行跳跃 / 列跳跃，以及 >1/3 的比例
  (B) 是否存在"扫过"的行（带内 min<=1/3 且 max>=2/3）
  (C) 扫过的行中，是否真的落进 [1/3,2/3]（9.D 的实例化）
"""
import sys, collections
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e11_corner.py").resolve(), encoding="utf-8").read().split("def main():")[0])

third = Fraction(1, 3); twothird = Fraction(2, 3)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    st = collections.Counter()
    maxjump = Fraction(0); maxjump_case = None
    sweeps = 0; sweep_ok = 0
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
                # D 块带
                band = [(k, l) for k in range(1, iD+1) for l in range(1, jD+1)
                        if not (L[Q[A[k-1]]][Q[B[l-1]]] or L[Q[B[l-1]]][Q[A[k-1]]])]
                if not band: continue
                st['cases'] += 1
                # 行跳跃
                for k in range(1, iD+1):
                    ls = sorted(l for (kk, l) in band if kk == k)
                    vals = [pval(L, n, Q[A[k-1]], Q[B[l-1]]) for l in ls]
                    for i in range(len(vals)-1):
                        d = vals[i+1]-vals[i]
                        if d > third: st['row_jump_gt'] += 1
                        if d > maxjump: maxjump, maxjump_case = d, (n, z, rel, part, k, ls[i], ls[i+1])
                    # 扫过？
                    if vals and min(vals) <= third and max(vals) >= twothird:
                        sweeps += 1
                        if any(third <= v <= twothird for v in vals): sweep_ok += 1
                # 列跳跃
                for l in range(1, jD+1):
                    ks = sorted(k for (k, ll) in band if ll == l)
                    vals = [pval(L, n, Q[A[k-1]], Q[B[l-1]]) for k in ks]
                    for i in range(len(vals)-1):
                        d = vals[i]-vals[i+1]     # 沿 k 递减
                        if d > third: st['col_jump_gt'] += 1
                        if d > maxjump: maxjump, maxjump_case = d, (n, z, rel, part, 'col', l, ks[i], ks[i+1])
    print(f"n={n} 案例 {st['cases']}")
    print(f"  行跳跃 > 1/3 的次数: {st['row_jump_gt']}")
    print(f"  列跳跃 > 1/3 的次数: {st['col_jump_gt']}")
    print(f"  全局最大跳跃: {maxjump}   case={maxjump_case}")
    print(f"  扫过的行(带内 min<=1/3 且 max>=2/3): {sweeps}  其中落进区间: {sweep_ok}")


if __name__ == '__main__':
    main()
