"""测试 9.D 的候选证明路线（n=6 穷举 fail+multicut+DxD 案例）：

J1  经典单调性：l -> p_P(a_k,b_l) 单调不减；k -> p_P(a_k,b_l) 单调不增
J2  p_I(x,y) 沿盒子（积序）单调？（若真，则 p_P ∈ [min_I, max_I]）
J3  存在 DxD 不可比对对，其 p_I ∈ [1/3,2/3] 对**所有** I ∈ box 成立
J4  该对只依赖 k,l 的"中位"位置？
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
    m = len(Ilist)
    eI = count_ext(Il, m)
    L2 = closure(m, [(i, j) for i in range(m) for j in range(m) if Il[i][j]] + [(idx[x], idx[y])])
    return Fraction(count_ext(L2, m), eI)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    st = collections.Counter()
    j2_viol = []
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
                A, B = part; s_, t_ = an['s'], an['t']
                Q = [v for v in range(n) if v != z]
                dp = [(k, l) for k in range(1, iD+1) for l in range(1, jD+1)
                      if not (L[Q[A[k-1]]][Q[B[l-1]]] or L[Q[B[l-1]]][Q[A[k-1]]])]
                if not dp: continue
                st['cases'] += 1
                # J1: l -> p_P(a_k,b_l) 单调不减 (在 D 块内可比/不可比都看)
                for k in range(1, iD+1):
                    prev = None
                    for l in range(1, jD+2):
                        if l <= jD:
                            p = pval(L, n, Q[A[k-1]], Q[B[l-1]])
                        else:
                            p = Fraction(1)  # 哨兵：b 超出 D
                        if prev is not None and p < prev:
                            st['J1_viol'] += 1
                        prev = p
                # J2: p_I(x,y) 沿盒子单调（对每个 DxD 对检查）
                for (k, l) in dp:
                    x, y = Q[A[k-1]], Q[B[l-1]]
                    vals = {}
                    for (i, j), w in an['mu'].items():
                        Ilist = [Q[a] for a in A[:i]] + [Q[b] for b in B[:j]]
                        vals[(i, j)] = p_within(L, n, Ilist, x, y)
                    # 单调：i,j 增大 => 值非减？（取所有可比对检查）
                    keys = sorted(vals)
                    mono_up = True
                    for (i1, j1) in keys:
                        for (i2, j2) in keys:
                            if i1 <= i2 and j1 <= j2 and vals[(i1, j1)] > vals[(i2, j2)]:
                                mono_up = False
                    if mono_up: st['J2_mono_up'] += 1
                    else: st['J2_not'] += 1
                    # J3: 该对全 I 平衡？
                    if all(third <= v <= twothird for v in vals.values()):
                        st['J3_this_ok'] += 1
                # J3 案例级：是否存在这样的对
                case_j3 = False
                for (k, l) in dp:
                    x, y = Q[A[k-1]], Q[B[l-1]]
                    vals = {}
                    for (i, j), w in an['mu'].items():
                        Ilist = [Q[a] for a in A[:i]] + [Q[b] for b in B[:j]]
                        vals[(i, j)] = p_within(L, n, Ilist, x, y)
                    if vals and all(third <= v <= twothird for v in vals.values()):
                        case_j3 = True
                if case_j3: st['J3_case_ok'] += 1
                else: st['J3_case_no'] += 1
    print(f"n={n} 9.D 案例数: {st['cases']}")
    print(f"  J1 违反(l-单调): {st['J1_viol']}")
    print(f"  J2 p_I 沿盒子单调: {st['J2_mono_up']}  非单调: {st['J2_not']}")
    print(f"  J3 单对全 I 平衡: {st['J3_this_ok']}")
    print(f"  J3 案例级存在: {st['J3_case_ok']}  不存在: {st['J3_case_no']}")


if __name__ == '__main__':
    main()
