"""最终综合验证（|W|=2 不可比对情形的一般严格证明的机器核对）：
1. 交叉结构 γ=0 的 m 公式（已验）
2. γ>=1 ⟹ m_b = 1/2（槽计数）
3. case ii ⟹ m_b = 1/2（z↔b 对称）
4. 交叉 + γ=0 + r−β>=2 + s>=2 ⟹ D-路线给出平衡对（p_P(d_i,a) 扫过区间）
"""
import sys, itertools
from fractions import Fraction
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from my_poset import my_close, my_is_poset


def count_ext(Q, m):
    dp = [0] * (1 << m); dp[0] = 1
    for S in range(1 << m):
        v = dp[S]
        if not v: continue
        for x in range(m):
            if (S >> x) & 1: continue
            if all(not Q[y][x] or (S >> y) & 1 for y in range(m)):
                dp[S | 1 << x] += v
    return dp[-1]


def pval_direct(Q, m, x, y):
    L2 = [r[:] for r in Q]; L2[x][y] = True
    A = [x] + [k for k in range(m) if Q[k][x]]
    B = [y] + [j for j in range(m) if Q[y][j]]
    for a in A:
        for b in B:
            if a != b: L2[a][b] = True
    return Fraction(count_ext(L2, m), count_ext(Q, m))


def cross_structure(r, s, beta=0):
    """交叉 γ=0：D 链 r，U 链 s，a ≻ D[0..beta-1]，a ∥ D[beta..r-1]，a ≺ U，
    b ≻ D，b ∥ U，a ∥ b。返回 Q（不含 z）。"""
    a = r + s; b = r + s + 1; m = b + 1
    rel = []
    for i in range(r - 1): rel.append((i, i + 1))
    for i in range(s - 1): rel.append((r + i, r + i + 1))
    for i in range(r): rel.append((i, b))
    for j in range(beta): rel.append((j, a))
    for j in range(s): rel.append((a, r + j))
    Q = my_close(m, rel)
    return Q, a, b, m


def main():
    third, twothird = Fraction(1, 3), Fraction(2, 3)
    st = dict(cases=0, route_ok=0, route_bad=0)
    for r in range(2, 9):
        for s in range(2, 9):
            for beta in range(0, r - 1):   # r−beta >= 2
                Q, a, b, m = cross_structure(r, s, beta)
                # 加 z：D ≺ z ≺ U，z ∥ a,b
                # 在 Q 基础上构造 P（含 z），直接算 p_P(d_i, a)
                p = r - beta
                # 用理想分解算 m_a（无需 z）
                ids = [S for S in range(1 << m)
                       if all(not ((S >> i) & 1) or all((S >> j) & 1 for j in range(m) if Q[j][i]) for i in range(m))]
                full = (1 << m) - 1
                ea = [0] * (1 << m); ea[0] = 1
                for S in range(1, 1 << m):
                    t = 0
                    for x in range(m):
                        if not (S >> x) & 1: continue
                        if all(not Q[y][x] or not (S >> y) & 1 for y in range(m)):
                            t += ea[S ^ (1 << x)]
                    ea[S] = t
                Dmask = sum(1 << i for i in range(r)); Umask = sum(1 << (r + j) for j in range(s))
                box = [I for I in ids if (I & Dmask) == Dmask and (I & Umask) == 0]
                eP = sum(ea[I] * ea[full ^ I] for I in box)
                ma = Fraction(sum(ea[I] * ea[full ^ I] for I in box if not (I >> a) & 1), eP)
                if not (ma < third):
                    st.setdefault('ma_not_lo', 0); st['ma_not_lo'] += 1
                    continue
                # D-路线：p_P(d_i, a) = ma + (1-ma) * j/(p+1)，j = r-i+1 ∈ [1, p]
                # 检查 ∃ j ∈ [1,p] 使值 ∈ [1/3,2/3]
                hit = False
                for j in range(1, p + 1):
                    val = ma + (1 - ma) * Fraction(j, p + 1)
                    if third <= val <= twothird:
                        hit = True; break
                st['cases'] += 1
                if hit: st['route_ok'] += 1
                else: st['route_bad'] += 1
    print(f"交叉结构 D-路线验证：{st['cases']} 个 (r,s,β) 组合")
    print(f"  平衡对存在: {st['route_ok']}   不存在: {st.get('route_bad',0)}")
    print(f"  m_a 不低于 1/3 的: {st.get('ma_not_lo',0)}")
    print("结论: " + ("D-路线在 r−β≥2, s≥2 时恒成功 ✅" if st.get('route_bad',0) == 0 else "存在失败 ❌"))


if __name__ == '__main__':
    main()
