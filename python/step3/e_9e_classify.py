"""结构分类：|W|=2 不可比对、D/U 链、width≤2 的全部可能结构，
对每个结构计算 m_a, m_b，判定哪些满足 fail（都 ∉ [1/3,2/3]）。

结构参数化（W={a,b}，a∥b）：
  对 d∈D：d≺a 或 d∥a；d≺b 或 d∥b；不同时∥（width）。
  对 u∈U：a≺u 或 a∥u；b≺u 或 b∥u；不同时∥。
  {d≺a},{d≺b} 是 D 前缀；{a≺u},{b≺u} 是 U 后缀。
  前缀覆盖 D ⟹ 其一是全 D；后缀覆盖 U ⟹ 其一是全 U。WLOG 取 b≻全D。

四种组合（b≻D 固定后）：
  A1: a≺全U（交叉），带 a 的 D-前缀 β、b 的 U-后缀 γ
  A2: b≺全U（b 完全夹中间），带 a 的 D-前缀 β、a 的 U-后缀 γ
再加 D×U 交叉关系（可选）。
"""
import sys, itertools, collections
from fractions import Fraction
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from my_poset import my_close, my_is_poset


def chain_sort(Q, els):
    return sorted(els, key=lambda v: sum(1 for u in els if Q[u][v]))


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


def mvalues(r, s, cross=True, beta=0, gamma=0):
    """构造结构并返回 (m_a, m_b)。cross: a≺U；否则 b≺U。
    beta = a 的 D-前缀长（a ≻ D[0..beta-1]）；gamma = 后缀参数。"""
    # 元素：D=0..r-1, U=r..r+s-1, a=r+s, b=r+s+1
    a = r + s; b = r + s + 1
    m = b + 1
    rel = []
    for i in range(r - 1): rel.append((i, i + 1))
    for i in range(s - 1): rel.append((r + i, r + i + 1))
    for i in range(r): rel.append((i, b))            # b ≻ D
    for j in range(beta): rel.append((j, a))         # a ≻ D[0..beta-1]
    if cross:
        for j in range(s): rel.append((a, r + j))     # a ≺ U
        for j in range(gamma): rel.append((b, r + j)) # b ≺ U[0..gamma-1]
    else:
        for j in range(s): rel.append((b, r + j))     # b ≺ U
        for j in range(gamma, s): rel.append((a, r + j))  # a ≺ U[gamma..s-1]
    Q = my_close(m, rel)
    if not my_is_poset(Q, m): return None
    # 宽度检查
    for x, y, z in itertools.combinations(range(m), 3):
        if not (Q[x][y] or Q[y][x] or Q[x][z] or Q[z][x] or Q[y][z] or Q[z][y]):
            return None
    ids = [S for S in range(1 << m)
           if all(not ((S >> i) & 1) or all((S >> j) & 1 for j in range(m) if Q[j][i]) for i in range(m))]
    full = (1 << m) - 1
    e_all = [0] * (1 << m); e_all[0] = 1
    for S in range(1, 1 << m):
        t = 0
        for x in range(m):
            if not (S >> x) & 1: continue
            if all(not Q[y][x] or not (S >> y) & 1 for y in range(m)):
                t += e_all[S ^ (1 << x)]
        e_all[S] = t
    Dmask = sum(1 << i for i in range(r)); Umask = sum(1 << (r + j) for j in range(s))
    box = [I for I in ids if (I & Dmask) == Dmask and (I & Umask) == 0]
    eP = sum(e_all[I] * e_all[full ^ I] for I in box)
    if eP == 0 or len(box) < 2: return None
    ma = Fraction(sum(e_all[I] * e_all[full ^ I] for I in box if not (I >> a) & 1), eP)
    mb = Fraction(sum(e_all[I] * e_all[full ^ I] for I in box if not (I >> b) & 1), eP)
    return ma, mb


def main():
    third, twothird = Fraction(1, 3), Fraction(2, 3)
    res = collections.Counter()
    fails = []
    for r in range(1, 7):
        for s in range(1, 7):
            for cross in (True, False):
                for beta in range(r + 1):
                    for gamma in range(s + 1):
                        rv = mvalues(r, s, cross, beta, gamma)
                        if rv is None: continue
                        ma, mb = rv
                        fail = (ma < third or ma > twothird) and (mb < third or mb > twothird)
                        key = (cross, beta, gamma)
                        if fail:
                            res['FAIL'] += 1
                            if len(fails) < 20: fails.append((r, s, cross, beta, gamma, str(ma), str(mb)))
                        else:
                            res['NONFAIL'] += 1
    print("结构枚举结果：")
    for k in sorted(res): print(f"  {k}  {res[k]}")
    print("fail 的结构（应只有 cross, β=0, γ=0 且 r,s≥2）:")
    for f in fails: print("   ", f)


if __name__ == '__main__':
    main()
