"""链盒（w1≺w2）结构分类：枚举 (α1, α2, ε1, ε2)，算 m1,m2，找 fail 条件。
参数：
  α1 = w1 的 D-不可比后缀起点（w1 ∥ d_{α1}..d_r，w1 ≻ d_1..d_{α1−1}），α1 ∈ [1, r+1]
  α2 = w2 的 D-不可比后缀起点（α2 ≥ α1）
  ε1 = w1 ≺ u_{ε1}..u_s（w1 的 U-上界后缀起点），ε1 ∈ [1, s+1]（ε1=s+1 ⟹ w1 与 U 全不可比）
  ε2 = w2 的 U-上界后缀起点
"""
import sys, itertools, collections
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


def chainbox(r, s, a1, a2, e1, e2):
    """元素：D=0..r-1, U=r..r+s-1, w1=r+s, w2=r+s+1。返回 Q 或 None（无效）。"""
    w1 = r + s; w2 = r + s + 1; m = w2 + 1
    rel = []
    for i in range(r - 1): rel.append((i, i + 1))
    for i in range(s - 1): rel.append((r + i, r + i + 1))
    rel.append((w1, w2))                          # w1 ≺ w2
    for j in range(a1 - 1): rel.append((j, w1))   # w1 ≻ d_1..d_{a1−1}
    for j in range(a2 - 1): rel.append((j, w2))   # w2 ≻ d_1..d_{a2−1}
    for j in range(e1 - 1, s): rel.append((w1, r + j))   # w1 ≺ u_{e1}..u_s
    for j in range(e2 - 1, s): rel.append((w2, r + j))   # w2 ≺ u_{e2}..u_s
    Q = my_close(m, rel)
    if not my_is_poset(Q, m): return None
    for x, y, z in itertools.combinations(range(m), 3):
        if not (Q[x][y] or Q[y][x] or Q[x][z] or Q[z][x] or Q[y][z] or Q[z][y]):
            return None
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
    if eP == 0 or len(box) < 2: return None
    m1 = Fraction(sum(ea[I] * ea[full ^ I] for I in box if not (I >> w1) & 1), eP)
    m2 = Fraction(sum(ea[I] * ea[full ^ I] for I in box if not (I >> w2) & 1), eP)
    return m1, m2


def main():
    third, twothird = Fraction(1, 3), Fraction(2, 3)
    fails = collections.Counter()
    examples = []
    for r in range(1, 5):
        for s in range(1, 5):
            for a1 in range(1, r + 2):
                for a2 in range(a1, r + 2):
                    for e1 in range(1, s + 2):
                        for e2 in range(1, s + 2):
                            rv = chainbox(r, s, a1, a2, e1, e2)
                            if rv is None: continue
                            m1, m2 = rv
                            if (m1 < third or m1 > twothird) and (m2 < third or m2 > twothird):
                                fails[(r - a1 + 1, s - e1 + 1)] += 1
                                if len(examples) < 12:
                                    examples.append((r, s, a1, a2, e1, e2, str(m1), str(m2)))
    print("fail 案例按 (p1=r−α1+1, q1=s−ε1+1) 分组:")
    for k in sorted(fails, key=lambda z: -fails[z]): print(f"   {k}  {fails[k]}")
    print("示例:")
    for e in examples: print("   ", e)


if __name__ == '__main__':
    main()
