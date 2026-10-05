"""验证一般"交叉结构"的 m 值公式（|W|=2 不可比对情形）：

结构：D = d_1≺…≺d_r 链，U = u_1≺…≺u_s 链，W = {a, b}
  a ∥ D（全不可比），a ≺ U（全低于 U）
  b ≻ D（全高于 D），b ∥ U（全不可比）
  a ∥ b
  z：D ≺ z ≺ U，z ∥ a, z ∥ b

推导：
  e(D)=1, e(D∪{a})=r+1, e(D∪{b})=1, e(D∪{a,b})=r+2
  e(Q∖D)=s+2, e(Q∖(D∪{a}))=s+1, e(Q∖(D∪{b}))=1, e(Q∖(D∪{a,b}))=1
  E = rs + 2r + 2s + 6
  m_a = (s+3)/E,  m_b = (rs+2r+2s+3)/E

本脚本对 r,s ∈ [1..6] 直接构造并核对。"""
import sys, itertools, collections
from fractions import Fraction
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from my_poset import my_close, my_is_poset, my_count_ext


def build(r, s):
    """返回 (Q, Dmask, Umask, Wlist, z_idx) —— Q 是 P∖z（不含 z），D/U 链 + 两窗口元。"""
    # 元素编号：D = 0..r-1, U = r..r+s-1, a = r+s, b = r+s+1
    a = r + s; b = r + s + 1
    m = b + 1                       # Q 的元素数（不含 z）
    rel = []
    for i in range(r - 1): rel.append((i, i + 1))             # D 链
    for i in range(s - 1): rel.append((r + i, r + i + 1))      # U 链
    for i in range(r): rel.append((i, b))                      # b ≻ D
    for j in range(s): rel.append((a, r + j))                  # a ≺ U
    Q = my_close(m, rel)
    Dmask = sum(1 << i for i in range(r))
    Umask = sum(1 << (r + j) for j in range(s))
    return Q, Dmask, Umask, [a, b], m


def main():
    st = collections.Counter()
    for r in range(1, 7):
        for s in range(1, 7):
            Q, Dmask, Umask, W, m = build(r, s)
            ids = [S for S in range(1 << m)
                   if all(not ((S >> i) & 1) or all((S >> j) & 1 for j in range(m) if Q[j][i]) for i in range(m))]
            full = (1 << m) - 1
            e_all = [0] * (1 << m); e_all[0] = 1
            for S in range(1, 1 << m):
                tot = 0
                for x in range(m):
                    if not (S >> x) & 1: continue
                    if all(not Q[y][x] or not (S >> y) & 1 for y in range(m)):
                        tot += e_all[S ^ (1 << x)]
                e_all[S] = tot
            box = [I for I in ids if (I & Dmask) == Dmask and (I & Umask) == 0]
            eP = sum(e_all[I] * e_all[full ^ I] for I in box)
            m = {}
            for w in W:
                ev = sum(e_all[I] * e_all[full ^ I] for I in box if not (I >> w) & 1)
                m[w] = Fraction(ev, eP)
            a, b = W
            E = r * s + 2 * r + 2 * s + 6
            fa = Fraction(s + 3, E); fb = Fraction(r * s + r + 2 * s + 3, E)
            if eP != E:
                st['E_BAD'] += 1
                print(f"E MISMATCH r={r} s={s}: eP={eP} E={E}")
            if m[a] != fa:
                st['ma_BAD'] += 1
                print(f"m_a MISMATCH r={r} s={s}: direct={m[a]} formula={fa}")
            if m[b] != fb:
                st['mb_BAD'] += 1
                print(f"m_b MISMATCH r={r} s={s}: direct={m[b]} formula={fb}")
            st['cases'] += 1
    print(f"交叉结构公式验证：{st['cases']} 个 (r,s) 对")
    print(f"  E 不一致 {st['E_BAD']}，m_a 不一致 {st['ma_BAD']}，m_b 不一致 {st['mb_BAD']}")
    print("结论: " + ("全部通过 ✅" if st['E_BAD'] == st['ma_BAD'] == st['mb_BAD'] == 0 else "存在不一致 ❌"))


if __name__ == '__main__':
    main()
