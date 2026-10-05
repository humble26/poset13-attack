"""验证链盒情形（w1 ≺ w2）的显式公式：

  e(I3) = (α2−α1)(r−α2+2) + (r−α2+1)(r−α2+2)/2
  #(d_i ≺ w1) = 分段：
     i ≥ α2:  (r−i+1)(r−i+2)/2
     i < α2:  (α2−1−i)(r−α2+2) + (r−α2+1)(r−α2+2)/2

其中 I3 = D∪{w1,w2}，D 链长 r，w1 不可比后缀 [α1, r]，w2 不可比后缀 [α2, r]，α2 ≥ α1。
推导：先插 w1（槽 s ∈ [α1−1, r]），再插 w2（须在 max(pos(d_{α2−1}), pos(w1)) 之后）。
"""
import sys, itertools, collections
from fractions import Fraction
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from my_poset import my_close, my_is_poset, my_induced, my_count_ext


def width_le2(Q, mm):
    for a, b, c in itertools.combinations(range(mm), 3):
        if not (Q[a][b] or Q[b][a] or Q[a][c] or Q[c][a] or Q[b][c] or Q[c][b]): return False
    return True


def chain_sort(Q, els):
    return sorted(els, key=lambda v: sum(1 for u in els if Q[u][v]))


def formula(r, a1, a2, i):
    """返回 (e(I3), #(d_i ≺ w1))，1-based 链指标，i ∈ [α1, r]。
    e3 = (α2−α1)(r−α2+2) + (r−α2+2)(r−α2+3)/2
    i ≥ α2:  cnt = (r−i+1)(r−i+2)/2
    i < α2:  cnt = (α2−1−i)(r−α2+2) + (r−α2+2)(r−α2+3)/2
    """
    e3 = (a2 - a1) * (r - a2 + 2) + (r - a2 + 2) * (r - a2 + 3) // 2
    if i >= a2:
        cnt = (r - i + 1) * (r - i + 2) // 2
    else:
        cnt = (a2 - 1 - i) * (r - a2 + 2) + (r - a2 + 2) * (r - a2 + 3) // 2
    return e3, cnt


def main():
    # 直接构造小实例验证：D 链 + w1 + w2，给定 α1, α2
    st = collections.Counter()
    for r in range(1, 6):
        for a1 in range(1, r + 2):
            for a2 in range(a1, r + 2):
                # 构造偏序：d_1..d_r 链；w1 ≺ w2；d_j ≺ w1 (j < α1)；d_j ≺ w2 (j < α2)
                els = r + 2          # d_0..d_{r-1}, w1=idx r, w2=idx r+1
                rel = [(i, i + 1) for i in range(r - 1)]        # D 链
                rel.append((r, r + 1))                            # w1 ≺ w2
                for j in range(a1 - 1): rel.append((j, r))       # d_j ≺ w1
                for j in range(a2 - 1): rel.append((j, r + 1))   # d_j ≺ w2
                L = my_close(els, rel)
                assert my_is_poset(L, els)
                e3_direct = my_count_ext(L, els)
                e3_formula, _ = formula(r, a1, a2, a1)
                if e3_direct != e3_formula:
                    st['e3_BAD'] += 1
                    if st['e3_BAD'] <= 3:
                        print(f"e3 MISMATCH r={r} a1={a1} a2={a2}: direct={e3_direct} formula={e3_formula}")
                else:
                    st['e3_ok'] += 1
                # 逐 i 验证 #(d_i ≺ w1)
                for i in range(a1, r + 1):
                    d = i - 1
                    # d ≺ w1 的扩展数：加关系后计数
                    L2 = my_close(els, rel + [(d, r)])
                    cnt_direct = my_count_ext(L2, els)
                    _, cnt_formula = formula(r, a1, a2, i)
                    if cnt_direct != cnt_formula:
                        st['cnt_BAD'] += 1
                        if st['cnt_BAD'] <= 5:
                            print(f"cnt MISMATCH r={r} a1={a1} a2={a2} i={i}: direct={cnt_direct} formula={cnt_formula}")
                    else:
                        st['cnt_ok'] += 1
    print(f"e(I3) 公式: ok {st['e3_ok']} / bad {st['e3_BAD']}")
    print(f"#(d_i≺w1) 公式: ok {st['cnt_ok']} / bad {st['cnt_BAD']}")
    print("结论: " + ("全部通过 ✅" if st['e3_BAD'] == 0 and st['cnt_BAD'] == 0 else "存在不一致 ❌"))


if __name__ == '__main__':
    main()
