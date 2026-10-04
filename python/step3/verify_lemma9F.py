"""验证引理 9.F（D-带）在 n=6,7 目标类实例上成立：
(i)  {l : b_l ∥ x} 是 B 链上的连续区间（行凸）；对偶列凸
(ii) l -> p_P(x,b_l) 单调不减；k -> p_P(a_k,y) 单调不增
(iii) p_P(x,b_l)=0 当 b_l ≺ x；=1 当 x ≺ b_l（区间外钉死）
"""
import random, sys, collections
from fractions import Fraction
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from my_poset import (my_close, my_is_poset, my_count_ext, my_width,
                      my_pval, my_rand_poset)

third = Fraction(1, 3); twothird = Fraction(2, 3)
n = int(sys.argv[1]) if len(sys.argv) > 1 else 6
trials = int(sys.argv[2]) if len(sys.argv) > 2 else 200000
rng = random.Random(20261004)
m = n - 1


def chain_parts(Q, mm):
    """枚举链划分 (A,B)，A,B 为元素列表（已按链序排列）。"""
    els = list(range(mm)); out = []
    for r in range(mm + 1):
        for S in __import__("itertools").combinations(els, r):
            T = [e for e in els if e not in S]
            def chain(C):
                C = list(C)
                return all(Q[a][b] or Q[b][a] for i, a in enumerate(C) for b in C[i+1:])
            if chain(S) and chain(T):
                A = sorted(S, key=lambda v: sum(1 for u in S if Q[u][v]))
                B = sorted(T, key=lambda v: sum(1 for u in T if Q[u][v]))
                out.append((A, B))
    return out


st = collections.Counter()
for _ in range(trials):
    Q = my_rand_poset(m, rng.uniform(0.25, 0.6), rng)
    if not my_is_poset(Q, m) or my_width(Q, m) > 2: continue
    parts = chain_parts(Q, m)
    if not parts: continue
    A, B = rng.choice(parts)
    if len(A) < 1 or len(B) < 1: continue
    # 构造 P：随机 D,U
    DS = [S for S in range(1 << m) if all(not ((S >> i) & 1) or all((S >> j) & 1 for j in range(m) if Q[j][i]) for i in range(m))]
    US = [S for S in range(1 << m) if all(not ((S >> i) & 1) or all((S >> j) & 1 for j in range(m) if Q[i][j]) for i in range(m))]
    D = rng.choice(DS); U = rng.choice(US)
    if D & U: continue
    rel = [(i, j) for i in range(m) for j in range(m) if Q[i][j]]
    rel += [(d, m) for d in range(m) if (D >> d) & 1]
    rel += [(m, u) for u in range(m) if (U >> u) & 1]
    # 用 Q 的索引作为 P 的前 m 个元素（z = m）
    L = my_close(n, rel)
    if not my_is_poset(L, n) or my_width(L, n) < 2: continue
    st['instances'] += 1
    # (i)/(ii)/(iii) 对每个 x ∈ A
    t, s = len(B), len(A)
    for x in A:
        incomp = [l for l in range(t) if not (L[x][B[l]] or L[B[l]][x])]
        below = [l for l in range(t) if L[B[l]][x]]
        above = [l for l in range(t) if L[x][B[l]]]
        # (i) 区间性
        if incomp and incomp != list(range(min(incomp), max(incomp) + 1)):
            st['F_i_viol'] += 1
        # (ii) 单调性（在 B 全链上）
        vals = [my_pval(L, n, x, B[l]) for l in range(t)]
        if any(vals[i] > vals[i+1] for i in range(t-1)):
            st['F_ii_viol'] += 1
        # (iii) 钉死
        for l in below:
            if my_pval(L, n, x, B[l]) != 0: st['F_iii0_viol'] += 1
        for l in above:
            if my_pval(L, n, x, B[l]) != 1: st['F_iii1_viol'] += 1
    # 对偶：每个 y ∈ B
    for y in B:
        incomp = [k for k in range(s) if not (L[y][A[k]] or L[A[k]][y])]
        if incomp and incomp != list(range(min(incomp), max(incomp) + 1)):
            st['F_dual_i_viol'] += 1
        vals = [my_pval(L, n, A[k], y) for k in range(s)]
        if any(vals[i] < vals[i+1] for i in range(s-1)):  # k -> p 单调不增
            st['F_dual_ii_viol'] += 1
print(f"n={n} trials={trials}  实例 {st['instances']}")
for k in sorted(st):
    if k != 'instances': print(f"  {k:16s} {st[k]}")
viol = sum(v for k, v in st.items() if k != 'instances')
print("引理 9.F 机器验证：" + ("全部通过 ✅" if viol == 0 else f"存在 {viol} 处违反 ❌"))
