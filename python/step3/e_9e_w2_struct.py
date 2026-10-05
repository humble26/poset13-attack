"""不可比对 |W|=2 案例的结构统一性检验：
对全部 844 例，记录 (|D|, |U|, eP, μ向量, 各 w 与 D/U 的不可比情况)。
检验"一个 w 与 D 全不可比、另一个与 U 全不可比"是否恒成立。
"""
import sys, itertools, collections
from fractions import Fraction
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from my_poset import my_close, my_is_poset


def width_le2(Q, mm):
    for a, b, c in itertools.combinations(range(mm), 3):
        if not (Q[a][b] or Q[b][a] or Q[a][c] or Q[c][a] or Q[b][c] or Q[c][b]): return False
    return True


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    mm = n - 1
    PAIRS = list(itertools.combinations(range(mm), 2))
    full = (1 << mm) - 1
    third, twothird = Fraction(1, 3), Fraction(2, 3)
    sizes = collections.Counter()
    eP_hist = collections.Counter()
    mu_hist = collections.Counter()
    struct = collections.Counter()
    for mask in range(1 << len(PAIRS)):
        rel = [PAIRS[k] for k in range(len(PAIRS)) if (mask >> k) & 1]
        Q = my_close(mm, rel)
        if not my_is_poset(Q, mm) or not width_le2(Q, mm): continue
        ids = [S for S in range(1 << mm)
               if all(not ((S >> i) & 1) or all((S >> j) & 1 for j in range(mm) if Q[j][i]) for i in range(mm))]
        e_all = [0] * (1 << mm); e_all[0] = 1
        for S in range(1, 1 << mm):
            tot = 0
            for x in range(mm):
                if not (S >> x) & 1: continue
                if all(not Q[y][x] or not (S >> y) & 1 for y in range(mm)):
                    tot += e_all[S ^ (1 << x)]
            e_all[S] = tot
        wgt = {I: e_all[I] * e_all[full ^ I] for I in ids}
        USM = [S for S in range(1 << mm)
               if all(not ((S >> i) & 1) or all((S >> j) & 1 for j in range(mm) if Q[i][j]) for i in range(mm))]
        for D in ids:
            for U in USM:
                if D & U: continue
                Wl = [v for v in range(mm) if not ((D >> v) & 1) and not ((U >> v) & 1)]
                if len(Wl) != 2: continue
                w1, w2 = Wl
                if Q[w1][w2] or Q[w2][w1]: continue
                Dl = [v for v in range(mm) if (D >> v) & 1]
                Ul = [v for v in range(mm) if (U >> v) & 1]
                if any(not (Q[a][b] or Q[b][a]) for i, a in enumerate(Dl) for b in Dl[i+1:]): continue
                if any(not (Q[a][b] or Q[b][a]) for i, a in enumerate(Ul) for b in Ul[i+1:]): continue
                box = [I for I in ids if (I & D) == D and (I & U) == 0]
                eP = sum(wgt[I] for I in box)
                if eP == 0 or len(box) < 2: continue
                fail = True
                m = {}
                for v in Wl:
                    ev = sum(wgt[I] for I in box if not (I >> v) & 1)
                    m[v] = Fraction(ev, eP)
                    if third <= m[v] <= twothird: fail = False
                if not fail: continue
                sizes[(len(Dl), len(Ul))] += 1
                eP_hist[eP] += 1
                # μ 向量（按 W 子集）
                muvec = []
                for I in sorted(box):
                    sub = tuple(w for w in Wl if (I >> w) & 1)
                    muvec.append((sub, str(Fraction(wgt[I], eP))))
                mu_hist[tuple(muvec)] += 1
                # 各 w 与 D/U 的不可比情况
                for w in Wl:
                    incD = sum(1 for d in Dl if not (Q[d][w] or Q[w][d]))
                    incU = sum(1 for u in Ul if not (Q[u][w] or Q[w][u]))
                    struct[(incD, incU)] += 1
    print(f"n={n} 不可比对 |W|=2 案例总数检查")
    print("(|D|,|U|) 分布:", dict(sizes))
    print("eP 分布:", dict(eP_hist))
    print("各 w 的 (|D不可比|, |U不可比|) 分布:")
    for k in sorted(struct, key=lambda z: -struct[z]): print(f"    {k}  {struct[k]}")
    print("μ 向量分布（前 8）:")
    for k in sorted(mu_hist, key=lambda z: -mu_hist[z])[:8]: print(f"    {k}  {mu_hist[k]}")


if __name__ == '__main__':
    main()
