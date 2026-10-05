"""|W|=2 不可比对情形：验证关键数值断言，寻找可证结构。

W = {w1, w2}, w1 ∥ w2, box 4 点。m_k = μ(w_k ∉ I) = p_P(z, w_k)。
在失败下 m1, m2 ∈ [1/3,2/3] 之外。检验：
 (A) 见证对是否总可用 m 的"低侧"选择确定
 (B) m1 与 m2 的联合分布（是否有 m1+m2=1 或类似刚性关系）
 (C) p_P(d_i, w1) 的仿射性/单调性沿 D-后缀
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


def chain_sort(Q, els):
    return sorted(els, key=lambda v: sum(1 for u in els if Q[u][v]))


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    mm = n - 1
    PAIRS = list(itertools.combinations(range(mm), 2))
    full = (1 << mm) - 1
    third, twothird = Fraction(1, 3), Fraction(2, 3)
    st = collections.Counter()
    mjoint = collections.Counter()
    affine_ok = affine_bad = 0
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
        pred = [0] * mm
        for i in range(mm):
            for j in range(mm):
                if Q[j][i]: pred[i] |= 1 << j
        exts = []
        def rec(pl, sq):
            if len(sq) == mm: exts.append(tuple(sq)); return
            for x in range(mm):
                if (pl >> x) & 1 or (pred[x] & ~pl): continue
                rec(pl | 1 << x, sq + [x])
        rec(0, [])
        poslists = []
        for seq in exts:
            pos = [0] * mm
            for i, v in enumerate(seq): pos[v] = i
            poslists.append(pos)
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
                gs = 0; gaps = []
                for pos in poslists:
                    t = max([pos[v] for v in Dl], default=-1)
                    nu = min([pos[v] for v in Ul], default=mm)
                    g = nu - t
                    gaps.append(g)
                    if g > 0: gs += g
                if gs == 0: continue
                def pP(a, b):
                    num = 0
                    for pos, g in zip(poslists, gaps):
                        if g > 0 and pos[a] < pos[b]: num += g
                    return Fraction(num, gs)
                st['cases'] += 1
                mjoint[(str(m[w1]), str(m[w2]))] += 1
                # p_P(d_i, w1) 沿 D-后缀仿射？
                Dchain = chain_sort(Q, Dl)
                inc1 = [i for i, v in enumerate(Dchain) if not (Q[v][w1] or Q[w1][v])]
                if len(inc1) >= 3:
                    vals = [pP(Dchain[i], w1) for i in inc1]
                    if vals[0] - 2 * vals[1] + vals[2] == 0: affine_ok += 1
                    else: affine_bad += 1
    print(f"n={n} 不可比对 |W|=2 案例: {st['cases']}")
    print("(m1,m2) 联合分布（前 15）:")
    for k in sorted(mjoint, key=lambda z: -mjoint[z])[:15]: print(f"    {k}  {mjoint[k]}")
    print(f"仿射性 p(d_i,w1): ok {affine_ok} / bad {affine_bad}")


if __name__ == '__main__':
    main()
