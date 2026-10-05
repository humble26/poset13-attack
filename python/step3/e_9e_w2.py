"""|W| = 2 的 9.E 情形刻画（全穷举）：
W = {w1, w2}。box ⊆ {D, D∪{w1}, D∪{w2}, D∪{w1,w2}}（理想性筛选后可能更少）。

要刻画：
 (a) 平衡对落在哪类（DW / WU / WW）—— 特别是 (w1,w2) 本身平衡吗？
 (b) p_P(d_i, w) 作为 i 的函数是否**仿射**（9.G 的闭式是仿射的特例）？
 (c) 若仿射：拟合端点，找 (μ, N) 的闭式。
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
    st = collections.Counter()
    affine_ok = 0; affine_bad = 0
    boxshapes = collections.Counter()
    witness = collections.Counter()
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
                Dl = [v for v in range(mm) if (D >> v) & 1]
                Ul = [v for v in range(mm) if (U >> v) & 1]
                # D,U 链状（9.E 残余）
                if any(not (Q[a][b] or Q[b][a]) for i, a in enumerate(Dl) for b in Dl[i+1:]): continue
                if any(not (Q[a][b] or Q[b][a]) for i, a in enumerate(Ul) for b in Ul[i+1:]): continue
                box = [I for I in ids if (I & D) == D and (I & U) == 0]
                eP = sum(wgt[I] for I in box)
                if eP == 0 or len(box) < 2: continue
                # fail
                isfail = True
                for v in Wl:
                    ev = sum(wgt[I] for I in box if not (I >> v) & 1)
                    if 3 * ev >= eP and 3 * ev <= 2 * eP: isfail = False; break
                if not isfail: continue
                st['cases'] += 1
                key = tuple(sorted(('1' if (I >> w1) & 1 else '0') + ('1' if (I >> w2) & 1 else '0') for I in box))
                boxshapes[key] += 1
                # M-帧
                gs = 0; gaps = []
                for pos in poslists:
                    t = max([pos[v] for v in Dl], default=-1)
                    nu = min([pos[v] for v in Ul], default=mm)
                    g = nu - t
                    gaps.append(g)
                    if g > 0: gs += g
                if gs == 0: continue
                # 平衡对分类
                def pP(a, b):
                    num = 0
                    for pos, g in zip(poslists, gaps):
                        if g > 0 and pos[a] < pos[b]: num += g
                    return Fraction(num, gs)
                cats = set()
                for (a, b) in PAIRS:
                    if Q[a][b] or Q[b][a]: continue
                    p = pP(a, b)
                    if Fraction(1, 3) <= p <= Fraction(2, 3):
                        def c(t):
                            return 'D' if (D >> t) & 1 else ('U' if (U >> t) & 1 else 'W')
                        cats.add(c(a) + c(b))
                for x in cats: witness[x] += 1
                # 仿射性：p_P(d_i, w1) 沿 D 的不可比后缀是否仿射
                Dchain = chain_sort(Q, Dl)
                inc1 = [i for i, v in enumerate(Dchain) if not (Q[v][w1] or Q[w1][v])]
                if len(inc1) >= 3:
                    vals = [(i, pP(Dchain[i], w1)) for i in inc1]
                    (i0, p0), (i1, p1), (i2, p2) = vals[0], vals[1], vals[2]
                    # 仿射检验：p0 - 2 p1 + p2 == 0 ?
                    if p0 - 2 * p1 + p2 == 0: affine_ok += 1
                    else: affine_bad += 1
    print(f"n={n} |W|=2 残余案例: {st['cases']}")
    print("平衡对类别:", dict(witness))
    print("box 形状分布（w1w2 ∈ {00,10,01,11} 的多重集）:", dict(boxshapes))
    print(f"仿射性（p(d_i,w1) 二阶差分为零）: ok {affine_ok} / bad {affine_bad}")


if __name__ == '__main__':
    main()
