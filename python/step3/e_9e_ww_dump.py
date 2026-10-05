"""dump 1152 个"仅 WW 见证"案例（全部 |W|=4）的窗口结构：
W 的内部偏序、box 形状、m 值与窗口元素对应、WW 见证对的 p 值。"""
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
    shown = collections.Counter()
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
                if len(Wl) != 4: continue
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
                bal_dw = any(third <= pP(d, w) <= twothird
                             for d in Dl for w in Wl if not (Q[d][w] or Q[w][d]))
                bal_wu = any(third <= pP(w, u) <= twothird
                             for w in Wl for u in Ul if not (Q[u][w] or Q[w][u]))
                if bal_dw or bal_wu: continue     # 只看 1152 例
                mkey = ",".join(sorted(str(m[w]) for w in Wl))
                if shown[mkey] >= 2: continue
                shown[mkey] += 1
                print("=" * 70)
                print(f"mask={mask} D={Dl} U={Ul} W={Wl}  eP={eP}  |box|={len(box)}")
                print(f"  m: {[(w, str(m[w])) for w in Wl]}")
                # W 内部偏序
                Winn = [(a, b) for (a, b) in itertools.combinations(Wl, 2) if Q[a][b]]
                print(f"  W 内部关系: {Winn}")
                # box 形状（按 W 子集）
                print(f"  box 点（W 子集）: {sorted(tuple(w for w in Wl if (I >> w) & 1) for I in box)}")
                print(f"  box μ: {[(tuple(w for w in Wl if (I >> w) & 1), str(Fraction(wgt[I], eP))) for I in sorted(box)]}")
                # WW 见证对
                for (a, b) in itertools.combinations(Wl, 2):
                    if Q[a][b] or Q[b][a]: continue
                    p = pP(a, b)
                    print(f"  WW 对 ({a},{b}): p={p}  平衡={third <= p <= twothird}")
                # D/W/U 与 W 的关系
                for w in Wl:
                    below = [d for d in Dl if Q[d][w]]
                    above = [u for u in Ul if Q[w][u]]
                    print(f"  w={w}: D 中低于 w: {below}；U 中高于 w: {above}")


if __name__ == '__main__':
    main()
