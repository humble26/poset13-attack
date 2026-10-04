"""9.E 机制刻画（全穷举）：残余情形（D、U 均为链）下，
平衡 Q-对落在哪些 D/W/U 类别里？位置有无规律？

类别按端点所属：D / W / U（无需链划分，仅用 D、U 两个集合即可定义）。
D×D 与 U×U 在残余情形无不可比对对，故类别 ∈ {DW, WW, WU, DU}（有序）。

同时记录：
  - 每个 9.E 案例里"全部"平衡对及其类别
  - 平衡对的最紧 margin
  - |W|、|box|、切割 (i*, j*) 的代理量（在此用 box 形状）
"""
import sys, itertools, collections
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
    PIDX = {p: k for k, p in enumerate(PAIRS)}
    full = (1 << mm) - 1
    cat = collections.Counter()
    sizes = collections.Counter()
    margins = collections.Counter()
    ncase = 0
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
        w = {I: e_all[I] * e_all[full ^ I] for I in ids}
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
        omasks = []
        for seq in exts:
            pos = [0] * mm
            for i, v in enumerate(seq): pos[v] = i
            poslists.append(pos)
            mk = 0
            for k, (a, b) in enumerate(PAIRS):
                if pos[a] < pos[b]: mk |= 1 << k
            omasks.append(mk)
        USM = [S for S in range(1 << mm)
               if all(not ((S >> i) & 1) or all((S >> j) & 1 for j in range(mm) if Q[i][j]) for i in range(mm))]
        for D in ids:
            for U in USM:
                if D & U: continue
                Dl = [v for v in range(mm) if (D >> v) & 1]
                Ul = [v for v in range(mm) if (U >> v) & 1]
                dd = [(x, y) for i, x in enumerate(Dl) for y in Dl[i+1:] if not (Q[x][y] or Q[y][x])]
                uu = [(x, y) for i, x in enumerate(Ul) for y in Ul[i+1:] if not (Q[x][y] or Q[y][x])]
                if dd or uu: continue                     # 只看 9.E 残余情形
                box = [I for I in ids if (I & D) == D and (I & U) == 0]
                eP = sum(w[I] for I in box)
                if eP == 0 or len(box) < 2: continue
                isfail = True
                for v in range(mm):
                    if (D >> v) & 1 or (U >> v) & 1: continue
                    ev = sum(w[I] for I in box if not (I >> v) & 1)
                    if 3 * ev >= eP and 3 * ev <= 2 * eP: isfail = False; break
                if not isfail: continue
                allp = [(a, b) for (a, b) in PAIRS if not (Q[a][b] or Q[b][a])]
                if not allp: continue
                gs = 0; gaps = []
                for pos in poslists:
                    t = max([pos[v] for v in Dl], default=-1)
                    nu = min([pos[v] for v in Ul], default=mm)
                    g = nu - t
                    gaps.append(g)
                    if g > 0: gs += g
                if gs == 0: continue
                ncase += 1
                W = [v for v in range(mm) if not ((D >> v) & 1) and not ((U >> v) & 1)]
                sizes[len(W)] += 1
                best = None
                for (a, b) in allp:
                    num = 0
                    for pos, g in zip(poslists, gaps):
                        if g > 0 and pos[a] < pos[b]: num += g
                    if 3 * num >= gs and 3 * num <= 2 * gs:
                        def c(t):
                            return 'D' if (D >> t) & 1 else ('U' if (U >> t) & 1 else 'W')
                        cat[c(a) + c(b)] += 1
                        cat[c(b) + c(a)] += 1
                        mar = min(3 * num - gs, 2 * gs - 3 * num)
                        if best is None or mar > best: best = mar
                if best is not None:
                    margins[float(best / gs)] += 1
    print(f"n={n} 9.E 残余案例: {ncase}")
    print("平衡对的 D/W/U 类别分布（有序对，双向计入）:")
    for k in sorted(cat, key=lambda z: -cat[z]): print(f"   {k}: {cat[k]}")
    print("|W| 分布:", dict(sizes))
    print("最紧 margin:", min(margins) if margins else None)
    print("margin 最小 6 个:", sorted(margins)[:6])


if __name__ == '__main__':
    main()
