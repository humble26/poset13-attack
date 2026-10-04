"""边际分析（n=6 严格穷举）：9.D 案例里"最好的那对"离区间边界多近？
margin = max_{(x,y) ∈ D×D 不可比对} min(p_P − 1/3, 2/3 − p_P)
margin 越小 => 9.D 越勉强 => 反例最可能出现在那里。
同时记录最紧案例的完整数据。"""
import sys, itertools, collections
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from my_poset import my_close, my_is_poset


def width_le2(Q, mm):
    for a, b, c in itertools.combinations(range(mm), 3):
        if not (Q[a][b] or Q[b][a] or Q[a][c] or Q[c][a] or Q[b][c] or Q[c][b]): return False
    return True


def main():
    n = 6; mm = n - 1
    PAIRS = list(itertools.combinations(range(mm), 2))
    PIDX = {p: k for k, p in enumerate(PAIRS)}
    full = (1 << mm) - 1
    margins = collections.Counter()
    tight = []
    minmargin = 9
    seen = 0
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
        def rec(placed, seq):
            if len(seq) == mm: exts.append(tuple(seq)); return
            for x in range(mm):
                if (placed >> x) & 1 or (pred[x] & ~placed): continue
                rec(placed | 1 << x, seq + [x])
        rec(0, [])
        poslists = []
        for seq in exts:
            pos = [0] * mm
            for i, v in enumerate(seq): pos[v] = i
            poslists.append(pos)
        DSM = ids
        USM = [S for S in range(1 << mm)
               if all(not ((S >> i) & 1) or all((S >> j) & 1 for j in range(mm) if Q[i][j]) for i in range(mm))]
        for D in DSM:
            for U in USM:
                if D & U: continue
                Dl = [v for v in range(mm) if (D >> v) & 1]
                Ul = [v for v in range(mm) if (U >> v) & 1]
                dd = [(x, y) for i, x in enumerate(Dl) for y in Dl[i+1:] if not (Q[x][y] or Q[y][x])]
                if not dd: continue
                box = [I for I in ids if (I & D) == D and (I & U) == 0]
                eP = sum(w[I] for I in box)
                if eP == 0 or len(box) < 2: continue
                isfail = True
                for v in range(mm):
                    if (D >> v) & 1 or (U >> v) & 1: continue
                    ev = sum(w[I] for I in box if not (I >> v) & 1)
                    if 3 * ev >= eP and 3 * ev <= 2 * eP: isfail = False; break
                if not isfail: continue
                # 只关心 D×D 对：p = Σ_{pos} gap·[x<_pos y] / Σ_{pos} gap
                gs = 0
                gaps = []
                for pos in poslists:
                    t = max([pos[v] for v in Dl], default=-1)
                    nu = min([pos[v] for v in Ul], default=mm)
                    g = nu - t
                    gaps.append(g)
                    if g > 0: gs += g
                if gs == 0: continue
                best = None
                for (x, y) in dd:
                    num = 0
                    for pos, g in zip(poslists, gaps):
                        if g > 0 and pos[x] < pos[y]: num += g
                    if 3 * num >= gs and 3 * num <= 2 * gs:
                        mar = min(3 * num - gs, 2 * gs - 3 * num)
                        if best is None or mar > best: best = mar
                if best is None:
                    margins['BAD'] += 1
                    tight.append((mask, D, U, str(dd)))
                    continue
                seen += 1
                margins[float(best / gs)] += 1
                if best / gs < minmargin:
                    minmargin = best / gs
                    tight.append((n, mask, D, U, str(best) + '/' + str(gs), [(x, y) for (x, y) in dd]))
                    if len(tight) > 500: tight.pop(0)
    print(f"n=6 9.D 案例 {seen}; 其中 BAD（无反例但需确认）: {margins['BAD']}")
    print(f"最小 margin: {minmargin}")
    print("margin 最小的 10 个:")
    for k in sorted([k for k in margins if isinstance(k, float)])[:10]:
        print(f"   margin≈{k:.5f}  案例数 {margins[k]}")
    print("最紧案例:")
    for t in tight[-3:]: print("   ", t)


if __name__ == '__main__':
    main()
