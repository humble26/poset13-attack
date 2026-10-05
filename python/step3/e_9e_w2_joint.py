"""|W|=2 逐案例联合分类：
对每个案例记录 (m1侧, m2侧, 见证类别集合)，看情形分析需要哪些分支。
m_k = μ(w_k ∉ I) = p_P(z, w_k)。
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
    joint = collections.Counter()
    third, twothird = Fraction(1, 3), Fraction(2, 3)
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
                if any(not (Q[a][b] or Q[b][a]) for i, a in enumerate(Dl) for b in Dl[i+1:]): continue
                if any(not (Q[a][b] or Q[b][a]) for i, a in enumerate(Ul) for b in Ul[i+1:]): continue
                box = [I for I in ids if (I & D) == D and (I & U) == 0]
                eP = sum(wgt[I] for I in box)
                if eP == 0 or len(box) < 2: continue
                m = {}
                fail = True
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
                cats = frozenset()
                for (a, b) in PAIRS:
                    if Q[a][b] or Q[b][a]: continue
                    if third <= pP(a, b) <= twothird:
                        def c(t):
                            return 'D' if (D >> t) & 1 else ('U' if (U >> t) & 1 else 'W')
                        cats |= {c(a) + c(b)}
                # m 侧别
                def side(v):
                    return 'lo' if m[v] < third else 'hi'
                # w1≺w2? w2≺w1? ∥?
                rel12 = 'lt' if Q[w1][w2] else ('gt' if Q[w2][w1] else 'inc')
                joint[(rel12, side(w1) + side(w2), tuple(sorted(cats)))] += 1
    print(f"n={n} |W|=2 联合分类（rel12, m1m2侧, 见证类别）:")
    for k in sorted(joint, key=lambda z: -joint[z]):
        print(f"   {k}  {joint[k]}")


if __name__ == '__main__':
    main()
