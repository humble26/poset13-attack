"""链盒（w1 ≺ w2）情形的 9.E 完整验证 + 结论提炼。

在真实偏序集上（不止链），验证链盒情形的守恒：
  W = {w1, w2}, w1 ≺ w2, D,U 链, fail, 多切割
  ⟹ 存在平衡对（且落在 D×W 或 W×W）。

并检验核心数值断言：沿 w1 的 D-不可比后缀 [α1, r]，p_P(d_i, w1) 从 <1/3 单调升到 >2/3
（或端点落区间），从而必扫过 [1/3,2/3]。
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
    sweep_bad = []
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
                if not Q[w1][w2] and not Q[w2][w1]: continue   # 只看链盒（可比的 w 对）
                if Q[w2][w1]: w1, w2 = w2, w1                 # 归一化 w1 ≺ w2
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
                st['chainbox_cases'] += 1
                # 是否有平衡对
                anybal = any(third <= pP(a, b) <= twothird
                             for (a, b) in PAIRS if not (Q[a][b] or Q[b][a]))
                if anybal: st['chainbox_bal'] += 1
                else:
                    st['chainbox_BAD'] += 1
                    if len(sweep_bad) < 3: sweep_bad.append(('nobal', mask, D, U))
                # w1 的 D-不可比后缀扫描
                Dchain = chain_sort(Q, Dl)
                inc = [i for i, v in enumerate(Dchain) if not (Q[v][w1] or Q[w1][v])]
                if inc:
                    vals = [pP(Dchain[i], w1) for i in inc]
                    lo = min(vals); hi = max(vals)
                    sweep_ok = (lo < third and hi > twothird) or (third <= lo <= twothird) or (third <= hi <= twothird)
                    if sweep_ok: st['w1_sweep_ok'] += 1
                    else:
                        st['w1_sweep_fail'] += 1
                        if len(sweep_bad) < 8: sweep_bad.append(('sweep', mask, D, U, [str(v) for v in vals]))
    print(f"n={n} 链盒（|W|=2 可比对）残余案例: {st['chainbox_cases']}")
    print(f"  有平衡对: {st['chainbox_bal']}   无: {st['chainbox_BAD']}")
    print(f"  w1 不可比后缀扫过区间: {st['w1_sweep_ok']}   失败: {st['w1_sweep_fail']}")
    for b in sweep_bad: print("  BAD:", b)


if __name__ == '__main__':
    main()
