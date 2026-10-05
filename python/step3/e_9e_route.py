"""路线引理 R 的检验（任意 |W| 的 9.E 残余案例，全穷举）：

R（锐利版）: 对每个 w ∈ W：
   m(w) < 1/3  ⟹  ∃ d ∈ D, d ∥ w,  p_P(d, w) ∈ [1/3, 2/3]
   m(w) > 2/3  ⟹  ∃ u ∈ U, u ∥ w,  p_P(w, u) ∈ [1/3, 2/3]
（若 R 成立，则 9.E 的证明统一为：fail ⟹ 每个 w 的 m 都在区间外 ⟹ 任取一个
  lo 的 w 走 D-路线或 hi 的 w 走 U-路线。）

R'（弱化版）: m(w) < 1/3 ⟹ 存在平衡对落在 D×W（不指定 w）。

同时记录 R 失败时 m 的分布——看是否卡在 (1/4, 1/3) 这种"不够低"的区域。
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
    st = collections.Counter()
    bad_R = []
    bad_Rp = []
    m_lo_hist = collections.Counter()
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
                st['cases'] += 1
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
                # R 锐利版
                for w in Wl:
                    if m[w] < third:
                        m_lo_hist[m[w]] += 1
                        cands = [d for d in Dl if not (Q[d][w] or Q[w][d])]
                        if cands:
                            st['lo_with_cand'] += 1
                            if any(third <= pP(d, w) <= twothird for d in cands):
                                st['R_lo_ok'] += 1
                            else:
                                st['R_lo_BAD'] += 1
                                if len(bad_R) < 5:
                                    bad_R.append(('lo', mask, D, U, w, str(m[w]),
                                                  [str(pP(d, w)) for d in cands]))
                    elif m[w] > twothird:
                        cands = [u for u in Ul if not (Q[u][w] or Q[w][u])]
                        if cands:
                            st['hi_with_cand'] += 1
                            if any(third <= pP(w, u) <= twothird for u in cands):
                                st['R_hi_ok'] += 1
                            else:
                                st['R_hi_BAD'] += 1
                                if len(bad_R) < 10:
                                    bad_R.append(('hi', mask, D, U, w, str(m[w]),
                                                  [str(pP(w, u)) for u in cands]))
                # R' 弱化版：案例级
                has_lo = any(m[w] < third for w in Wl)
                has_hi = any(m[w] > twothird for w in Wl)
                bal_dw = any(third <= pP(d, w) <= twothird
                             for d in Dl for w in Wl if not (Q[d][w] or Q[w][d]))
                bal_wu = any(third <= pP(w, u) <= twothird
                             for w in Wl for u in Ul if not (Q[u][w] or Q[w][u]))
                if has_lo and not bal_dw and not bal_wu:
                    st['Rp_BAD'] += 1
                    if len(bad_Rp) < 5: bad_Rp.append(('lo-nodw', mask, D, U, [str(m[w]) for w in Wl]))
                if has_hi and not bal_wu and not bal_dw:
                    st['Rp_BAD'] += 1
                if not (bal_dw or bal_wu):
                    st['no_witness_at_all'] += 1
                    st[f'nowit|W|={len(Wl)}'] += 1
                    st[f'nowit_mvals|{",".join(sorted(str(m[w]) for w in Wl))}'] += 1
    print(f"n={n} 9.E 残余案例: {st['cases']}")
    print(f"R 锐利版: lo 有候选 {st['lo_with_cand']}，成立 {st['R_lo_ok']}，失败 {st['R_lo_BAD']}")
    print(f"          hi 有候选 {st['hi_with_cand']}，成立 {st['R_hi_ok']}，失败 {st['R_hi_BAD']}")
    print(f"R' 案例级失败: {st['Rp_BAD']}   完全无见证: {st['no_witness_at_all']}")
    print("无 DW/WU 见证案例的 |W| 分布:")
    for k in sorted(st):
        if k.startswith('nowit|W|'): print(f"    {k}  {st[k]}")
    print("无 DW/WU 见证案例的 m 值组合（前 10）:")
    cnt = 0
    for k in sorted(st, key=lambda z: -st[z]):
        if k.startswith('nowit_mvals|') and cnt < 10:
            print(f"    {k[len('nowit_mvals|'):]}  {st[k]}"); cnt += 1
    print("m(w)<1/3 的取值分布（前 12）:")
    for k in sorted(m_lo_hist)[:12]: print(f"    m={k}  {m_lo_hist[k]}")
    for b in bad_R: print("  R-BAD:", b)
    for b in bad_Rp: print("  R'-BAD:", b)


if __name__ == '__main__':
    main()
