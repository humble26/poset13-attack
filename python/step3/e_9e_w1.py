"""|W| = 1 的 9.E 情形：验证候选闭式。

设 D、U 为链（残余情形），W = {w}，box = {D, D∪{w}}（多切割），
m := μ(D) = p_P(z,w)（条带公式），q := 1 − m。

猜测：对 d ∈ D ∩ {与 w 不可比}，设 D = (d_1 ≺ … ≺ d_r)，
不可比区间为 [α, β]（1-based），N = β − α + 2，则

    p_P(d_i, w) = m + (1 − m)·(β − i + 1)/N      （i ∈ [α, β]）

即：w 在 D∪{w} 中有 N 个等概率插入槽，p 是"槽位比例的仿射像"。
本脚本用 M-帧引擎核对（n=7 抽目标类，只看 |W| = 1 且 D,U 链状且 fail 且多切割）。
"""
import sys, itertools, collections
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from my_poset import my_close, my_is_poset


def width_le2(Q, mm):
    for a, b, c in itertools.combinations(range(mm), 3):
        if not (Q[a][b] or Q[b][a] or Q[a][c] or Q[c][a] or Q[b][c] or Q[c][b]): return False
    return True


def chain_indices(Q, S):
    """S 是链时返回其元素按序列表，否则 None。"""
    els = list(S)
    if all(Q[a][b] or Q[b][a] for i, a in enumerate(els) for b in els[i+1:]):
        return sorted(els, key=lambda v: sum(1 for u in els if Q[u][v]))
    return None


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    mm = n - 1
    PAIRS = list(itertools.combinations(range(mm), 2))
    full = (1 << mm) - 1
    st = collections.Counter()
    examples = []
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
                W = [v for v in range(mm) if not ((D >> v) & 1) and not ((U >> v) & 1)]
                if len(W) != 1: continue
                w = W[0]
                Dl = [v for v in range(mm) if (D >> v) & 1]
                Ul = [v for v in range(mm) if (U >> v) & 1]
                if not chain_indices(Q, Dl) or not chain_indices(Q, Ul): continue   # D,U 链状
                dd = [(x, y) for i, x in enumerate(Dl) for y in Dl[i+1:] if not (Q[x][y] or Q[y][x])]
                uu = [(x, y) for i, x in enumerate(Ul) for y in Ul[i+1:] if not (Q[x][y] or Q[y][x])]
                if dd or uu: continue
                box = [I for I in ids if (I & D) == D and (I & U) == 0]
                eP = sum(wgt[I] for I in box)
                if eP == 0 or len(box) != 2: continue
                isfail = True
                ev = sum(wgt[I] for I in box if not (I >> w) & 1)
                if 3 * ev >= eP and 3 * ev <= 2 * eP: isfail = False
                if not isfail: continue
                m_num, m_den = ev, eP          # m = mu(D) = ev/eP
                st['cases'] += 1
                # 交叉分类：w 是否与 D / U 有不可比元素；m 落在哪一侧
                has_d = any(not (Q[v][w] or Q[w][v]) for v in Dl)
                has_u = any(not (Q[v][w] or Q[w][v]) for v in Ul)
                mr = 'lo' if 3 * ev < eP else ('hi' if 3 * ev > 2 * eP else 'mid')
                st[f'cfg|hasD={has_d}|hasU={has_u}|m={mr}'] += 1
                # D-路线可用 = has_d 且 m<1/3；U-路线可用（对偶）= has_u 且 m>2/3
                if not ((has_d and mr == 'lo') or (has_u and mr == 'hi')):
                    st['NOT_COVERED_BY_ROUTES'] += 1
                # M-帧：真实 p_P(d,w)
                gs = 0; gaps = []
                for pos in poslists:
                    t = max([pos[v] for v in Dl], default=-1)
                    nu = min([pos[v] for v in Ul], default=mm)
                    g = nu - t
                    gaps.append(g)
                    if g > 0: gs += g
                # D 链序
                Dchain = chain_indices(Q, Dl)
                # —— U-路线公式核对：p_P(w,u_j) = (1−m) + m(j−α'+1)/N'
                Uchain = chain_indices(Q, Ul)
                Uinc = [j for j, v in enumerate(Uchain) if not (Q[v][w] or Q[w][v])]
                if Uinc and Uinc == list(range(min(Uinc), max(Uinc) + 1)):
                    a2, b2 = min(Uinc), max(Uinc)
                    N2 = (b2 - a2 + 1) + 1
                    for j in range(a2, b2 + 1):
                        u = Uchain[j]
                        num = 0
                        for pos, g in zip(poslists, gaps):
                            if g > 0 and pos[w] < pos[u]: num += g
                        pn = (m_den - m_num) * N2 + m_num * (j - a2 + 1)
                        pd = m_den * N2
                        from fractions import Fraction
                        if Fraction(num, gs) != Fraction(pn, pd):
                            st['UFORMULA_MISMATCH'] += 1
                        else:
                            st['uformula_ok'] += 1
                # 与 w 不可比的 D 元素索引区间
                inc = [i for i, v in enumerate(Dchain) if not (Q[v][w] or Q[w][v])]
                if not inc:
                    st['no_d_incomp'] += 1
                    continue
                al, be = min(inc), max(inc)      # 0-based
                if inc != list(range(al, be + 1)):
                    st['interval_violation'] += 1
                    continue
                N = (be - al + 1) + 1            # = β−α+2（1-based 记法）
                # —— 关键推论验算 ——
                # 若 w 与 D 有不可比、但与 U 全可比（¬hasU），应恰有 m = 1/(N+1)
                if not has_u:
                    st['check_m_num1'] += 1
                    if m_num * (N + 1) == m_den: st['m_eq_1_over_N1'] += 1
                    else: st['m_neq_1_over_N1'] += 1
                # 又：不可比区间应是 D 的后缀（be == len(Dchain)−1）
                if be != len(Dchain) - 1: st['not_suffix'] += 1
                else: st['is_suffix'] += 1
                for i in range(al, be + 1):
                    d = Dchain[i]
                    num = 0
                    for pos, g in zip(poslists, gaps):
                        if g > 0 and pos[d] < pos[w]: num += g
                    pred_num = m_num * N + (m_den - m_num) * (be - i + 1)
                    pred_den = m_den * N
                    # 约分后比较
                    from fractions import Fraction
                    if Fraction(num, gs) != Fraction(pred_num, pred_den):
                        st['FORMULA_MISMATCH'] += 1
                        if len(examples) < 4:
                            examples.append((Dchain, w, i, (num, gs), (pred_num, pred_den)))
                    else:
                        st['formula_ok'] += 1
    print(f"n={n} |W|=1 残余案例: {st['cases']}")
    print(f"  公式吻合: {st['formula_ok']}   不吻合: {st['FORMULA_MISMATCH']}")
    print(f"  U-路线公式吻合: {st['uformula_ok']}   不吻合: {st['UFORMULA_MISMATCH']}")
    print(f"  D 无与 w 不可比元素: {st['no_d_incomp']}")
    print(f"  不可比集非区间: {st['interval_violation']}")
    print("  交叉分类（hasD/hasU/m 侧）:")
    for k in sorted(st):
        if k.startswith('cfg|'): print(f"     {k}  {st[k]}")
    print(f"  两条路线都覆盖不到: {st['NOT_COVERED_BY_ROUTES']}")
    print(f"  [推论] ¬hasU 案例 {st['check_m_num1']} 个：m == 1/(N+1) 成立 {st['m_eq_1_over_N1']}，不成立 {st['m_neq_1_over_N1']}")
    print(f"  [结构] D 不可比区间是 D 的后缀：{st['is_suffix']} 是 / {st['not_suffix']} 否")
    for e in examples: print("   MISMATCH:", e)


if __name__ == '__main__':
    main()
