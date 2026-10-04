"""n=7 全穷举证伪引擎（9.D / 9.D对偶 / 9.E）。

思路（避免逐个 DP）：
  Q = 6 元 width<=2 偏序，枚举 2^15 个 mask（恒等排列为线性扩展，WLOG），
  再枚举所有 (D, U)  ——  D 下集、U 上集、D∩U=∅，等价于 z 的全部合法位置。
  P = Q + z（z 记作第 6 号元素）。

  * e(P) 与 e(P+z≺v) 用引理 4.1（理想分解）算：
        e(P)          = Σ_{I ⊇ D, I∩U=∅}      e(I)·e(Q∖I)
        e(P+z≺v)      = Σ_{I ⊇ D, I∩(U∪{v})=∅} e(I)·e(Q∖I)
    e(I)、e(Q∖I) 对每个 Q 预计算一次（I 取遍 Q 的理想）。
  * 任意 Q-对 (x,y) 的 p_P 用引理 C'（M-帧重加权）算：
        p_P(x,y) = Σ_{M ∈ ext(Q)} gap(M)·[x <_M y] / Σ_M gap(M)
    gap(M) = max(0, ν(M) − t(M))，t = D 在 M 中的最大位次，ν = U 的最小位次。
    每个 Q 预计算 ext(Q) 及每条的 21-bit "x 先于 y" 指示掩码。

自检：--selftest 用直接 DP（my_pval）对拍 M-帧公式。
"""
import sys, time, itertools, collections
from fractions import Fraction
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from my_poset import my_close, my_is_poset, my_count_ext, my_width, my_pval, my_induced

m = 6
NPAIRS = [(i, j) for i, j in itertools.combinations(range(m), 2)]
PAIRIDX = {p: k for k, p in enumerate(NPAIRS)}


def ideals_of(Q, mm):
    """Q 的全部理想（down-closed bitmask）。"""
    return [S for S in range(1 << mm)
            if all(not ((S >> i) & 1) or all((S >> j) & 1 for j in range(mm) if Q[j][i])
                   for i in range(mm))]


def downsets(Q, mm):
    return ideals_of(Q, mm)


def upsets(Q, mm):
    return [S for S in range(1 << mm)
            if all(not ((S >> i) & 1) or all((S >> j) & 1 for j in range(mm) if Q[i][j])
                   for i in range(mm))]


def extensions(Q, mm):
    """枚举 Q 的全部线性扩展（作为 tuple），回溯。"""
    pred = [0] * mm
    for i in range(mm):
        for j in range(mm):
            if Q[j][i]: pred[i] |= 1 << j
    out = []

    def rec(placed, seq):
        if len(seq) == mm:
            out.append(tuple(seq)); return
        for x in range(mm):
            if (placed >> x) & 1: continue
            if pred[x] & ~placed: continue
            rec(placed | 1 << x, seq + [x])
    rec(0, [])
    return out


def order_masks(ext):
    """每条扩展 -> 21-bit 掩码：第 k 位为 1 当且仅当 NPAIRS[k] 的先者在后者的前面。"""
    masks = []
    for seq in ext:
        pos = {v: i for i, v in enumerate(seq)}
        mk = 0
        for k, (a, b) in enumerate(NPAIRS):
            if pos[a] < pos[b]: mk |= 1 << k
        masks.append(mk)
    return masks


def p_from_gap(gap_sum, pair_sums, eP, k, sign_rev=False):
    """返回 (num, den) 形式的 p 值，已约成与 eP 同分母；sign_rev 表示取 1-p。"""
    num = pair_sums[k]
    if sign_rev: num = gap_sum - num
    return num, gap_sum


def main():
    t0 = time.time()
    st = collections.Counter()
    bad = []
    nmax = 1 << len(NPAIRS)
    for mask in range(nmax):
        rel = [NPAIRS[k] for k in range(len(NPAIRS)) if (mask >> k) & 1]
        Q = my_close(m, rel)
        if not my_is_poset(Q, m) or my_width(Q, m) > 2: continue
        st['Q'] += 1
        # 预计算：理想 + e(I)e(Q∖I)
        ids = ideals_of(Q, m)
        w = {}
        for I in ids:
            Il = [v for v in range(m) if (I >> v) & 1]
            Cl = [v for v in range(m) if not (I >> v) & 1]
            eI = my_count_ext(my_induced(Q, Il), len(Il)) if Il else 1
            eC = my_count_ext(my_induced(Q, Cl), len(Cl)) if Cl else 1
            w[I] = eI * eC
        # 预计算：ext(Q) 与顺序掩码
        exts = extensions(Q, m)
        omasks = order_masks(exts)
        DS = downsets(Q, m); US = upsets(Q, m)
        for D in DS:
            for U in US:
                if D & U: continue
                # box ideals
                box = [I for I in ids if (I & D) == D and (I & U) == 0]
                eP = sum(w[I] for I in box)
                if eP == 0: continue
                Wmask = ((1 << m) - 1) & ~D & ~U
                Wv = [v for v in range(m) if (Wmask >> v) & 1]
                # fail 判定
                isfail = True
                for v in Wv:
                    ev = sum(w[I] for I in box if not (I >> v) & 1)
                    # p = ev/eP ∈ [1/3,2/3] ?
                    if 3 * ev >= eP and 3 * ev <= 2 * eP:
                        isfail = False; break
                if not isfail: continue
                st['fail'] += 1
                if len(box) < 2: continue      # 单切割，已由 9.B 关闭
                st['fail_multicut'] += 1
                # M-帧：gap(M) 与 pair 累积
                gap_sum = 0
                pair_sum = [0] * len(NPAIRS)
                for seq, om in zip(exts, omasks):
                    pos = {v: i for i, v in enumerate(seq)}
                    t = max([pos[v] for v in range(m) if (D >> v) & 1], default=-1)
                    nu = min([pos[v] for v in range(m) if (U >> v) & 1], default=m)
                    g = nu - t
                    if g <= 0: continue
                    gap_sum += g
                    for k in range(len(NPAIRS)):
                        if (om >> k) & 1: pair_sum[k] += g
                if gap_sum != eP:
                    st['GAP_MISMATCH'] += 1
                    continue
                # 9.D : D×D 不可比对
                dl = [v for v in range(m) if (D >> v) & 1]
                dd = [(x, y) for i, x in enumerate(dl) for y in dl[i + 1:] if not (Q[x][y] or Q[y][x])]
                ul = [v for v in range(m) if (U >> v) & 1]
                uu = [(x, y) for i, x in enumerate(ul) for y in ul[i + 1:] if not (Q[x][y] or Q[y][x])]
                ddok = None
                if dd:
                    st['9D_avail'] += 1
                    ok = False
                    for (x, y) in dd:
                        k = PAIRIDX[(min(x, y), max(x, y))]
                        num = pair_sum[k] if x < y else gap_sum - pair_sum[k]
                        if 3 * num >= gap_sum and 3 * num <= 2 * gap_sum: ok = True; break
                    if ok: st['9D_ok'] += 1
                    else:
                        st['9D_BAD'] += 1
                        if len(bad) < 5: bad.append(('9D', mask, rel, D, U))
                if uu:
                    st['9Ddual_avail'] += 1
                    ok = False
                    for (x, y) in uu:
                        k = PAIRIDX[(min(x, y), max(x, y))]
                        num = pair_sum[k] if x < y else gap_sum - pair_sum[k]
                        if 3 * num >= gap_sum and 3 * num <= 2 * gap_sum: ok = True; break
                    if ok: st['9Ddual_ok'] += 1
                    else:
                        st['9Ddual_BAD'] += 1
                        if len(bad) < 10: bad.append(('9D*', mask, rel, D, U))
                if not dd and not uu:
                    st['9E_residual'] += 1
                    ok = False
                    for k, (a, b) in enumerate(NPAIRS):
                        if Q[a][b] or Q[b][a]: continue
                        num = pair_sum[k]
                        if 3 * num >= gap_sum and 3 * num <= 2 * gap_sum: ok = True; break
                    if ok: st['9E_ok'] += 1
                    else:
                        st['9E_BAD'] += 1
                        if len(bad) < 15: bad.append(('9E', mask, rel, D, U))
        if mask % 4096 == 0:
            print(f"  ... mask {mask}/{nmax}  {time.time()-t0:.0f}s  {dict(st)}", flush=True)
    print("=" * 60)
    for k in sorted(st): print(f"  {k:18s} {st[k]}")
    if bad:
        print("!!! 反例:")
        for b in bad: print("   ", b)
    else:
        print("零反例")
    print(f"总耗时 {time.time()-t0:.0f}s")


def selftest():
    """对拍 M-帧公式 vs 直接 DP。"""
    import random
    rng = random.Random(1)
    checked = 0; mism = 0
    while checked < 120:
        mask = rng.randrange(1 << len(NPAIRS))
        rel = [NPAIRS[k] for k in range(len(NPAIRS)) if (mask >> k) & 1]
        Q = my_close(m, rel)
        if not my_is_poset(Q, m) or my_width(Q, m) > 2: continue
        DS = downsets(Q, m); US = upsets(Q, m)
        D = rng.choice(DS); U = rng.choice(US)
        if D & U: continue
        # 构造 n=7 的 P
        n = 7
        relP = [(i, j) for i in range(m) for j in range(m) if Q[i][j]]
        relP += [(v, m) for v in range(m) if (D >> v) & 1]
        relP += [(m, v) for v in range(m) if (U >> v) & 1]
        P = my_close(n, relP)
        if not my_is_poset(P, n): continue
        eP = my_count_ext(P, n)
        exts = extensions(Q, m); omasks = order_masks(exts)
        gap_sum = 0; pair_sum = [0] * len(NPAIRS)
        for seq, om in zip(exts, omasks):
            pos = {v: i for i, v in enumerate(seq)}
            t = max([pos[v] for v in range(m) if (D >> v) & 1], default=-1)
            nu = min([pos[v] for v in range(m) if (U >> v) & 1], default=m)
            g = nu - t
            if g <= 0: continue
            gap_sum += g
            for k in range(len(NPAIRS)):
                if (om >> k) & 1: pair_sum[k] += g
        if gap_sum != eP:
            print("gap_sum != eP", mask, D, U, gap_sum, eP); mism += 1; continue
        for k, (a, b) in enumerate(NPAIRS):
            if Q[a][b] or Q[b][a]: continue
            from fractions import Fraction as F
            p_direct = my_pval(P, n, a, b)
            p_frame = F(pair_sum[k], gap_sum)
            if p_direct != p_frame:
                mism += 1
                if mism <= 5: print("P-MISMATCH", mask, D, U, (a, b), p_direct, p_frame)
        checked += 1
    print(f"自检：{checked} 实例，p 不一致 {mism} —— " + ("通过 ✅" if mism == 0 else "失败 ❌"))


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--selftest':
        selftest()
    else:
        main()
