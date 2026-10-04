"""n 全穷举证伪引擎（9.D / 9.D对偶 / 9.E）—— 通用版。

枚举 Q = P∖z（n-1 元，width<=2）的全部实例：
  Q 用 2^(m(m-1)/2) 个 mask（恒等排列为线性扩展，WLOG），筛 poset + width<=2；
  再枚举全部 (D,U)（D 下集、U 上集、D∩U=∅）。
  P = Q + z（z = 第 m 号元素）。

  * e(P)、e(P+z≺v) 走理想分解（引理 4.1），e(I)、e(Q∖I) 每个 Q 预计算一次。
  * Q-对 p_P 走 M-帧重加权（引理 C'）：width(Q)<=2 ⟹ |ext(Q)| <= C(m, floor(m/2))，很小。
  * width(Q) <= 2 用"无三元反链"判定（比枚举全部子集快得多）。

用法: python exhaustive_n.py 7 | 8
"""
import sys, time, itertools, collections
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from my_poset import my_close, my_is_poset, my_count_ext, my_width, my_induced


def width_le2(Q, mm):
    """width(Q) <= 2  <=>  不存在三元反链。"""
    for a, b, c in itertools.combinations(range(mm), 3):
        if not (Q[a][b] or Q[b][a] or Q[a][c] or Q[c][a] or Q[b][c] or Q[c][b]):
            return False
    return True


def ideals_of(Q, mm):
    return [S for S in range(1 << mm)
            if all(not ((S >> i) & 1) or all((S >> j) & 1 for j in range(mm) if Q[j][i])
                   for i in range(mm))]


def extensions(Q, mm):
    pred = [0] * mm
    for i in range(mm):
        for j in range(mm):
            if Q[j][i]: pred[i] |= 1 << j
    out = []
    def rec(placed, seq):
        if len(seq) == mm: out.append(tuple(seq)); return
        for x in range(mm):
            if (placed >> x) & 1 or (pred[x] & ~placed): continue
            rec(placed | 1 << x, seq + [x])
    rec(0, [])
    return out


def main():
    n = int(sys.argv[1])
    lo = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    hi = int(sys.argv[3]) if len(sys.argv) > 3 else None
    mm = n - 1
    PAIRS = list(itertools.combinations(range(mm), 2))
    PIDX = {p: k for k, p in enumerate(PAIRS)}
    t0 = time.time()
    st = collections.Counter()
    bad = []
    nmask = 1 << len(PAIRS)
    if hi is None: hi = nmask
    for mask in range(lo, hi):
        rel = [PAIRS[k] for k in range(len(PAIRS)) if (mask >> k) & 1]
        Q = my_close(mm, rel)
        if not my_is_poset(Q, mm) or not width_le2(Q, mm): continue
        st['Q'] += 1
        ids = ideals_of(Q, mm)
        # e_all[S] = 诱导子偏序 Q|S 的线性扩展数（单次子集 DP）
        full = (1 << mm) - 1
        e_all = [0] * (1 << mm)
        e_all[0] = 1
        for S in range(1, 1 << mm):
            tot = 0
            for x in range(mm):
                if not (S >> x) & 1: continue
                okmin = True
                for y in range(mm):
                    if Q[y][x] and ((S >> y) & 1): okmin = False; break
                if okmin: tot += e_all[S ^ (1 << x)]
            e_all[S] = tot
        w = {I: e_all[I] * e_all[full ^ I] for I in ids}
        exts = extensions(Q, mm)
        omasks = []
        poslists = []
        for seq in exts:
            pos = [0] * mm
            for i, v in enumerate(seq): pos[v] = i
            poslists.append(pos)
            mk = 0
            for k, (a, b) in enumerate(PAIRS):
                if pos[a] < pos[b]: mk |= 1 << k
            omasks.append(mk)
        DS = ideals_of(Q, mm)
        US = [S for S in range(1 << mm)
              if all(not ((S >> i) & 1) or all((S >> j) & 1 for j in range(mm) if Q[i][j])
                     for i in range(mm))]
        for D in DS:
            for U in US:
                if D & U: continue
                box = [I for I in ids if (I & D) == D and (I & U) == 0]
                eP = 0
                for I in box: eP += w[I]
                if eP == 0: continue
                Wv = [v for v in range(mm) if not ((D >> v) & 1) and not ((U >> v) & 1)]
                isfail = True
                for v in Wv:
                    ev = 0
                    for I in box:
                        if not (I >> v) & 1: ev += w[I]
                    if 3 * ev >= eP and 3 * ev <= 2 * eP: isfail = False; break
                if not isfail: continue
                st['fail'] += 1
                if len(box) < 2: continue
                st['fail_multicut'] += 1
                gap_sum = 0
                ps = [0] * len(PAIRS)
                Dl = [v for v in range(mm) if (D >> v) & 1]
                Ul = [v for v in range(mm) if (U >> v) & 1]
                for pos, om in zip(poslists, omasks):
                    t = -1
                    for v in Dl:
                        p = pos[v]
                        if p > t: t = p
                    nu = mm                      # U=∅ 时 ν 取 m（插入槽 k ∈ [t+1, m]）
                    for v in Ul:
                        p = pos[v]
                        if p < nu: nu = p
                    g = nu - t
                    if g <= 0: continue
                    gap_sum += g
                    for k in range(len(PAIRS)):
                        if (om >> k) & 1: ps[k] += g
                if gap_sum != eP:
                    st['GAP_MISMATCH'] += 1; continue
                dl = Dl
                dd = [(x, y) for i, x in enumerate(dl) for y in dl[i+1:] if not (Q[x][y] or Q[y][x])]
                ul = Ul
                uu = [(x, y) for i, x in enumerate(ul) for y in ul[i+1:] if not (Q[x][y] or Q[y][x])]
                def bal(x, y):
                    k = PIDX[(min(x, y), max(x, y))]
                    num = ps[k] if x < y else gap_sum - ps[k]
                    return 3 * num >= gap_sum and 3 * num <= 2 * gap_sum
                if dd:
                    st['9D_avail'] += 1
                    if any(bal(x, y) for (x, y) in dd): st['9D_ok'] += 1
                    else:
                        st['9D_BAD'] += 1
                        if len(bad) < 5: bad.append(('9D', mask, D, U))
                if uu:
                    st['9Ddual_avail'] += 1
                    if any(bal(x, y) for (x, y) in uu): st['9Ddual_ok'] += 1
                    else:
                        st['9Ddual_BAD'] += 1
                        if len(bad) < 10: bad.append(('9D*', mask, D, U))
                if not dd and not uu:
                    st['9E_residual'] += 1
                    allp = [(a, b) for (a, b) in PAIRS if not (Q[a][b] or Q[b][a])]
                    if not allp:
                        # Q 是链：无 Q-对；此时 fail 若发生即守恒引理本身的反例
                        st['Qchain_residual'] += 1
                        continue
                    if any(bal(a, b) for (a, b) in allp):
                        st['9E_ok'] += 1
                    else:
                        st['9E_BAD'] += 1
                        if len(bad) < 15: bad.append(('9E', mask, D, U))
        if st['Q'] % 2000 == 0 and mask % 65536 == 0:
            print(f"  ... mask {mask}/{nmask} Q={st['Q']} fm={st['fail_multicut']} {time.time()-t0:.0f}s", flush=True)
    print(f"=== 分片 [{lo},{hi}) ===")
    for k in sorted(st): print(f"  {k:18s} {st[k]}")
    print("反例:" if bad else "零反例")
    for b in bad: print("   ", b)
    print(f"总耗时 {time.time()-t0:.0f}s")
    print("SUMMARY|%d|%d|%s" % (lo, hi, dict(st)))


if __name__ == '__main__':
    main()
