"""9.D / 9.E 的 n=7,8 定向压力测试（独立实现版）。

与 gen_9de.py 的区别：全部计算走 my_poset.py 的独立实现，
不复用仓库代码，从而构成对论文结论的**独立复核**。

生成：Q = 随机 width-2 偏序（6 或 7 元），z = 第 n 个元素，
      随机下集 D / 上集 U（D∩U=∅）决定 z 的连接。
检验：
  9.D : fail + |box|>=2 + D×D 有不可比对  =>  某 D×D 对平衡
  9.D*: U×U 对偶
  9.E : fail + |box|>=2 + D、U 均链状     =>  存在平衡 Q-对
"""
import random, sys, collections
from fractions import Fraction
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from my_poset import (my_close, my_is_poset, my_count_ext, my_width,
                      my_pval, my_rand_poset, my_induced)

third = Fraction(1, 3); twothird = Fraction(2, 3)


def downsets(Q, m):
    return [S for S in range(1 << m)
            if all(not ((S >> i) & 1) or all((S >> j) & 1 for j in range(m) if Q[j][i])
                   for i in range(m))]


def upsets(Q, m):
    return [S for S in range(1 << m)
            if all(not ((S >> i) & 1) or all((S >> j) & 1 for j in range(m) if Q[i][j])
                   for i in range(m))]


def main():
    n = int(sys.argv[1]); trials = int(sys.argv[2])
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 20261004
    rng = random.Random(seed)
    m = n - 1
    st = collections.Counter(); ex = []
    for _ in range(trials):
        Q = my_rand_poset(m, rng.uniform(0.25, 0.6), rng)
        if not my_is_poset(Q, m) or my_width(Q, m) > 2: continue
        st['width2'] += 1
        DS = downsets(Q, m); US = upsets(Q, m)
        for _ in range(3):
            D = rng.choice(DS); U = rng.choice(US)
            if D & U: continue
            rel = [(i, j) for i in range(m) for j in range(m) if Q[i][j]]
            rel += [(d, m) for d in range(m) if (D >> d) & 1]
            rel += [(m, u) for u in range(m) if (U >> u) & 1]
            L = my_close(n, rel)
            if not my_is_poset(L, n) or my_width(L, n) < 2: continue
            # fail: 所有 z-对不平衡
            if any(third <= my_pval(L, n, m, v) <= twothird
                   for v in range(m) if not (L[v][m] or L[m][v])): continue
            st['fail'] += 1
            # 盒子
            Qidx = [v for v in range(n) if v != m]
            eP = my_count_ext(L, n)
            bm = {}
            for S in range(1 << m):
                if not all(not ((S >> i) & 1) or all((S >> j) & 1 for j in range(m) if Q[j][i])
                           for i in range(m)): continue
                if (S & D) != D or (S & U): continue
                Ilist = [Qidx[k] for k in range(m) if (S >> k) & 1]
                Clist = [Qidx[k] for k in range(m) if not (S >> k) & 1]
                eI = my_count_ext(my_induced(L, Ilist), len(Ilist)) if Ilist else 1
                eC = my_count_ext(my_induced(L, Clist), len(Clist)) if Clist else 1
                bm[S] = Fraction(eI * eC, eP)
            if not bm or len(bm) < 2: continue
            st['fail_multicut'] += 1
            dl = [v for v in range(m) if (D >> v) & 1]
            ul = [v for v in range(m) if (U >> v) & 1]
            dp = [(x, y) for i, x in enumerate(dl) for y in dl[i+1:] if not (L[x][y] or L[y][x])]
            up = [(x, y) for i, x in enumerate(ul) for y in ul[i+1:] if not (L[x][y] or L[y][x])]
            if dp:
                st['9D_avail'] += 1
                if any(third <= my_pval(L, n, x, y) <= twothird for (x, y) in dp): st['9D_ok'] += 1
                else:
                    st['9D_BAD'] += 1
                    if len(ex) < 6: ex.append(('9D', n, rel, D, U, [str(my_pval(L, n, x, y)) for (x, y) in dp]))
            if up:
                st['9Ddual_avail'] += 1
                if any(third <= my_pval(L, n, x, y) <= twothird for (x, y) in up): st['9Ddual_ok'] += 1
                else:
                    st['9Ddual_BAD'] += 1
                    if len(ex) < 12: ex.append(('9D*', n, rel, D, U, [str(my_pval(L, n, x, y)) for (x, y) in up]))
            if not dp and not up:
                st['9E_residual'] += 1
                allp = [(x, y) for i, x in enumerate(range(m)) for y in range(i+1, m) if not (L[x][y] or L[y][x])]
                if any(third <= my_pval(L, n, x, y) <= twothird for (x, y) in allp): st['9E_ok'] += 1
                else:
                    st['9E_BAD'] += 1
                    if len(ex) < 18: ex.append(('9E', n, rel, D, U, [str(my_pval(L, n, x, y)) for (x, y) in allp]))
    print(f"[independent] n={n} trials={trials} seed={seed}")
    for k in sorted(st): print(f"  {k:16s} {st[k]}")
    if ex:
        print("!! 反例样本:")
        for e in ex: print("  ", e)
    else:
        print("零反例")


if __name__ == '__main__':
    main()
