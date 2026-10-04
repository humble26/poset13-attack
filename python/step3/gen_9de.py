"""9.D / 9.E 的 n=7,8 定向压力测试。

论文证据：n<=6 穷举 + 部分采样。本脚本直接生成目标类实例（width-2 的 Q +
随机 z-关系），在 n=7,8 上大量搜索 9.D / 9.E 的反例。

生成方法：
  Q: 随机偏序（随机线性序见证 + 独立加关系 + 传递闭包），筛 width(Q) <= 2
  z: 随机选择下集 D（down-closed）与上集 U（up-closed），D∩U=∅，加 z 的连接
  -> P = Q + z，自动满足 target 类

检验：
  9.D: fail + |box|>=2 + DxD 有不可比对  =>  某 DxD 对平衡
  9.E: fail + |box|>=2 + D、U 均链状     =>  存在平衡 Q-对
"""
import random, sys, collections
from fractions import Fraction

third = Fraction(1, 3); twothird = Fraction(2, 3)


def close(n, rel):
    less = [[False] * n for _ in range(n)]
    for (i, j) in rel: less[i][j] = True
    for k in range(n):
        for i in range(n):
            if less[i][k]:
                ri, rk = less[i], less[k]
                for j in range(n):
                    if rk[j]: ri[j] = True
    return less


def is_poset(less, n):
    for i in range(n):
        if less[i][i]: return False
    return True


def count_ext(less, n):
    dp = [0] * (1 << n); dp[0] = 1
    for S in range(1 << n):
        v = dp[S]
        if not v: continue
        for x in range(n):
            if (S >> x) & 1: continue
            ok = True
            for y in range(n):
                if less[y][x] and not (S >> y) & 1: ok = False; break
            if ok: dp[S | 1 << x] += v
    return dp[-1]


def width(less, n):
    best = 0
    for S in range(1, 1 << n):
        c = bin(S).count('1')
        if c <= best: continue
        els = [x for x in range(n) if (S >> x) & 1]
        if all(not (less[a][b] or less[b][a]) for i, a in enumerate(els) for b in els[i+1:]):
            best = c
    return best


def pval(less, n, x, y):
    L2 = [row[:] for row in less]; L2[x][y] = True
    # 加边 x≺y 的传递闭包增量：{x}∪Pred(x) 全体 → {y}∪Succ(y) 全体
    A = [x] + [k for k in range(n) if less[k][x]]
    B = [y] + [j for j in range(n) if less[y][j]]
    for a in A:
        for b in B:
            if a != b: L2[a][b] = True
    return Fraction(count_ext(L2, n), count_ext(less, n))


def rand_poset(m, p, rng):
    perm = list(range(m)); rng.shuffle(perm)
    less = [[False] * m for _ in range(m)]
    for i in range(m):
        for j in range(i + 1, m):
            if rng.random() < p: less[perm[i]][perm[j]] = True
    # 传递闭包
    for k in range(m):
        for i in range(m):
            if less[i][k]:
                for j in range(m):
                    if less[k][j]: less[i][j] = True
    return less


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    trials = int(sys.argv[2]) if len(sys.argv) > 2 else 200000
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 20261004
    rng = random.Random(seed)
    m = n - 1
    st = collections.Counter()
    ex = []

    def downsets(less, m):
        res = []
        for S in range(1 << m):
            if all((not ((S >> i) & 1)) or all((S >> j) & 1 for j in range(m) if less[j][i])
                   for i in range(m)): res.append(S)
        return res

    def upsets(less, m):
        res = []
        for S in range(1 << m):
            if all((not ((S >> i) & 1)) or all((S >> j) & 1 for j in range(m) if less[i][j])
                   for i in range(m)): res.append(S)
        return res

    for _ in range(trials):
        Q = rand_poset(m, rng.uniform(0.25, 0.6), rng)
        if not is_poset(Q, m): continue
        if width(Q, m) > 2: continue
        st['width2'] += 1
        DS = downsets(Q, m); US = upsets(Q, m)
        # 选若干 (D,U)
        for _ in range(3):
            D = rng.choice(DS); U = rng.choice(US)
            if D & U: continue
            # 构造 n 元 P（z = m）
            rel = [(i, j) for i in range(m) for j in range(m) if Q[i][j]]
            for d in range(m):
                if (D >> d) & 1: rel.append((d, m))
            for u in range(m):
                if (U >> u) & 1: rel.append((m, u))
            L = close(n, rel)
            if not is_poset(L, n): continue
            if width(L, n) < 2: continue
            # z-对是否全失败
            zp = []
            fail = True
            for v in range(m):
                if L[v][m] or L[m][v]: continue
                p = pval(L, n, m, v)
                zp.append(p)
                if third <= p <= twothird: fail = False; break
            if not fail: continue
            st['fail'] += 1
            # 盒子：Q 的理想 I ⊇ D 且 I∩U=∅
            boxmu = {}
            Qidx = [v for v in range(n) if v != m]
            eP = count_ext(L, n)
            for S in range(1 << m):
                # S 是 Q 的理想（down-closed in Q）？
                if not all((not ((S >> i) & 1)) or all((S >> j) & 1 for j in range(m) if Q[j][i])
                           for i in range(m)): continue
                Dm = D; Um = U
                if (S & Dm) != Dm or (S & Um): continue
                Ilist = [Qidx[k] for k in range(m) if (S >> k) & 1]
                Clist = [Qidx[k] for k in range(m) if not (S >> k) & 1]
                eI = count_ext(induced_rows(L, n, m, Ilist), len(Ilist)) if Ilist else 1
                eC = count_ext(induced_rows(L, n, m, Clist), len(Clist)) if Clist else 1
                boxmu[S] = Fraction(eI * eC, eP)
            if not boxmu: continue
            if len(boxmu) < 2: continue
            st['fail_multicut'] += 1
            # DxD 不可比对？
            dl = [v for v in range(m) if (D >> v) & 1]
            ul = [v for v in range(m) if (U >> v) & 1]
            dpairs = [(x, y) for i, x in enumerate(dl) for y in dl[i+1:] if not (L[x][y] or L[y][x])]
            upairs = [(x, y) for i, x in enumerate(ul) for y in ul[i+1:] if not (L[x][y] or L[y][x])]
            # 9.D
            if dpairs:
                st['9D_avail'] += 1
                ok = any(third <= pval(L, n, x, y) <= twothird for (x, y) in dpairs)
                if ok: st['9D_ok'] += 1
                else:
                    st['9D_BAD'] += 1
                    if len(ex) < 5: ex.append(('9D', n, rel, D, U, [str(pval(L, n, x, y)) for (x, y) in dpairs]))
            if upairs:
                st['9D_dual_avail'] += 1
                ok = any(third <= pval(L, n, x, y) <= twothird for (x, y) in upairs)
                if ok: st['9D_dual_ok'] += 1
                else:
                    st['9D_dual_BAD'] += 1
                    if len(ex) < 10: ex.append(('9D*', n, rel, D, U, [str(pval(L, n, x, y)) for (x, y) in upairs]))
            # 9.E 残例：D、U 均链状（无内部不可比对）
            if not dpairs and not upairs:
                st['9E_residual'] += 1
                allpairs = [(x, y) for i, x in enumerate(range(m)) for y in range(i+1, m)
                            if not (L[x][y] or L[y][x])]
                ok = any(third <= pval(L, n, x, y) <= twothird for (x, y) in allpairs)
                if ok: st['9E_ok'] += 1
                else:
                    st['9E_BAD'] += 1
                    if len(ex) < 15: ex.append(('9E', n, rel, D, U, [str(pval(L, n, x, y)) for (x, y) in allpairs]))
    print(f"n={n} trials={trials} seed={seed}")
    for k in sorted(st): print(f"  {k:18s} {st[k]}")
    print("反例样本:")
    for e in ex: print("  ", e)


def induced_rows(less, n, m, sub):
    k = len(sub)
    idx = {v: i for i, v in enumerate(sub)}
    out = [[False] * k for _ in range(k)]
    for i, a in enumerate(sub):
        for j, b in enumerate(sub):
            if less[a][b]: out[i][j] = True
    return out


if __name__ == '__main__':
    main()
