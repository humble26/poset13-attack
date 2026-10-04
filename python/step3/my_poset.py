"""独立实现的最小偏序集工具箱（刻意不复用仓库代码，用于交叉验证）。
命名全部加前缀，避免与仓库 e5_ideal / e9_anatomy / e11_corner 的名称冲突。"""
from fractions import Fraction


def my_close(n, rel):
    less = [[False] * n for _ in range(n)]
    for (i, j) in rel: less[i][j] = True
    for k in range(n):
        for i in range(n):
            if less[i][k]:
                for j in range(n):
                    if less[k][j]: less[i][j] = True
    return less


def my_is_poset(less, n):
    return all(not less[i][i] for i in range(n))


def my_count_ext(less, n):
    """线性扩展计数：DP over 下集。"""
    dp = [0] * (1 << n); dp[0] = 1
    for S in range(1 << n):
        v = dp[S]
        if not v: continue
        for x in range(n):
            if (S >> x) & 1: continue
            if all(not less[y][x] or ((S >> y) & 1) for y in range(n)):
                dp[S | 1 << x] += v
    return dp[-1]


def my_width(less, n):
    best = 0
    for S in range(1, 1 << n):
        c = bin(S).count('1')
        if c <= best: continue
        els = [x for x in range(n) if (S >> x) & 1]
        if all(not (less[a][b] or less[b][a]) for i, a in enumerate(els) for b in els[i+1:]):
            best = c
    return best


def my_pval(less, n, x, y):
    """p_P(x,y) = e(P+x≺y)/e(P)：显式加边 + 传递闭包。"""
    L2 = [row[:] for row in less]
    L2[x][y] = True
    A = [x] + [k for k in range(n) if less[k][x]]
    B = [y] + [j for j in range(n) if less[y][j]]
    for a in A:
        for b in B:
            if a != b: L2[a][b] = True
    return Fraction(my_count_ext(L2, n), my_count_ext(less, n))


def my_rand_poset(m, p, rng):
    perm = list(range(m)); rng.shuffle(perm)
    less = [[False] * m for _ in range(m)]
    for i in range(m):
        for j in range(i + 1, m):
            if rng.random() < p: less[perm[i]][perm[j]] = True
    for k in range(m):
        for i in range(m):
            if less[i][k]:
                for j in range(m):
                    if less[k][j]: less[i][j] = True
    return less


def my_induced(less, sub):
    k = len(sub)
    return [[less[a][b] for b in sub] for a in sub]
