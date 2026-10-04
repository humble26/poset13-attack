"""关键判别：9.D 是否需要 fail 假设？
对 P = Q+z, width(Q)<=2, multicut box, DxD 有不可比对：
  (i)  在**全部**实例上，是否某 DxD 对平衡？       （无需 fail）
  (ii) 只在 fail 实例上，是否某 DxD 对平衡？        （9.D 原陈述）
若 (i) 也零反例 => fail 多余，9.D 是纯盒子结构命题。"""
import random, sys, collections
from fractions import Fraction
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from my_poset import (my_close, my_is_poset, my_count_ext, my_width,
                      my_pval, my_rand_poset, my_induced)

third = Fraction(1, 3); twothird = Fraction(2, 3)


def downsets(Q, m):
    return [S for S in range(1 << m) if all(not ((S >> i) & 1) or all((S >> j) & 1 for j in range(m) if Q[j][i]) for i in range(m))]


def upsets(Q, m):
    return [S for S in range(1 << m) if all(not ((S >> i) & 1) or all((S >> j) & 1 for j in range(m) if Q[i][j]) for i in range(m))]


def main():
    n = int(sys.argv[1]); trials = int(sys.argv[2])
    rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
    m = n - 1
    st = collections.Counter(); ex = []
    for _ in range(trials):
        Q = my_rand_poset(m, rng.uniform(0.25, 0.6), rng)
        if not my_is_poset(Q, m) or my_width(Q, m) > 2: continue
        DS = downsets(Q, m); US = upsets(Q, m)
        for _ in range(3):
            D = rng.choice(DS); U = rng.choice(US)
            if D & U: continue
            rel = [(i, j) for i in range(m) for j in range(m) if Q[i][j]]
            rel += [(d, m) for d in range(m) if (D >> d) & 1]
            rel += [(m, u) for u in range(m) if (U >> u) & 1]
            L = my_close(n, rel)
            if not my_is_poset(L, n) or my_width(L, n) < 2: continue
            fails = not any(third <= my_pval(L, n, m, v) <= twothird
                            for v in range(m) if not (L[v][m] or L[m][v]))
            # 盒子
            Qidx = [v for v in range(n) if v != m]
            eP = my_count_ext(L, n)
            nb = 0; dpairs = []
            for S in range(1 << m):
                if not all(not ((S >> i) & 1) or all((S >> j) & 1 for j in range(m) if Q[j][i]) for i in range(m)): continue
                if (S & D) != D or (S & U): continue
                nb += 1
            if nb < 2: continue
            dl = [v for v in range(m) if (D >> v) & 1]
            dpairs = [(x, y) for i, x in enumerate(dl) for y in dl[i+1:] if not (L[x][y] or L[y][x])]
            if not dpairs: continue
            ok = any(third <= my_pval(L, n, x, y) <= twothird for (x, y) in dpairs)
            if fails:
                st['fail_cases'] += 1
                if ok: st['fail_9D_ok'] += 1
                else:
                    st['fail_9D_BAD'] += 1
                    if len(ex) < 5: ex.append(('fail', n, rel, D, U, [str(my_pval(L, n, x, y)) for (x, y) in dpairs]))
            else:
                st['nofail_cases'] += 1
                if ok: st['nofail_9D_ok'] += 1
                else:
                    st['nofail_9D_BAD'] += 1
                    if len(ex) < 10: ex.append(('NOFAIL', n, rel, D, U, [str(my_pval(L, n, x, y)) for (x, y) in dpairs]))
    print(f"n={n} trials={trials}")
    for k in sorted(st): print(f"  {k:16s} {st[k]}")
    for e in ex: print("  ", e)


if __name__ == '__main__':
    main()
