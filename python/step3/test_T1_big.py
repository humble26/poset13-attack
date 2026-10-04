"""把更强的猜想 T1 拿到 n=7,8 上压力测试：
T1: fail + multicut + D×D 有不可比对  =>  存在 (x,y) ∈ D×D 不可比对，
    使 p_D(x,y) ∈ [1/3,2/3] **且** p_P(x,y) ∈ [1/3,2/3]。
（注意 T1 ⟹ 9.D；若 T1 在 n=7/8 也零反例，则它是比 9.D 更好的攻击目标。）
"""
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
    rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 7)
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
            if any(third <= my_pval(L, n, m, v) <= twothird
                   for v in range(m) if not (L[v][m] or L[m][v])): continue   # 需 fail
            # 多切割
            nb = sum(1 for S in range(1 << m)
                     if all(not ((S >> i) & 1) or all((S >> j) & 1 for j in range(m) if Q[j][i]) for i in range(m))
                     and (S & D) == D and not (S & U))
            if nb < 2: continue
            dl = [v for v in range(m) if (D >> v) & 1]
            band = [(x, y) for i, x in enumerate(dl) for y in dl[i+1:] if not (L[x][y] or L[y][x])]
            if not band: continue
            st['cases'] += 1
            # D 作为独立偏序集
            Dl = my_induced(L, dl)
            pos = {v: i for i, v in enumerate(dl)}
            ok = False; okp = []
            for (x, y) in band:
                pd = my_pval(Dl, len(dl), pos[x], pos[y])
                pp = my_pval(L, n, x, y)
                if third <= pd <= twothird:
                    okp.append((x, y, pd, pp))
                    if third <= pp <= twothird: ok = True
            if ok: st['T1_ok'] += 1
            else:
                st['T1_BAD'] += 1
                if len(ex) < 6: ex.append((n, rel, D, U, [(x, y, str(pd), str(pp)) for (x, y, pd, pp) in okp]))
    print(f"[T1] n={n} trials={trials}  9.D 案例 {st['cases']}")
    print(f"  T1 成立: {st['T1_ok']}   反例: {st['T1_BAD']}")
    for e in ex: print("  BAD:", e)


if __name__ == '__main__':
    main()
