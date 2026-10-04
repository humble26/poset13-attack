"""测关键下界/上界猜想（fail 案例，D×D 带内）：
 L1  p_P(x,y) >= 1/3 对**所有** D×D 不可比对对成立？
 U1  带内 min p_P <= 2/3 ？（若 L1 且 U1 => 9.D 立即成立）
 L1' 对偶：p_P <= 2/3 对所有 U×U 不可比对对成立？
同时给出带内 min/max p_P 的分布。
"""
import sys, collections, random
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
    rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 3)
    m = n - 1
    st = collections.Counter(); ex = []
    minv = Fraction(9); maxofmin = Fraction(0)
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
                   for v in range(m) if not (L[v][m] or L[m][v])): continue
            nb = sum(1 for S in range(1 << m)
                     if all(not ((S >> i) & 1) or all((S >> j) & 1 for j in range(m) if Q[j][i]) for i in range(m))
                     and (S & D) == D and not (S & U))
            if nb < 2: continue
            dl = [v for v in range(m) if (D >> v) & 1]
            dp = [(x, y) for i, x in enumerate(dl) for y in dl[i+1:] if not (L[x][y] or L[y][x])]
            if not dp: continue
            st['cases'] += 1
            vals = [my_pval(L, n, x, y) for (x, y) in dp]
            st['pairs'] += len(vals)
            if all(v >= third for v in vals): st['L1_all_ge'] += 1
            else:
                st['L1_viol'] += 1
                if len(ex) < 6: ex.append(('L1', n, rel, D, U, [str(v) for v in vals]))
            mn = min(vals)
            if mn < minv: minv = mn
            if mn > maxofmin: maxofmin = mn
            if mn <= twothird: st['U1_min_le'] += 1
            else: st['U1_min_gt'] += 1
            # 对偶 U×U
            ul = [v for v in range(m) if (U >> v) & 1]
            up = [(x, y) for i, x in enumerate(ul) for y in ul[i+1:] if not (L[x][y] or L[y][x])]
            if up:
                uv = [my_pval(L, n, x, y) for (x, y) in up]
                if all(v <= twothird for v in uv): st['L1dual_all_le'] += 1
                else: st['L1dual_viol'] += 1
    print(f"[下界/上界] n={n} trials={trials}  案例 {st['cases']}")
    print(f"  L1 所有 D×D 对 p_P>=1/3 : {st['L1_all_ge']}   违反 : {st['L1_viol']}")
    print(f"  U1 带内 min p_P <= 2/3  : {st['U1_min_le']}   > 2/3 : {st['U1_min_gt']}")
    print(f"  L1' 所有 U×U 对 p_P<=2/3: {st['L1dual_all_le']}  违反 : {st['L1dual_viol']}")
    print(f"  带内 min p_P 范围: [{minv}, {maxofmin}]")
    for e in ex: print("  VIOL:", e)


if __name__ == '__main__':
    main()
