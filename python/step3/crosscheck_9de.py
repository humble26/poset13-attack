"""交叉验证：独立实现 my_poset vs 仓库 analyze/anatomy/pval。
对同一批实例，比较 (a) fail 判定 (b) 全部不可比对对的 p 值 (c) box 权重和。"""
import random, sys
from fractions import Fraction
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
# 先加载仓库代码
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e11_corner.py").resolve(), encoding="utf-8").read().split("def main():")[0])
# 再导入独立实现（不覆盖仓库名）
from my_poset import (my_close, my_is_poset, my_count_ext, my_width,
                      my_pval, my_rand_poset, my_induced)

third = Fraction(1, 3); twothird = Fraction(2, 3)
rng = random.Random(31337)
n = 6; m = n - 1
checked = 0; bad_p = 0; bad_fail = 0; bad_box = 0
for _ in range(60000):
    Q = my_rand_poset(m, rng.uniform(0.25, 0.6), rng)
    if not my_is_poset(Q, m) or my_width(Q, m) > 2: continue
    DS = [S for S in range(1 << m) if all(not ((S >> i) & 1) or all((S >> j) & 1 for j in range(m) if Q[j][i]) for i in range(m))]
    US = [S for S in range(1 << m) if all(not ((S >> i) & 1) or all((S >> j) & 1 for j in range(m) if Q[i][j]) for i in range(m))]
    for _ in range(2):
        D = rng.choice(DS); U = rng.choice(US)
        if D & U: continue
        rel = [(i, j) for i in range(m) for j in range(m) if Q[i][j]]
        rel += [(d, m) for d in range(m) if (D >> d) & 1]
        rel += [(m, u) for u in range(m) if (U >> u) & 1]
        # 独立实现的 L
        Lmy = my_close(n, rel)
        if not my_is_poset(Lmy, n) or my_width(Lmy, n) < 2: continue
        # 仓库实现的 L（同 rel，用仓库 closure）
        Lrepo = closure(n, rel)
        if not is_poset(Lrepo, n): continue
        # (a) 全对 p 值一致
        for x in range(n):
            for y in range(n):
                if x == y: continue
                if Lmy[x][y] or Lmy[y][x]: continue
                a = my_pval(Lmy, n, x, y); b = pval(Lrepo, n, x, y)
                if a != b:
                    bad_p += 1
                    if bad_p <= 3: print("P-MISMATCH", rel, x, y, a, b)
        # (b) fail 判定一致
        fmy = not any(third <= my_pval(Lmy, n, m, v) <= twothird for v in range(m) if not (Lmy[v][m] or Lmy[m][v]))
        r = analyze(Lrepo, n, m)
        frepo = not r['zpair']
        if fmy != frepo:
            bad_fail += 1
            if bad_fail <= 3: print("FAIL-MISMATCH", rel, fmy, frepo)
        # (c) box 权重和 = 1（corner_check 内部断言只在 fail 情形成立）
        if frepo:
            Ql, _ = induced(Lrepo, n, [v for v in range(n) if v != m])
            for part in chain_partitions(Ql, m):
                cc = corner_check(Lrepo, n, m, part)
                if cc is None: continue
                an = cc['an']
                if sum(an['mu'].values()) != 1:
                    bad_box += 1
                break
        checked += 1
        if checked >= 150: break
    if checked >= 150: break
print(f"交叉验证：实例 {checked}，p 值不一致 {bad_p}，fail 不一致 {bad_fail}，box 非概率 {bad_box}")
print("结论：" + ("全部一致 ✅" if bad_p == 0 and bad_fail == 0 and bad_box == 0 else "存在不一致 ❌"))
