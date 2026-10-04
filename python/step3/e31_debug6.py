"""E31: debug why n=6 contributed 0 generic-corner fail cases in E30,
then measure the (mu, WV) joint distribution at the corner in fail cases (n=5+6).
"""
import random
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e11_corner.py").resolve(), encoding="utf-8").read().split("def main():")[0])

third = Fraction(1, 3)
random.seed(7)

# --- part 1: why 0 generic at n=6 ---
n = 6
npairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
cnt = dict(fail=0, generic=0, skipped=0)
reasons = {}
for _ in range(4000):
    mask = random.randrange(1 << len(npairs))
    rel = [npairs[k] for k in range(len(npairs)) if (mask >> k) & 1]
    L = closure(n, rel)
    if not is_poset(L, n): continue
    if antichain_width(L, n) < 2: continue
    z = random.randrange(n)
    Qz = [v for v in range(n) if v != z]
    Ql, _ = induced(L, n, Qz)
    if antichain_width(Ql, n - 1) > 2: continue
    r = analyze(L, n, z)
    if r['zpair']: continue
    cnt['fail'] += 1
    done = False
    for part in chain_partitions(Ql, n - 1):
        cc = corner_check(L, n, z, part)
        if cc is None: continue
        an = cc['an']
        ist, jst = cc['istar'], cc['jstar']
        s_, t_ = an['s'], an['t']
        if ist >= 1 and jst >= 1 and ist + 1 <= s_ and jst + 1 <= t_:
            cnt['generic'] += 1
            done = True
            break
        else:
            key = (ist == 0, jst == 0, ist + 1 > s_, jst + 1 > t_)
            reasons[key] = reasons.get(key, 0) + 1
            done = True
            break
    if not done:
        cnt['skipped'] += 1
print("part1 (n=6 debug):", cnt)
print("  edge-hit reasons (ist==0, jst==0, ist+1>s, jst+1>t):", reasons)
