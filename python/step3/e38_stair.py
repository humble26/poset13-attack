"""E38: residual-case sharp candidate: the staircase-corner pairs
S1 = (a_{iD}, b_{jD+1}) and S2 = (b_{jD}, a_{iD+1}) — is at least one
(existing & incomparable) balanced in residual fail cases?
"""
import random
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e11_corner.py").resolve(), encoding="utf-8").read().split("def main():")[0])

third = Fraction(1, 3); twothird = Fraction(2, 3)
random.seed(424242)

stats = dict(res=0, cand_exist=0, cand_ok=0, cand_bad=0, nocand=0)
fails = []
for n in (5, 6):
    npairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    if n == 5:
        it = [(mask, z) for mask in range(1 << len(npairs)) for z in range(n)]
    else:
        it = [(random.randrange(1 << len(npairs)), random.randrange(n)) for _ in range(25000)]
    for mask, z in it:
        rel = [npairs[k] for k in range(len(npairs)) if (mask >> k) & 1]
        L = closure(n, rel)
        if not is_poset(L, n): continue
        if antichain_width(L, n) < 2: continue
        Qz = [v for v in range(n) if v != z]
        Ql, _ = induced(L, n, Qz)
        if antichain_width(Ql, n - 1) > 2: continue
        r = analyze(L, n, z)
        if r['zpair']: continue
        for part in chain_partitions(Ql, n - 1):
            cc = corner_check(L, n, z, part)
            if cc is None: continue
            an = cc['an']
            iD, jD, iU, jU = an['iD'], an['jD'], an['iU'], an['jU']
            if len(an['mu']) < 2: continue
            s_, t_ = an['s'], an['t']
            Q = [v for v in range(n) if v != z]
            A, B = part
            dinc = any(not (L[Q[A[k-1]]][Q[B[l-1]]] or L[Q[B[l-1]]][Q[A[k-1]]])
                       for k in range(1, iD + 1) for l in range(1, jD + 1)) if (iD >= 1 and jD >= 1) else False
            uinc = any(not (L[Q[A[k-1]]][Q[B[l-1]]] or L[Q[B[l-1]]][Q[A[k-1]]])
                       for k in range(iU + 1, s_ + 1) for l in range(jU + 1, t_ + 1)) if (s_ - iU >= 1 and t_ - jU >= 1) else False
            if dinc or uinc: continue  # covered by 9.8
            stats['res'] += 1
            cands = []
            if iD >= 1 and jD + 1 <= t_:
                x, y = Q[A[iD - 1]], Q[B[jD]]
                if not (L[x][y] or L[y][x]): cands.append(pval(L, n, x, y))
            if jD >= 1 and iD + 1 <= s_:
                x, y = Q[B[jD - 1]], Q[A[iD]]
                if not (L[x][y] or L[y][x]): cands.append(pval(L, n, x, y))
            if not cands:
                stats['nocand'] += 1
                continue
            stats['cand_exist'] += 1
            if any(third <= p_ <= twothird for p_ in cands): stats['cand_ok'] += 1
            else:
                stats['cand_bad'] += 1
                if len(fails) < 4:
                    fails.append((n, z, part, [str(p_) for p_ in cands],
                                  {k: str(v) for k, v in sorted(an['mu'].items())},
                                  {k: str(v) for k, v in sorted(an['MAv'].items())},
                                  {k: str(v) for k, v in sorted(an['MBv'].items())}))
print(stats)
for f in fails: print("  FAIL:", f)
