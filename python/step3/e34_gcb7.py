"""E34: capture n=7 GCB counterexamples (both corner cells unbalanced) for the record."""
import random
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e11_corner.py").resolve(), encoding="utf-8").read().split("def main():")[0])

third = Fraction(1, 3); twothird = Fraction(2, 3)
random.seed(2024)
n = 7
npairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
found = 0
tried = 0
while tried < 50000 and found < 2:
    tried += 1
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
    for part in chain_partitions(Ql, n - 1):
        cc = corner_check(L, n, z, part)
        if cc is None: continue
        an = cc['an']
        ist, jst = cc['istar'], cc['jstar']
        s_, t_ = an['s'], an['t']
        if ist == 0 or jst == 0 or ist + 1 > s_ or jst + 1 > t_: continue
        Q = [v for v in range(n) if v != z]
        A, B = part
        cells = {}
        okgen = True
        for (nm, kk, ll) in (('m1', ist, jst + 1), ('m2', ist + 1, jst),
                             ('c1', ist, jst), ('c2', ist + 1, jst + 1)):
            x, y = Q[A[kk - 1]], Q[B[ll - 1]]
            if L[x][y] or L[y][x]: okgen = False; break
            cells[nm] = pval(L, n, x, y)
        if not okgen: continue
        if third <= cells['c1'] <= twothird or third <= cells['c2'] <= twothird: continue
        found += 1
        print('n=7 GCB-counterexample #%d: rel=%s z=%d' % (found, rel, z))
        print('  A=%s B=%s ist=%d jst=%d' % (A, B, ist, jst))
        print('  m1=%s m2=%s c1=%s c2=%s' % (cells['m1'], cells['m2'], cells['c1'], cells['c2']))
        print('  mu=%s' % {k: str(v) for k, v in sorted(an['mu'].items())})
        print('  MAv=%s' % {k: str(v) for k, v in sorted(an['MAv'].items())})
        print('  MBv=%s' % {k: str(v) for k, v in sorted(an['MBv'].items())})
        # where is the balanced pair? (S1 check)
        bal = []
        for a1 in range(an['s']):
            for b1 in range(an['t']):
                x, y = Q[A[a1]], Q[B[b1]]
                if L[x][y] or L[y][x]: continue
                p = pval(L, n, x, y)
                if third <= p <= twothird:
                    bal.append((a1 + 1, b1 + 1, str(p)))
        print('  balanced Q-pairs (k,l,p):', bal)
print('done, found', found)
