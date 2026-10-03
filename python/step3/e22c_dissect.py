"""E22c: dissect S3-BAD cases — print the full p-matrix of the first BADs."""
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e11_corner.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e22b_s3.py").resolve(), encoding="utf-8").read().split("random.seed(314)")[0].replace("def pmatrix", "def pmatrixOLD"))
exec(open((Path(__file__).parent / "e22b_s3.py").resolve(), encoding="utf-8").read().split("def pmatrix(less")[1].split("random.seed")[0].join(["def pmatrix(less", ""]))

third = Fraction(1, 3); twothird = Fraction(2, 3)
n = 5
npairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
found = 0
for mask in range(1 << len(npairs)):
    rel = [npairs[k] for k in range(len(npairs)) if (mask >> k) & 1]
    L = closure(n, rel)
    if not is_poset(L, n): continue
    if antichain_width(L, n) < 2: continue
    z = 0
    Qz = [v for v in range(n) if v != z]
    Ql, _ = induced(L, n, Qz)
    if antichain_width(Ql, n - 1) > 2: continue
    r = analyze(L, n, z)
    if r['zpair']: continue
    for part in chain_partitions(Ql, n - 1):
        if tuple(part[0]) != (0, 1, 2) or tuple(part[1]) != (3,): continue
        cc = corner_check(L, n, z, part)
        an = cc['an']
        pm = pmatrix(L, n, z, part, an)
        ok, where = s3_holds(pm, an)
        if ok: continue
        print('rel=', rel, 'istar=', cc['istar'], 'jstar=', cc['jstar'])
        print('MAv=', {k: str(v) for k, v in sorted(an['MAv'].items())})
        print('MBv=', {k: str(v) for k, v in sorted(an['MBv'].items())})
        print('mu=', {k: str(v) for k, v in sorted(an['mu'].items())})
        for k in range(1, an['s'] + 1):
            row = []
            for l in range(1, an['t'] + 1):
                rl, p = pm[(k, l)]
                row.append((rl, str(p) if p is not None else ''))
            print(' row', k, row)
        found += 1
    if found >= 3: break
