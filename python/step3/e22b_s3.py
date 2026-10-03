"""E22b: S3 turn-structure verification with full aggregates (n=5 exhaustive, n=6 sampled)."""
import random
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e11_corner.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e20_walk.py").resolve(), encoding="utf-8").read().split("def partA()")[0])

third = Fraction(1, 3); twothird = Fraction(2, 3)


def cellp(pm, k, l):
    if (k, l) not in pm: return None
    rel, p = pm[(k, l)]
    if rel == 'xy': return Fraction(1)
    if rel == 'yx': return Fraction(0)
    return p


def s3_holds(pm, an):
    s, t = an['s'], an['t']
    for k in range(1, s + 1):
        for l in range(1, t + 1):
            if (k, l) not in pm: continue
            rel, p = pm[(k, l)]
            if rel != 'inc' or not (third <= p <= twothird): continue
            highs = lows = 0
            for (kk, ll) in ((k, l + 1), (k - 1, l), (k, l - 1), (k + 1, l)):
                q = cellp(pm, kk, ll)
                if q is None: continue
                if q > twothird: highs += 1
                elif q < third: lows += 1
            if highs >= 1 and lows >= 1:
                return True, (k, l)
    return False, None


def pmatrix(less, n, z, part, an):
    A, B = part
    s, t = an['s'], an['t']
    Q = [v for v in range(n) if v != z]
    out = {}
    for a1 in range(s):
        for b1 in range(t):
            gx, gy = Q[A[a1]], Q[B[b1]]
            k, l = a1 + 1, b1 + 1
            if less[gx][gy]: rel = 'xy'
            elif less[gy][gx]: rel = 'yx'
            else: rel = 'inc'
            p = pval(less, n, gx, gy) if rel == 'inc' else None
            out[(k, l)] = (rel, p)
    return out


random.seed(314)


def check(L, n, z, bucket):
    Qz = [v for v in range(n) if v != z]
    Ql, _ = induced(L, n, Qz)
    if antichain_width(Ql, n - 1) > 2: return
    r = analyze(L, n, z)
    if r['zpair']: return
    for part in chain_partitions(Ql, n - 1):
        cc = corner_check(L, n, z, part)
        if cc is None: continue
        an = cc['an']
        pm = pmatrix(L, n, z, part, an)
        bucket['cases'] += 1
        ok, where = s3_holds(pm, an)
        if ok: bucket['ok'] += 1
        else: bucket['bad'].append((n, z, part))


b5 = dict(cases=0, ok=0, bad=[]); b6 = dict(cases=0, ok=0, bad=[])
n = 5
npairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
for mask in range(1 << len(npairs)):
    rel = [npairs[k] for k in range(len(npairs)) if (mask >> k) & 1]
    L = closure(n, rel)
    if not is_poset(L, n): continue
    if antichain_width(L, n) < 2: continue
    for z in range(n):
        check(L, n, z, b5)
n = 6
npairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
tried = 0
while tried < 3000:
    tried += 1
    mask = random.randrange(1 << len(npairs))
    rel = [npairs[k] for k in range(len(npairs)) if (mask >> k) & 1]
    L = closure(n, rel)
    if not is_poset(L, n): continue
    if antichain_width(L, n) < 2: continue
    z = random.randrange(n)
    check(L, n, z, b6)
print('n5 exhaustive: cases', b5['cases'], 'ok', b5['ok'], 'bad', len(b5['bad']))
for x in b5['bad'][:3]: print('  BAD', x)
print('n6 sampled: cases', b6['cases'], 'ok', b6['ok'], 'bad', len(b6['bad']))
for x in b6['bad'][:3]: print('  BAD', x)
