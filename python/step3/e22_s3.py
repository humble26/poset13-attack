"""E22: verify claim S3 (turn structure / landing cell).
In a fail case (all z-pairs unbalanced), for EVERY chain partition:
there is an incomparable cross pair (a_k, b_l) with p in [1/3,2/3] such that
among its existing 4-neighbours in the p-matrix (comparable cells count as p=1/0)
there is one with p > 2/3 AND one with p < 1/3.
Run: python e22_s3.py  (n=5 exhaustive, n=6 and n=7 sampled)
"""
import random
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e12_c1.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e20_walk.py").resolve(), encoding="utf-8").read().split("def partA()")[0])

third = Fraction(1, 3); twothird = Fraction(2, 3)


def cellp(pm, k, l, less, n, z, part):
    """p-value of cell (k,l); comparable cells give 1 or 0."""
    if (k, l) not in pm: return None
    rel, p = pm[(k, l)]
    if rel == 'xy': return Fraction(1)
    if rel == 'yx': return Fraction(0)
    return p


def s3_holds(less, n, z, part, an, pm):
    A, B = part
    s, t = an['s'], an['t']
    for k in range(1, s + 1):
        for l in range(1, t + 1):
            if (k, l) not in pm: continue
            rel, p = pm[(k, l)]
            if rel != 'inc' or not (third <= p <= twothird): continue
            highs = lows = 0
            for (kk, ll) in ((k, l + 1), (k - 1, l), (k, l - 1), (k + 1, l)):
                q = cellp(pm, kk, ll, less, n, z, part)
                if q is None: continue
                if q > twothird: highs += 1
                elif q < third: lows += 1
            if highs >= 1 and lows >= 1:
                return True, (k, l)
    return False, None


def main():
    agg = dict(n5=dict(cases=0, ok=0, bad=[]), n6=dict(cases=0, ok=0, bad=[]),
               n7=dict(cases=0, ok=0, bad=[]))
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
            ok, where = s3_holds(L, n, z, part, an, pm)
            if ok: bucket['ok'] += 1
            else: bucket['bad'].append((n, z, part, where))

    n = 5
    npairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for mask in range(1 << len(npairs)):
        rel = [npairs[k] for k in range(len(npairs)) if (mask >> k) & 1]
        L = closure(n, rel)
        if not is_poset(L, n): continue
        if antichain_width(L, n) < 2: continue
        for z in range(n):
            check(L, n, z, agg['n5'])
    for n, nsamp in ((6, 6000), (7, 4000)):
        npairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        tried = 0
        while tried < nsamp:
            tried += 1
            mask = random.randrange(1 << len(npairs))
            rel = [npairs[k] for k in range(len(npairs)) if (mask >> k) & 1]
            L = closure(n, rel)
            if not is_poset(L, n): continue
            if antichain_width(L, n) < 2: continue
            z = random.randrange(n)
            check(L, n, z, agg['n%d' % n])
    for key in ('n5', 'n6', 'n7'):
        b = agg[key]
        print(key, "cases:", b['cases'], "ok:", b['ok'], "bad:", len(b['bad']))
        for bb in b['bad'][:4]:
            print("   BAD:", bb)


if __name__ == '__main__':
    main()
