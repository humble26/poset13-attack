"""E17: consolidated verification of session-3 claims on fail cases
(P = Q+z, width(Q)<=2, all z-pairs unbalanced):
  G1 heavy layer: some box column or row has mass > 1/3            [theorem track]
  G2 squeeze: at theta-jump columns c (a_{c+1} in window, colmass>1/3,
     M_A(c)<1/3): 2/3 - colmass[c] < M_A(c) < 1/3                   [theorem track]
  S1 refined conservation: some balanced Q-pair lies in a SAME-SIDE block:
     (a_k,b_l) inc with (k<=i* and l<=j*) or (k>=i*+1 and l>=j*+1)   [target law]
  S2 mixed corner orientation: p(a_{i*},b_{j*+1})>2/3, p(a_{i*+1},b_{j*})<1/3
     when those cells exist and are incomparable                     [target law]
Run: python e17_final.py nmax [nsamp6]
"""
import random, sys
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e12_c1.py").resolve(), encoding="utf-8").read().split("def main():")[0])

third = Fraction(1, 3); twothird = Fraction(2, 3)
agg = dict(failparts=0, G1_ok=0, G1_bad=0, G2_ok=0, G2_bad=0,
           S1_ok=0, S1_none=0, S1_outside=0, S2a=0, S2a_ok=0, S2b=0, S2b_ok=0)
bad = []

def check_case(L, n, z):
    Qz = [v for v in range(n) if v != z]
    Ql, _ = induced(L, n, Qz)
    if antichain_width(Ql, n - 1) > 2: return
    r = analyze(L, n, z)
    if r['zpair']: return
    for part in chain_partitions(Ql, n - 1):
        cc = corner_check(L, n, z, part)
        if cc is None: continue
        an = cc['an']; agg['failparts'] += 1
        istar, jstar = cc['istar'], cc['jstar']
        # G1
        heavy = any(w > third for w in an['colmass'].values()) or \
                any(w > third for w in an['rowmass'].values())
        if heavy: agg['G1_ok'] += 1
        else: agg['G1_bad'] += 1; bad.append(('G1', n, z, part))
        # G2 squeeze
        for c, m0, cm in an['squeezes']:
            if twothird - cm < m0 < third: agg['G2_ok'] += 1
            else: agg['G2_bad'] += 1; bad.append(('G2', n, z, part, c))
        # S1
        A, B = part; s, t = an['s'], an['t']
        found_in = False; found_any = False
        for (ba, bb, pp, cxA, cyB, lx, hx, ly, hy, k, l) in an['bal_pairs']:
            found_any = True
            if (k <= istar and l <= jstar) or (k >= istar + 1 and l >= jstar + 1):
                found_in = True; break
        if found_in: agg['S1_ok'] += 1
        elif found_any:
            agg['S1_outside'] += 1
            bad.append(('S1-outside', n, z, part, istar, jstar,
                        [(k, l, str(pp)) for (_, _, pp, *_ , k, l) in an['bal_pairs']]))
        else:
            agg['S1_none'] += 1; bad.append(('S1-none', n, z, part))
        # S2
        Q = [v for v in range(n) if v != z]
        cells = cc['cells']
        c01 = cells.get((0, 1)); c10 = cells.get((1, 0))
        if c01 and c01[0] == 'inc':
            agg['S2a'] += 1
            if c01[1] > twothird: agg['S2a_ok'] += 1
            else: bad.append(('S2a', n, z, part, str(c01[1])))
        if c10 and c10[0] == 'inc':
            agg['S2b'] += 1
            if c10[1] < third: agg['S2b_ok'] += 1
            else: bad.append(('S2b', n, z, part, str(c10[1])))


nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
nsamp = int(sys.argv[2]) if len(sys.argv) > 2 else 0
for n in range(4, nmax + 1):
    npairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for mask in range(1 << len(npairs)):
        rel = [npairs[k] for k in range(len(npairs)) if (mask >> k) & 1]
        L = closure(n, rel)
        if not is_poset(L, n): continue
        if antichain_width(L, n) < 2: continue
        for z in range(n):
            check_case(L, n, z)
if nsamp:
    random.seed(2026)
    n = nmax + 1 if nsamp else 6
    n = 6
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
        check_case(L, n, z)
print(agg)
print("violations:", len(bad))
for b in bad[:10]:
    print(b)
