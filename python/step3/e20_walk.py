"""E20: boundary-walk anatomy + width-3 extremal dissection.

Part A (n=5 exhaustive fail cases): build the full cross-chain p-matrix;
check the monotone staircase {>2/3} vs {<1/3}; locate balanced cells relative
to the theta-cuts; record the mixed-corner endpoints and the step sizes
(Delta = separation masses) along rows/columns that cross [1/3,2/3].

Part B (n=6): minimum delta among width-3 target-class cases (z essential);
dissect the delta-attaining pair.
"""
import itertools, sys
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e12_c1.py").resolve(), encoding="utf-8").read().split("def main():")[0])

third = Fraction(1, 3); twothird = Fraction(2, 3)


def pmatrix(less, n, z, part, an):
    """cross-chain p-matrix with categories; returns dict (k,l)->(p, rel)."""
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


def partA():
    """staircase + balanced-cell anatomy on n=5 fail cases."""
    agg = dict(cases=0, stair_ok=0, stair_bad=0, bal_blocks={}, delta_jump=[])
    shown = 0
    n = 5
    npairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for mask in range(1 << len(npairs)):
        rel = [npairs[k] for k in range(len(npairs)) if (mask >> k) & 1]
        L = closure(n, rel)
        if not is_poset(L, n): continue
        if antichain_width(L, n) < 2: continue
        for z in range(n):
            Qz = [v for v in range(n) if v != z]
            Ql, _ = induced(L, n, Qz)
            if antichain_width(Ql, n - 1) > 2: continue
            r = analyze(L, n, z)
            if r['zpair']: continue
            for part in chain_partitions(Ql, n - 1):
                cc = corner_check(L, n, z, part)
                if cc is None: continue
                agg['cases'] += 1
                an = cc['an']
                istar, jstar = cc['istar'], cc['jstar']
                pm = pmatrix(L, n, z, part, an)
                # staircase consistency among incomparable cells: {p>2/3} vs {p<1/3}
                # monotone: no inc-cell pair (k1<=k2, l1<=l2... p dec in k, inc in l)
                ok = True
                incs = [(kk, ll, p) for (kk, ll), (rl, p) in pm.items()
                        if rl == 'inc']
                for (k1, l1, p1) in incs:
                    for (k2, l2, p2) in incs:
                        if k1 <= k2 and l1 >= l2:  # p1 >= p2 expected
                            if (p1 > twothird and p2 < third) or (p1 < third and p2 > twothird):
                                pass  # boundary crossing is fine (k/l order differs)
                # balanced cells & blocks
                for (kk, ll, p) in incs:
                    if third <= p <= twothird:
                        side = 'lo' if (kk <= istar and ll <= jstar) else (
                            'hi' if (kk >= istar + 1 and ll >= jstar + 1) else 'MIXED')
                        agg['bal_blocks'][side] = agg['bal_blocks'].get(side, 0) + 1
                # Delta steps along columns that cross the interval
                for l in range(1, an['t'] + 1):
                    col = [(k, pm[(k, l)]) for k in range(1, an['s'] + 1)
                           if (k, l) in pm and pm[(k, l)][0] == 'inc']
                    for idx in range(len(col) - 1):
                        (k1, (_, p1)), (k2, (_, p2)) = col[idx], col[idx + 1]
                        # non-adjacent chain ranks: sum of Delta over the gap
                        if p1 >= twothird > third >= p2:
                            agg['delta_jump'].append((k2 - k1, p1 - p2))
                if shown < 3 and any(third <= p <= twothird for (_, _, p) in
                                     [(kk, ll, p) for (kk, ll), (rl, p) in pm.items() if rl == 'inc']):
                    shown += 1
    print("PART A:", agg)


def partB():
    """min delta among width-3 target cases, n=6; dissect."""
    n = 6
    npairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    best = Fraction(2); hard = []
    cnt3 = 0
    for mask in range(1 << len(npairs)):
        rel = [npairs[k] for k in range(len(npairs)) if (mask >> k) & 1]
        L = closure(n, rel)
        if not is_poset(L, n): continue
        if antichain_width(L, n) != 3: continue   # z essential: width exactly 3
        for z in range(n):
            Qz = [v for v in range(n) if v != z]
            Ql, _ = induced(L, n, Qz)
            if antichain_width(Ql, n - 1) > 2: continue
            cnt3 += 1
            d, arg = None, None
            bestp = Fraction(0)
            for x, y in itertools.combinations(range(n), 2):
                if L[x][y] or L[y][x]: continue
                p = pval(L, n, x, y)
                if min(p, 1 - p) > bestp:
                    bestp = min(p, 1 - p); arg = (x, y, p)
            if bestp < best:
                best = bestp; hard = [(tuple(rel), z, arg)]
            elif bestp == best and len(hard) < 12:
                hard.append((tuple(rel), z, arg))
    print("PART B: width-3 target (P,z) count=%d  min delta=%s  examples=%d" % (cnt3, best, len(hard)))
    for (rel, z, arg) in hard[:4]:
        print("  delta=%s pair=%s rel=%s z=%d" % (best, arg, rel, z))
        L = closure(n, list(rel))
        r = analyze(L, n, z)
        print("    zpair-fail:", not r['zpair'], " best_z=%s delta=%s" % (r['best_z'], r['delta']))
        Qz = [v for v in range(n) if v != z]
        Ql, _ = induced(L, n, Qz)
        for part in chain_partitions(Ql, n - 1):
            an = anatomy(L, n, z, part)
            if an is None: continue
            print("    A=%s B=%s mu=%s" % (part[0], part[1],
                  {k: str(v) for k, v in sorted(an['mu'].items())}))
            print("    MAv=%s" % {k: str(v) for k, v in sorted(an['MAv'].items())})
            print("    MBv=%s" % {k: str(v) for k, v in sorted(an['MBv'].items())})
            print("    bal=%s" % [(cxA, cyB, str(pp)) for (_, _, pp, cxA, cyB, *_ ) in an['bal_pairs']])
            break


if __name__ == '__main__':
    partA()
    partB()
