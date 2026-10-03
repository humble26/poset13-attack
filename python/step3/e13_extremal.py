"""E13: find and dissect extremal (min-delta) cases of the target class at n=6.

For each (P,z) in the target class (width(Q)<=2, width(P)>=2, P not chain... any):
compute delta(P); track the minima. Then dissect the extremal cases:
 - are they zpair-fail? which pair attains delta?
 - corner structure (i*, j*), heavy layer, balanced same-side pairs.
"""
import itertools, sys
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main()")[0])
exec(open((Path(__file__).parent / "e11_corner.py").resolve(), encoding="utf-8").read().split("def main():")[0])

third = Fraction(1, 3)


def delta_of(less, n):
    best = Fraction(0)
    arg = None
    for x, y in itertools.combinations(range(n), 2):
        if less[x][y] or less[y][x]: continue
        p = pval(less, n, x, y)
        m = min(p, 1 - p)
        if m > best:
            best = m; arg = (x, y, p)
    return best, arg


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    top = []
    for n in range(5, nmax + 1):
        npairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        all_cases = []
        for mask in range(1 << len(npairs)):
            rel = [npairs[k] for k in range(len(npairs)) if (mask >> k) & 1]
            L = closure(n, rel)
            if not is_poset(L, n): continue
            if antichain_width(L, n) < 2: continue
            for z in range(n):
                Qz = [v for v in range(n) if v != z]
                Ql, _ = induced(L, n, Qz)
                if antichain_width(Ql, n - 1) > 2: continue
                d, arg = delta_of(L, n)
                all_cases.append((d, tuple(rel), z, arg))
        mn = min(c[0] for c in all_cases)
        hard = [c for c in all_cases if c[0] == mn]
        print("n=%d target-class cases=%d  min delta=%s  #hard=%d" % (n, len(all_cases), mn, len(hard)))
        for (d, rel, z, arg) in hard[:4]:
            print("  hard: delta=%s arg(pair,p)=(%s) rel=%s z=%d" % (d, arg, rel, z))
            L = closure(n, rel)
            r = analyze(L, n, z)
            print("    zpair-fail:", not r['zpair'], " qpair-bal:", r['qpair'], " delta=%s best_z=%s" % (r['delta'], r['best_z']))
            Qz = [v for v in range(n) if v != z]
            Ql, _ = induced(L, n, Qz)
            for part in chain_partitions(Ql, n - 1):
                an = anatomy(L, n, z, part)
                if an is None: continue
                print("    A=%s B=%s iD=%d iU=%d jD=%d jU=%d" % (part[0], part[1], an['iD'], an['iU'], an['jD'], an['jU']))
                print("    MAv=%s" % {k: str(v) for k, v in sorted(an['MAv'].items())})
                print("    MBv=%s" % {k: str(v) for k, v in sorted(an['MBv'].items())})
                print("    mu=%s" % {k: str(v) for k, v in sorted(an['mu'].items())})
                print("    bal_pairs=%s" % [(cxA, cyB, str(pp)) for (_, _, pp, cxA, cyB, *_ ) in an['bal_pairs']])
                break


if __name__ == '__main__':
    main()
