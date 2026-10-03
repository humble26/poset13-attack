"""E21: how big can a single Delta-step get inside fail cases?
Delta_k(l) = P_L(b_l strictly between a_k and a_{k+1}) for chain-adjacent a's.
Record max Delta over fail cases (n=5 exhaustive, n=6 sampled), split by
whether the two neighbouring cells (k,l),(k+1,l) are both incomparable.
Also: same for the B-direction (row steps).
"""
import random
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e12_c1.py").resolve(), encoding="utf-8").read().split("def main():")[0])

third = Fraction(1, 3)


def deltas(less, n, z, part, an):
    """all column-step and row-step Delta values with neighbour comparability."""
    A, B = part
    s, t = an['s'], an['t']
    Q = [v for v in range(n) if v != z]
    out = []
    # column steps (A-direction): a_k, a_{k+1} chain-adjacent, b_l fixed
    for k in range(1, s):
        ak, ak1 = Q[A[k - 1]], Q[A[k]]
        for l in range(1, t + 1):
            bl = Q[B[l - 1]]
            if less[ak][bl] or less[bl][ak] or less[ak1][bl] or less[bl][ak1]:
                continue
            between = 0
            tot = count_ext(less, n)
            # count extensions with ak < bl < ak1 via DP over ranks? brute via pval:
            # P(ak<bl) - P(ak1<bl) is exactly Delta by the delta identity
            d = pval(less, n, ak, bl) - pval(less, n, ak1, bl)
            out.append(('col', k, l, d))
    for l in range(1, t):
        bl, bl1 = Q[B[l - 1]], Q[B[l]]
        for k in range(1, s + 1):
            ak = Q[A[k - 1]]
            if less[ak][bl] or less[bl][ak] or less[ak][bl1] or less[bl1][ak]:
                continue
            d = pval(less, n, bl, ak) - pval(less, n, bl1, ak)
            out.append(('row', k, l, d))
    return out


def main():
    random.seed(4242)
    mx = Fraction(0); mx_info = None
    hist = {}
    n_case = 0
    for n in (5,):
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
                    an = anatomy(L, n, z, part)
                    if an is None: continue
                    n_case += 1
                    for (tag, k, l, d) in deltas(L, n, z, part, an):
                        if d > mx:
                            mx = d; mx_info = (n, z, rel, part, tag, k, l)
                        key = d
                        hist[key] = hist.get(key, 0) + 1
    # n=6 sampled
    n = 6
    npairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    tried = 0
    while tried < 4000:
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
            an = anatomy(L, n, z, part)
            if an is None: continue
            n_case += 1
            for (tag, k, l, d) in deltas(L, n, z, part, an):
                if d > mx:
                    mx = d; mx_info = (n, z, rel, part, tag, k, l)
                hist[d] = hist.get(d, 0) + 1
    print("fail-partition instances:", n_case)
    print("max Delta:", mx, "at", mx_info)
    print("Delta histogram (top):", sorted(hist.items(), reverse=True)[:12])
    over = sum(v for k, v in hist.items() if k > third)
    print("steps with Delta > 1/3:", over)


if __name__ == '__main__':
    main()
