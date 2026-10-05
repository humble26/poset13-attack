"""E40: test sub-claims for Lemma 9.D.
S9D-1: for incomparable D-pairs: p_D in [1/3,2/3] => p_P in [1/3,2/3].
S9D-1strong: p_D in [1/3,2/3] => p_I in [1/3,2/3] for ALL box cuts I.
Counted over all incomparable D-pairs in fail+multi-cut cases.
"""
import random, collections
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e11_corner.py").resolve(), encoding="utf-8").read().split("def main():")[0])

third = Fraction(1, 3); twothird = Fraction(2, 3)
random.seed(8080)


def p_within(L, n, Ilist, x, y):
    Il, _ = induced(L, n, Ilist)
    pos = {v: k for k, v in enumerate(Ilist)}
    extsI = []
    def rec(placed, seq):
        if len(seq) == len(Ilist): extsI.append(list(seq)); return
        for v in Ilist:
            if (placed >> pos[v]) & 1: continue
            if all((placed >> pos[u]) & 1 for u in Ilist if L[u][v]):
                rec(placed | 1 << pos[v], seq + [v])
    rec(0, [])
    if not extsI: return None
    return Fraction(sum(1 for e in extsI if e.index(x) < e.index(y)), len(extsI))


def main():
    s1 = dict(ok=0, bad=0, badex=[])
    s1s = dict(ok=0, bad=0, badex=[])
    n_pair = 0
    for n in (5, 6):
        npairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        if n == 5:
            it = [(mask, z) for mask in range(1 << len(npairs)) for z in range(n)]
        else:
            it = [(random.randrange(1 << len(npairs)), random.randrange(n)) for _ in range(15000)]
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
                iD, jD = an['iD'], an['jD']
                if iD < 1 or jD < 1: continue
                if len(an['mu']) < 2: continue
                Q = [v for v in range(n) if v != z]
                A, B = part
                Dlist = sorted(set(A[:iD]) | set(B[:jD]))
                for k in range(1, iD + 1):
                    for l in range(1, jD + 1):
                        x, y = Q[A[k - 1]], Q[B[l - 1]]
                        if L[x][y] or L[y][x]: continue
                        pD = p_within(L, n, Dlist, x, y)
                        if not (third <= pD <= twothird): continue
                        n_pair += 1
                        pP = pval(L, n, x, y)
                        if third <= pP <= twothird: s1['ok'] += 1
                        else:
                            s1['bad'] += 1
                            if len(s1['badex']) < 3:
                                s1['badex'].append((n, z, (k, l), str(pD), str(pP)))
                        # strong: all box cuts
                        allin = True
                        for (i, j) in an['mu']:
                            Ilist = sorted(set(A[:i]) | set(B[:j]))
                            pI = p_within(L, n, Ilist, x, y)
                            if pI is None or not (third <= pI <= twothird):
                                allin = False; break
                        if allin: s1s['ok'] += 1
                        else:
                            s1s['bad'] += 1
                            if len(s1s['badex']) < 3:
                                s1s['badex'].append((n, z, (k, l), str(pD), str(pP)))
    print("S9D-1 pairs tested:", n_pair, "| ok:", s1['ok'], "bad:", s1['bad'])
    for f in s1['badex']: print("   S9D-1 FAIL:", f)
    print("S9D-1strong ok:", s1s['ok'], "bad:", s1s['bad'])
    for f in s1s['badex']: print("   STRONG FAIL:", f)


if __name__ == '__main__':
    main()
