"""dump：p_D 最接近 1/2 的 D-见证对**不存活**的案例（T3 失效），
看 T1 的存活对与不存活对的差别。"""
import sys
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e11_corner.py").resolve(), encoding="utf-8").read().split("def main():")[0])

third = Fraction(1, 3); twothird = Fraction(2, 3)


def p_within(less, n, Ilist, x, y):
    Il, idx = induced(less, n, Ilist); m = len(Ilist); eI = count_ext(Il, m)
    L2 = closure(m, [(i, j) for i in range(m) for j in range(m) if Il[i][j]] + [(idx[x], idx[y])])
    return Fraction(count_ext(L2, m), eI)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    shown = 0
    for mask in range(1 << (n * (n - 1) // 2)):
        npairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        rel = [npairs[k] for k in range(len(npairs)) if (mask >> k) & 1]
        L = closure(n, rel)
        if not is_poset(L, n) or antichain_width(L, n) < 2: continue
        for z in range(n):
            Qz = [v for v in range(n) if v != z]; Ql, _ = induced(L, n, Qz)
            if antichain_width(Ql, n - 1) > 2: continue
            r = analyze(L, n, z)
            if r['zpair']: continue
            for part in chain_partitions(Ql, n - 1):
                cc = corner_check(L, n, z, part)
                if cc is None: continue
                an = cc['an']; iD, jD = an['iD'], an['jD']
                if len(an['mu']) < 2 or iD < 1 or jD < 1: continue
                A, B = part; Q = [v for v in range(n) if v != z]
                band = [(k, l) for k in range(1, iD+1) for l in range(1, jD+1)
                        if not (L[Q[A[k-1]]][Q[B[l-1]]] or L[Q[B[l-1]]][Q[A[k-1]]])]
                if not band: continue
                Imin = [Q[a] for a in A[:iD]] + [Q[b] for b in B[:jD]]
                rows = [(k, l, p_within(L, n, Imin, Q[A[k-1]], Q[B[l-1]]),
                         pval(L, n, Q[A[k-1]], Q[B[l-1]])) for (k, l) in band]
                best = min(rows, key=lambda t: abs(t[2] - Fraction(1, 2)))
                if third <= best[3] <= twothird: continue
                shown += 1
                if shown > 4: continue
                print("=" * 66)
                print(f"n={n} z={z} A={A} B={B} iD={iD} jD={jD} box={[(k, str(v)) for k, v in sorted(an['mu'].items())]}")
                for (k, l, pd, pp) in rows:
                    print(f"  (a{k},b{l}) p_D={pd}  p_P={pp}  p_D∈区间={third<=pd<=twothird}  p_P∈区间={third<=pp<=twothird}")
    print(f"\n共 {shown} 个 T3 失效案例")


if __name__ == '__main__':
    main()
