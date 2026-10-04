"""诊断：n=6 fail+multicut 且 DxD 有不可比对的案例，
列出每个 DxD 对的 p_P，以及盒最小点(=D)、最大点处的 p_I。
目标：看"哪个对平衡"以及平衡从何而来。"""
import sys, collections
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e11_corner.py").resolve(), encoding="utf-8").read().split("def main():")[0])

third = Fraction(1, 3); twothird = Fraction(2, 3)


def p_within(less, n, Ilist, x, y):
    Il, idx = induced(less, n, Ilist)
    m = len(Ilist)
    eI = count_ext(Il, m)
    L2 = closure(m, [(i, j) for i in range(m) for j in range(m) if Il[i][j]] + [(idx[x], idx[y])])
    return Fraction(count_ext(L2, m), eI)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    hist = collections.Counter()
    stats = dict(cases=0, npairs=collections.Counter(), ok=0, bad=0)
    pD_in = 0; pD_tot = 0
    pDmax_in = 0
    shown = 0
    for mask in range(1 << (n * (n - 1) // 2)):
        npairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        rel = [npairs[k] for k in range(len(npairs)) if (mask >> k) & 1]
        L = closure(n, rel)
        if not is_poset(L, n) or antichain_width(L, n) < 2: continue
        for z in range(n):
            Qz = [v for v in range(n) if v != z]
            Ql, _ = induced(L, n, Qz)
            if antichain_width(Ql, n - 1) > 2: continue
            r = analyze(L, n, z)
            if r['zpair']: continue
            for part in chain_partitions(Ql, n - 1):
                cc = corner_check(L, n, z, part)
                if cc is None: continue
                an = cc['an']
                iD, jD, iU, jU = an['iD'], an['jD'], an['iU'], an['jU']
                if len(an['mu']) < 2: continue
                A, B = part; s_, t_ = an['s'], an['t']
                Q = [v for v in range(n) if v != z]
                if iD < 1 or jD < 1: continue
                dp = [(k, l, Q[A[k-1]], Q[B[l-1]]) for k in range(1, iD+1) for l in range(1, jD+1)
                      if not (L[Q[A[k-1]]][Q[B[l-1]]] or L[Q[B[l-1]]][Q[A[k-1]]])]
                if not dp: continue
                stats['cases'] += 1
                stats['npairs'][len(dp)] += 1
                pv = [(k, l, pval(L, n, x, y)) for (k, l, x, y) in dp]
                if any(third <= p <= twothird for (_, _, p) in pv): stats['ok'] += 1
                else: stats['bad'] += 1
                # p_I at min box point = D
                for (k, l, p) in pv:
                    x, y = Q[A[k-1]], Q[B[l-1]]
                    Imin = [Q[a] for a in A[:iD]] + [Q[b] for b in B[:jD]]
                    pD = p_within(L, n, Imin, x, y)
                    pD_tot += 1
                    if third <= pD <= twothird: pD_in += 1
                    # 该对是否平衡
                    bal = third <= p <= twothird
                    hist[(str(pD), bal)] += 1
                    if shown < 4 and len(dp) >= 2:
                        print(f"\nn={n} z={z} rel={rel} A={A} B={B}")
                        print(f"  iD={iD} jD={jD} iU={iU} jU={jU} box={[(k, str(v)) for k, v in sorted(an['mu'].items())]}")
                        for (kk, ll, pp) in pv:
                            xx, yy = Q[A[kk-1]], Q[B[ll-1]]
                            pdv = p_within(L, n, Imin, xx, yy)
                            print(f"    pair(a{kk},b{ll}) p_P={pp}  p_I(D)={pdv}  in={third<=pp<=twothird}")
                        shown += 1
    print("\nstats:", stats['cases'], "cases; 平衡:", stats['ok'], "不平衡:", stats['bad'])
    print("每案例 DxD 不可比对对数分布:", dict(stats['npairs']))
    print(f"p_D(x,y) 落在 [1/3,2/3] 的比例: {pD_in}/{pD_tot} = {Fraction(pD_in, max(pD_tot,1))}")
    print("(p_D, 该对是否平衡) 直方图 top:")
    for kk, v in hist.most_common(12): print("   ", kk, v)


if __name__ == '__main__':
    main()
