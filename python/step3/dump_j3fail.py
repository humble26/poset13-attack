"""dump J3 失效案例（不存在"全 I 平衡"的 DxD 对），展示平均机制。"""
import sys
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e11_corner.py").resolve(), encoding="utf-8").read().split("def main():")[0])

third = Fraction(1, 3); twothird = Fraction(2, 3)


def p_within(less, n, Ilist, x, y):
    Il, idx = induced(less, n, Ilist)
    m = len(Ilist); eI = count_ext(Il, m)
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
                if len(an['mu']) < 2 or iD < 1 or jD < 1: continue
                A, B = part; Q = [v for v in range(n) if v != z]
                dp = [(k, l) for k in range(1, iD+1) for l in range(1, jD+1)
                      if not (L[Q[A[k-1]]][Q[B[l-1]]] or L[Q[B[l-1]]][Q[A[k-1]]])]
                if not dp: continue
                # 计算每对的 p_I 向量
                info = []
                for (k, l) in dp:
                    x, y = Q[A[k-1]], Q[B[l-1]]
                    vals = {}
                    for (i, j), w in sorted(an['mu'].items()):
                        Ilist = [Q[a] for a in A[:i]] + [Q[b] for b in B[:j]]
                        vals[(i, j)] = p_within(L, n, Ilist, x, y)
                    pp = pval(L, n, x, y)
                    allin = all(third <= v <= twothird for v in vals.values())
                    info.append((k, l, pp, allin, vals))
                if any(t[3] for t in info): continue  # 有全 I 平衡的，跳过
                shown += 1
                if shown > 5: continue
                print("=" * 68)
                print(f"n={n} z={z} rel={rel} A={A} B={B} iD={iD} jD={jD} iU={an['iU']} jU={an['jU']}")
                print(f"  box mu: {[(kk, str(v)) for kk, v in sorted(an['mu'].items())]}")
                for (k, l, pp, allin, vals) in info:
                    vv = " ".join(f"{kk}:{vals[kk]}" for kk in sorted(vals))
                    print(f"  (a{k},b{l}) p_P={pp}  全I平衡={allin}   p_I: {vv}")
    print(f"\n共 {shown} 个 J3 失效案例")


if __name__ == '__main__':
    main()
