"""DIAG: 打印 9.D 案例的完整结构（box / mu / 每个 DxD 对的 p_I 分解）。
目标：看清"某个 DxD 对平衡"是怎么来的。"""
import itertools, sys
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e11_corner.py").resolve(), encoding="utf-8").read().split("def main():")[0])

third = Fraction(1, 3); twothird = Fraction(2, 3)


def p_within(less, n, Ilist, x, y):
    """p_I(x,y) = e(I+x<y)/e(I) 在子偏序 I 内。"""
    Il, idx = induced(less, n, Ilist)
    m = len(Ilist)
    eI = count_ext(Il, m)
    # add x<y
    rel2 = [(idx[a], idx[b]) for a in Ilist for b in Ilist if Il[idx[a]][idx[b]] if a != b]
    rel2 = [(i, j) for (i, j) in rel2 if not Il[i][j] or True]
    L2 = closure(m, [(i, j) for i in range(m) for j in range(m) if Il[i][j]] + [(idx[x], idx[y])])
    e2 = count_ext(L2, m)
    return Fraction(e2, eI)


def main():
    found = 0
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    for n in range(4, nmax + 1):
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
                    an = cc['an']
                    iD, jD, iU, jU = an['iD'], an['jD'], an['iU'], an['jU']
                    if len(an['mu']) < 2: continue
                    A, B = part
                    Q = [v for v in range(n) if v != z]
                    s_, t_ = an['s'], an['t']
                    if iD < 1 or jD < 1: continue
                    dpairs = []
                    for k in range(1, iD + 1):
                        for l in range(1, jD + 1):
                            x, y = Q[A[k-1]], Q[B[l-1]]
                            if not (L[x][y] or L[y][x]):
                                dpairs.append((k, l, x, y))
                    if not dpairs: continue
                    pv = [(k, l, pval(L, n, x, y)) for (k, l, x, y) in dpairs]
                    if any(third <= p <= twothird for (_, _, p) in pv):
                        found += 1
                        if found > 6: continue
                        print("=" * 72)
                        print(f"n={n} z={z} rel={rel}  A={A} B={B}  abort={r}")
                        print(f"  s={s_} t={t_} iD={iD} jD={jD} iU={iU} jU={jU}  eP={an['eP']}")
                        print(f"  box mu = {[(k, str(v)) for k, v in sorted(an['mu'].items())]}")
                        print(f"  MAv={[(k, str(v)) for k, v in sorted(an['MAv'].items())]}")
                        print(f"  MBv={[(k, str(v)) for k, v in sorted(an['MBv'].items())]}")
                        print(f"  DxD pairs: {[(k, l, str(p)) for (k, l, p) in pv]}")
                        # 对每个 DxD 对，列出 p_I 分解
                        for (k, l, p) in pv:
                            x, y = Q[A[k-1]], Q[B[l-1]]
                            terms = []
                            for (i, j), w in sorted(an['mu'].items()):
                                Ilist = [Q[a] for a in A[:i]] + [Q[b] for b in B[:j]]
                                pi = p_within(L, n, Ilist, x, y)
                                terms.append((i, j, str(w), str(pi)))
                            print(f"    ({k},{l}) p_P={p}  <-- sum mu*p_I:")
                            for (i, j, w, pi) in terms:
                                print(f"        I({i},{j}) mu={w} p_I={pi}")
    print(f"\n统一发现 {found} 个 9.D 案例（含平衡 DxD 对）")


if __name__ == '__main__':
    main()
