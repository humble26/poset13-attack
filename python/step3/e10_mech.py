"""E10: dissect the dominant fail configurations.

For each fail case and partition, dump:
  - window M-profiles and cut config
  - ALL cross-chain Q-pairs with p, categories (D/W lo/hi, U)
  - for balanced pairs: decomposition over box cuts (which cuts deliver balance)
  - p_P(x,y) vs p_Q(x,y) (reweighting shift)
Focus configs: (alllo,alllo) and (jump,jump).
"""
import itertools, sys
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main()")[0])

third = Fraction(1, 3)
twothird = Fraction(2, 3)


def full_pair_table(less, n, z, part, an):
    """all cross-chain pairs with p_P, p_Q, cats; and cut-decomposition for bal pairs"""
    A, B = part
    s, t = an['s'], an['t']
    Q = [v for v in range(n) if v != z]
    Ql, _ = induced(less, n, Q)
    eQ = count_ext(Ql, n - 1)
    rows = []
    for a1 in range(s):
        for b1 in range(t):
            gx, gy = Q[A[a1]], Q[B[b1]]
            k, l = a1 + 1, b1 + 1
            if less[gx][gy]:
                rel = 'xy'
            elif less[gy][gx]:
                rel = 'yx'
            else:
                rel = 'inc'
            pP = pval(less, n, gx, gy) if rel == 'inc' else None
            lx_, ly_ = Q.index(gx), Q.index(gy)
            pQ = Fraction(sum(1 for ext in extensions_of(Ql, n-1) if ext.index(lx_) < ext.index(ly_)), eQ) if rel == 'inc' else None
            cxA = 'D' if k <= an['iD'] else ('U' if k > an['iU'] else 'W')
            cyB = 'D' if l <= an['jD'] else ('U' if l > an['jU'] else 'W')
            fx = 'lo' if (cxA == 'D' or (cxA == 'W' and an['MAv'][k] < third)) else 'hi'
            fy = 'lo' if (cyB == 'D' or (cyB == 'W' and an['MBv'][l] < third)) else 'hi'
            rows.append(dict(k=k, l=l, gx=gx, gy=gy, rel=rel, pP=pP, pQ=pQ,
                             cxA=cxA, cyB=cyB, fx=fx, fy=fy))
    return rows


_ext_cache = {}
def extensions_of(less, m):
    key = id(less)
    if key in _ext_cache: return _ext_cache[key]
    res = []
    def rec(placed, seq):
        if len(seq) == m:
            res.append(list(seq)); return
        for x in range(m):
            if (placed >> x) & 1: continue
            if all((placed >> y) & 1 for y in range(m) if less[y][x]):
                rec(placed | 1 << x, seq + [x])
    rec(0, [])
    _ext_cache[key] = res
    return res


def cut_decomp(less, n, z, part, an, gx, gy):
    """p_P(gx,gy) split over box cuts by region class of (i,j) relative to pair."""
    A, B = part
    Q = [v for v in range(n) if v != z]
    xl, yl = Q.index(gx), Q.index(gy)
    mu = an['mu']
    out = {}
    for (i, j), w in sorted(mu.items()):
        Iset = set(A[:i]) | set(B[:j])  # Ql-indices
        x_in = xl in Iset; y_in = yl in Iset
        tag = ('in' if x_in else 'out') + ('in' if y_in else 'out')
        if tag == 'inout':
            contrib = Fraction(1)
        elif tag == 'outin':
            contrib = Fraction(0)
        elif tag == 'inin':
            xf, yf, tot = inter_count(less, [Q[v] for v in A[:i]], [Q[v] for v in B[:j]], gx, gy)
            contrib = Fraction(xf, tot)
        else:
            xf, yf, tot = inter_count(less, [Q[v] for v in A[i:]], [Q[v] for v in B[j:]], gx, gy)
            contrib = Fraction(xf, tot)
        out[(i, j)] = (w, tag, contrib, w * contrib)
    return out


def ideal_elems(A, B, i, j):
    return set(A[:i]) | set(B[:j])


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    want_cfg = sys.argv[2] if len(sys.argv) > 2 else 'alllo,alllo'
    shown = 0
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
                for (A, B) in chain_partitions(Ql, n - 1):
                    an = anatomy(L, n, z, (A, B))
                    if an is None: continue
                    sideA = 'alllo' if an['iHi'] is None else ('allhi' if an['iLo'] is None else 'jump')
                    sideB = 'alllo' if an['jHi'] is None else ('allhi' if an['jLo'] is None else 'jump')
                    if f"{sideA},{sideB}" != want_cfg: continue
                    if shown >= int(sys.argv[3] if len(sys.argv) > 3 else 4): continue
                    shown += 1
                    print("=" * 70)
                    print("n=%d z=%d rel=%s A=%s B=%s cfg=%s" % (n, z, rel, A, B, (sideA, sideB)))
                    Dg = sorted(v for v in Qz if L[v][z]); Ug = sorted(v for v in Qz if L[z][v])
                    print("  D=%s U=%s box_mu=%s" % (Dg, Ug, {kk: str(v) for kk, v in sorted(an['mu'].items())}))
                    print("  MAv=%s" % {kk: str(v) for kk, v in sorted(an['MAv'].items())})
                    print("  MBv=%s" % {kk: str(v) for kk, v in sorted(an['MBv'].items())})
                    rows = full_pair_table(L, n, z, (A, B), an)
                    for rr in rows:
                        if rr['rel'] != 'inc':
                            print("   pair (a%d,b%d)=(%d,%d) %s %s%s" % (rr['k'], rr['l'], rr['gx'], rr['gy'], rr['rel'], rr['cxA'], rr['cyB']))
                        else:
                            print("   pair (a%d,b%d)=(%d,%d) inc %s%s/%s%s pP=%s pQ=%s" % (
                                rr['k'], rr['l'], rr['gx'], rr['gy'], rr['cxA'], rr['fx'], rr['cyB'], rr['fy'], rr['pP'], rr['pQ']))
                    # decomposition for first balanced pair
                    for rr in rows:
                        if rr['rel'] == 'inc' and third <= rr['pP'] <= twothird:
                            dd = cut_decomp(L, n, z, (A, B), an, rr['gx'], rr['gy'])
                            print("   DECOMP pP(%d,%d) over cuts:" % (rr['gx'], rr['gy']))
                            for (ij, (w, tag, c, wc)) in dd.items():
                                print("      cut %s w=%s %s contrib=%s -> %s" % (ij, w, tag, c, wc))
                            break
                    if shown >= int(sys.argv[3] if len(sys.argv) > 3 else 4):
                        break
                if shown >= int(sys.argv[3] if len(sys.argv) > 3 else 4): break
            if shown >= int(sys.argv[3] if len(sys.argv) > 3 else 4): break
        if shown >= int(sys.argv[3] if len(sys.argv) > 3 else 4): break


if __name__ == '__main__':
    main()
