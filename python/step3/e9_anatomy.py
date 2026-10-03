"""E9: anatomy of zpair-fail cases (conservation mechanism).

For P = Q + z with width(Q) <= 2, in a FAIL case (no balanced z-pair):
  - verify structural facts: window continuity, M_A monotone, heavy layer,
    squeeze M_A(i0) in (2/3 - colmass[i0], 1/3) at jump columns;
  - find all balanced Q-pairs and classify them (D/W/U x D/W/U, low x high?);
  - for union-bound rescue pairs decompose p(y,x) = mu_L + W_in + W_out
    to locate the source of the missing UPPER bound p <= 2/3.
"""
import itertools, sys
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])

third = Fraction(1, 3)
twothird = Fraction(2, 3)


def chain_partitions(Ql, m):
    """all partitions of Ql (width<=2) into two chains, up to swap."""
    els = list(range(m))
    res = []
    for r in range(1, m):
        for S in itertools.combinations(els, r):
            Sset = set(S); Tset = set(els) - set(S)
            def is_chain(C):
                C = list(C)
                for a, b in itertools.combinations(C, 2):
                    if not (Ql[a][b] or Ql[b][a]): return False
                return True
            if is_chain(Sset) and is_chain(Tset):
                A = sorted(Sset, key=lambda v: sum(1 for u in Sset if Ql[u][v]))
                B = sorted(Tset, key=lambda v: sum(1 for u in Tset if Ql[u][v]))
                if tuple(A) <= tuple(B):
                    res.append((A, B))
    return res


def inter_count(less, elemsA, elemsB, x, y):
    """extensions of chain poset elemsA+elemsB with x before y (and y before x).
    DP over (ia, ib, state); prefixes must be down-closed in the sub-poset."""
    s, t = len(elemsA), len(elemsB)
    sub = list(elemsA) + list(elemsB)
    pos = {v: i for i, v in enumerate(sub)}
    pred = []
    for u in sub:
        msk = 0
        for v in sub:
            if v != u and less[v][u]: msk |= 1 << pos[v]
        pred.append(msk)
    xi, yi = pos[x], pos[y]
    ax = xi < s  # x on A chain?
    # state: 0 none, 1 only-x, 2 only-y, 3 both x-first, 4 both y-first
    dp = {(0, 0, 0): 1}
    for ia in range(s + 1):
        for ib in range(t + 1):
            for st in (0, 1, 2, 3, 4):
                v = dp.get((ia, ib, st))
                if not v: continue
                placed = 0
                for k in range(ia): placed |= 1 << k
                for k in range(ib): placed |= 1 << (s + k)
                if ia < s:
                    u = elemsA[ia]; up = pos[u]
                    if pred[up] & ~placed == 0:
                        nst = st
                        if up == xi: nst = {0: 1, 2: 4}.get(st, st)
                        if up == yi: nst = {0: 2, 1: 3}.get(st, st)
                        dp[(ia + 1, ib, nst)] = dp.get((ia + 1, ib, nst), 0) + v
                if ib < t:
                    u = elemsB[ib]; up = pos[u]
                    if pred[up] & ~placed == 0:
                        nst = st
                        if up == xi: nst = {0: 1, 2: 4}.get(st, st)
                        if up == yi: nst = {0: 2, 1: 3}.get(st, st)
                        dp[(ia, ib + 1, nst)] = dp.get((ia, ib + 1, nst), 0) + v
    tot = sum(dp.get((s, t, st), 0) for st in (0, 1, 2, 3, 4))
    xf = dp.get((s, t, 3), 0); yf = dp.get((s, t, 4), 0)
    assert xf + yf == tot, (xf, yf, tot)
    return xf, yf, tot


def anatomy(less, n, z, part):
    """full anatomy of one (P, z, chain-partition). Returns dict or None."""
    Q = [v for v in range(n) if v != z]
    Ql, idx = induced(less, n, Q); m = n - 1
    A, B = part
    s, t = len(A), len(B)
    # D = {v in Q : Q[v] < z}, U = {v in Q : z < Q[v]}  (Ql indices)
    D = set(); U = set()
    for v in range(m):
        g = Q[v]
        if less[g][z]: D.add(v)
        if less[z][g]: U.add(v)
    # chains as Ql-index lists A, B (already chain-sorted)
    def pref_count(chain, S):
        # S (set) must be prefix of chain; return length
        c = 0
        for v in chain:
            if v in S: c += 1
            else: break
        # sanity: prefix property
        if S & set(chain[c:]): return None
        return c
    iD = pref_count(A, D); jD = pref_count(B, D)
    # U is up-set: U cap A is suffix
    iU = s - sum(1 for v in A if v in U); jU = t - sum(1 for v in B if v in U)
    if iD is None or jD is None: return None
    if any(v in U for v in A[:iU]) or any(v in U for v in B[:jU]): return None
    eP = count_ext(less, n)
    # box: ideals of Ql with D<=S, ScapU=empty -> (i,j) grid
    box = {}  # (i,j) -> [eI, eF]
    for S in ideals(Ql, m):
        if (sum(1 for v in A if (S >> v) & 1) != len([v for v in A if (S >> v) & 1])): pass
        Dmask = 0; 
        for v in D: Dmask |= 1 << v
        Umask = 0
        for v in U: Umask |= 1 << v
        if (S & Dmask) != Dmask or (S & Umask): continue
        i = sum(1 for v in A if (S >> v) & 1)
        j = sum(1 for v in B if (S >> v) & 1)
        Ilist = [Q[k] for k in range(m) if (S >> k) & 1]
        Clist = [Q[k] for k in range(m) if not (S >> k) & 1]
        Il, _ = induced(less, n, Ilist); Cl, _ = induced(less, n, Clist)
        box[(i, j)] = [count_ext(Il, len(Ilist)), count_ext(Cl, len(Clist))]
    if not box: return None
    tot = sum(e * f for (e, f) in box.values())
    assert tot == eP, (tot, eP)
    mu = {k: Fraction(v[0] * v[1], eP) for k, v in box.items()}
    colmass = {}; rowmass = {}
    for (i, j), w in mu.items():
        colmass[i] = colmass.get(i, 0) + w
        rowmass[j] = rowmass.get(j, 0) + w
    # M_A(k) = mu(i < k)  (k = 1-based chain index of a_k)
    def MA(k): return sum(w for (i, j), w in mu.items() if i < k)
    def MB(l): return sum(w for (i, j), w in mu.items() if j < l)
    # verify strip formula on window elements
    for k in range(1, s + 1):
        g = Q[A[k - 1]]
        if less[g][z] or less[z][g]: continue
        assert MA(k) == pval(less, n, z, g), ("strip A", k)
    for l in range(1, t + 1):
        g = Q[B[l - 1]]
        if less[g][z] or less[z][g]: continue
        assert MB(l) == pval(less, n, z, g), ("strip B", l)
    # windows (1-based chain indices, contiguous (iD, iU])
    winA = [k for k in range(1, s + 1) if not (less[Q[A[k - 1]]][z] or less[z][Q[A[k - 1]]])]
    winB = [l for l in range(1, t + 1) if not (less[Q[B[l - 1]]][z] or less[z][Q[B[l - 1]]])]
    assert winA == list(range(iD + 1, iU + 1)), (winA, iD, iU)
    assert winB == list(range(jD + 1, jU + 1)), (winB, jD, jU)
    MAv = {0: Fraction(0)}
    for k in range(1, s + 2): MAv[k] = MA(k)
    MBv = {0: Fraction(0)}
    for l in range(1, t + 2): MBv[l] = MB(l)
    # monotonicity check
    for k in range(1, s + 1): assert MAv[k] <= MAv[k + 1]
    for l in range(1, t + 1): assert MBv[l] <= MBv[l + 1]
    iLo = max([k for k in winA if MAv[k] < third], default=None)
    iHi = min([k for k in winA if MAv[k] > twothird], default=None)
    jLo = max([l for l in winB if MBv[l] < third], default=None)
    jHi = min([l for l in winB if MBv[l] > twothird], default=None)
    # heavy layer check: jump column or terminal/first column heavy
    heavy = None
    for c, w in colmass.items():
        if w > third:
            heavy = ("col", c, w); break
    if heavy is None:
        for r, w in rowmass.items():
            if w > third:
                heavy = ("row", r, w); break
    # squeeze: consecutive window columns k-1,k both with a_{k-1},a_k in window:
    # colmass[k-1] > 1/3 forced jump -> MAv[k-1] in (2/3 - colmass, 1/3)
    squeezes = []
    for c in range(iD, iU + 1):
        # columns c and previous window cut: if a_{c+1} in window and MAv jumped
        if (c + 1) in winA and colmass.get(c, 0) > third and MAv[c] < third:
            squeezes.append((c, MAv[c], colmass[c]))
    # category of each Q element
    def catA(k):  # 1-based
        if k <= iD: return 'D'
        if k > iU: return 'U'
        return 'W'
    def catB(l):
        if l <= jD: return 'D'
        if l > jU: return 'U'
        return 'W'
    def lowA(k): return catA(k) == 'D' or (catA(k) == 'W' and MAv[k] < third)
    def highA(k): return catA(k) == 'U' or (catA(k) == 'W' and MAv[k] > twothird)
    def lowB(l): return catB(l) == 'D' or (catB(l) == 'W' and MBv[l] < third)
    def highB(l): return catB(l) == 'U' or (catB(l) == 'W' and MBv[l] > twothird)
    # all Q-pairs: balanced?
    bal_pairs = []
    for a1 in range(s):
        for b1 in range(t):
            gx, gy = Q[A[a1]], Q[B[b1]]
            if less[gx][gy] or less[gy][gx]: continue
            p = pval(less, n, gx, gy)
            if third <= p <= 1 - third:
                k, l = a1 + 1, b1 + 1
                bal_pairs.append((A[a1], B[b1], p, catA(k), catB(l),
                                  lowA(k), highA(k), lowB(l), highB(l), k, l))
    # union-bound rescue pairs R1=(a_iLo, b_jHi), R2=(b_jLo, a_iHi)
    rescues = []
    if iLo and jHi:
        gx, gy = Q[A[iLo - 1]], Q[B[jHi - 1]]
        if not (less[gx][gy] or less[gy][gx]):
            p = pval(less, n, gx, gy)
            rescues.append(('R1', A[iLo - 1], B[jHi - 1], p, iLo, jHi))
    if jLo and iHi:
        gx, gy = Q[B[jLo - 1]], Q[A[iHi - 1]]
        if not (less[gy][gx] or less[gx][gy]):
            p = pval(less, n, gy, gx)
            rescues.append(('R2', A[iHi - 1], B[jLo - 1], p, iHi, jLo))
    return dict(eP=eP, mu=mu, colmass=colmass, rowmass=rowmass, s=s, t=t,
                iD=iD, jD=jD, iU=iU, jU=jU, MAv=MAv, MBv=MBv,
                iLo=iLo, iHi=iHi, jLo=jLo, jHi=jHi, heavy=heavy,
                squeezes=squeezes, bal_pairs=bal_pairs, rescues=rescues,
                winA=winA, winB=winB)


def decompose_rev(less, n, z, part, an, k, l):
    """decompose p(b_l, a_k) = mu_L + W_in + W_out for x=a_k, y=b_l."""
    A, B = part
    s, t = an['s'], an['t']
    Q = [v for v in range(n) if v != z]
    gx, gy = Q[A[k - 1]], Q[B[l - 1]]
    mu = an['mu']; eP = an['eP']
    muL = muR = w_in = w_out = Fraction(0)
    for (i, j), w in mu.items():
        x_in = (k <= i); y_in = (l <= j)
        if x_in and not y_in: muR += w
        elif y_in and not x_in: muL += w
        elif x_in and y_in:
            xf, yf, tot = inter_count(less, A[:i], B[:j], gx, gy)
            assert tot > 0
            w_in += w * Fraction(yf, tot)
        else:
            xf, yf, tot = inter_count(less, A[i:], B[j:], gx, gy)
            assert tot > 0
            w_out += w * Fraction(yf, tot)
    return muL, w_in, w_out, muR


def main():
    stats = dict(cases=0, fail=0, parts=0, conserv_ok=0, conserv_bad=0)
    cat_hist = {}
    cutcfg = {}
    ub_rescue = dict(exists_bal=0, exists_incomp=0, cases_with=0)
    rev_ge_third = 0; rev_ge_third_cases = 0
    decomp_terms = dict(muL=0, w_in=0, w_out=0)
    p_rescue_min = Fraction(2); p_rescue_max = Fraction(0)
    squeeze_ok = 0; heavy_ok = 0; heavy_absent = 0
    examples = []
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
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
                stats['cases'] += 1
                if r['zpair']: continue
                stats['fail'] += 1
                parts = chain_partitions(Ql, n - 1)
                stats['parts'] += len(parts)
                case_ub = False; case_rev = False
                for (A, B) in parts:
                    an = anatomy(L, n, z, (A, B))
                    if an is None: continue
                    # heavy layer
                    if an['heavy'] is not None: heavy_ok += 1
                    else: heavy_absent += 1
                    for sq in an['squeezes']:
                        c, m0, cm = sq
                        if twothird - cm < m0 < third: squeeze_ok += 1
                        else: print("SQUEEZE-FAIL", n, z, sq)
                    # conservation
                    if an['bal_pairs']:
                        stats['conserv_ok'] += 1
                    else:
                        stats['conserv_bad'] += 1
                        print("CONSERV-BAD", n, z, rel, (A, B))
                    for (xa, yb, p, cxA, cyB, lx, hx, ly, hy, k, l) in an['bal_pairs']:
                        def flag(cat, low, high, mval):
                            if cat == 'D': return 'lo'
                            if cat == 'U': return 'hi'
                            return 'lo' if mval < third else 'hi'
                        fx = flag(cxA, lx, hx, an['MAv'][k])
                        fy = flag(cyB, ly, hy, an['MBv'][l])
                        cat_hist[(cxA, cyB, fx, fy)] = cat_hist.get((cxA, cyB, fx, fy), 0) + 1
                    sideA = 'alllo' if an['iHi'] is None else ('allhi' if an['iLo'] is None else 'jump')
                    sideB = 'alllo' if an['jHi'] is None else ('allhi' if an['jLo'] is None else 'jump')
                    cutcfg[(sideA, sideB)] = cutcfg.get((sideA, sideB), 0) + 1
                    if len(examples) < 8 and n >= 5:
                        examples.append((n, z, rel, (A, B), an))
                    # union-bound rescues
                    for (tag, gx, gy, p, k, l) in an['rescues']:
                        ub_rescue['exists_incomp'] += 1
                        if third <= p <= 1 - third:
                            ub_rescue['exists_bal'] += 1
                            case_ub = True
                            if p < p_rescue_min: p_rescue_min = p
                            if p > p_rescue_max: p_rescue_max = p
                            if len(examples) < 12 and n >= 5:
                                examples.append((n, z, rel, tag, (A, B), an, k, l, p))
                            # reverse decomposition
                            muL, w_in, w_out, muR = decompose_rev(L, n, z, (A, B), an, k, l)
                            rev = muL + w_in + w_out
                            if rev >= third:
                                rev_ge_third += 1
                                case_rev = True
                            # dominant term
                            terms = {'muL': muL, 'w_in': w_in, 'w_out': w_out}
                            dom = max(terms, key=terms.get)
                            decomp_terms[dom] += 1
                if case_ub: ub_rescue['cases_with'] += 1
                if case_rev: rev_ge_third_cases += 1
    print("stats:", stats)
    print("cut configs (sideA,sideB):", cutcfg)
    print("balanced-pair hist (catA,catB,flagx,flagy):")
    for kk in sorted(cat_hist): print("   ", kk, cat_hist[kk])
    print("ub_rescue:", ub_rescue)
    print("p(rescue) range:", p_rescue_min, p_rescue_max)
    print("reverse >= 1/3:", rev_ge_third, "cases:", rev_ge_third_cases)
    print("dominant term of p(y,x):", decomp_terms)
    print("heavy present:", heavy_ok, "absent:", heavy_absent, "squeeze_ok:", squeeze_ok)
    for (n, z, rel, part, an) in examples[:8]:
        print("--- example n=%d z=%d rel=%s" % (n, z, rel))
        print("    A=%s B=%s" % (part[0], part[1]))
        print("    box_mu=%s" % {kk: str(v) for kk, v in sorted(an['mu'].items())})
        print("    MAv=%s" % {kk: str(v) for kk, v in sorted(an['MAv'].items())})
        print("    MBv=%s" % {kk: str(v) for kk, v in sorted(an['MBv'].items())})
        print("    heavy=%s" % (an['heavy'],))
        print("    bal=%s" % [(cxA, cyB, str(pp)) for (_, _, pp, cxA, cyB, *_ ) in an['bal_pairs']])


if __name__ == '__main__':
    main()
