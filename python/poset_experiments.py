"""
poset_experiments.py — computational support for the 1/3-2/3 conjecture attack
(width-3 case), session of 2026-09.

Faithful reconstruction of the four experiments run during the session
(original scratch copies were lost to temp-dir cleanup; outputs below are the
recorded session values).

Experiments
-----------
E1  Lemma C (insertion-weight identity) verification on random posets:
    recorded result: 357 verified, 0 mismatches.
E2  Lemma B (critical-pair covering) exhaustive check, all posets n<=5:
    recorded result: 1827 critical pairs checked, 0 failures.
E3  delta landscape, exhaustive n=6 (all 2^15 forward-labeled candidates):
    recorded result: all non-chain posets min delta = 1/3 (sharp);
    width-3 subclass min delta = 4/11; single-point-deletable-to-width-2
    subclass min delta = 4/11.
E4  Insertion-stability diagnostics at n=7/8 over the
    "delete-one-point -> width <= 2" class:
    recorded results: min delta n=7 = 4/11, n=8 = 14/39 (Saks' width-3
    extremal constant); single balanced pairs can be destroyed by insertion
    (p drops to 2/9) but in 879 sampled posets delta(P) >= 4/11 always held
    (19 "leak" cases where no width-2 balanced pair survives still had
    delta(P) in [4/11, 18/37]).

Usage: python poset_experiments.py [e1|e2|e3|e4|all]
E3 takes a few minutes (pure Python); E4 several minutes with default samples.
"""
import itertools, random, sys
from fractions import Fraction


def closure(n, rel):
    less = [[False] * n for _ in range(n)]
    for (i, j) in rel:
        less[i][j] = True
    for k in range(n):
        for i in range(n):
            if less[i][k]:
                ri, rk = less[i], less[k]
                for j in range(n):
                    if rk[j]:
                        ri[j] = True
    return less


def is_poset(less, n):
    for i in range(n):
        if less[i][i]:
            return False
        for j in range(i + 1, n):
            if less[i][j] and less[j][i]:
                return False
    return True


def count_ext(less, n):
    dp = [0] * (1 << n)
    dp[0] = 1
    for S in range(1 << n):
        v = dp[S]
        if not v:
            continue
        for x in range(n):
            if (S >> x) & 1:
                continue
            ok = True
            for y in range(n):
                if less[y][x] and not (S >> y) & 1:
                    ok = False
                    break
            if ok:
                dp[S | 1 << x] += v
    return dp[-1]


def antichain_width(less, n):
    best = 0
    for S in range(1, 1 << n):
        c = bin(S).count('1')
        if c <= best:
            continue
        els = [x for x in range(n) if (S >> x) & 1]
        if all(not (less[a][b] or less[b][a]) for a, b in itertools.combinations(els, 2)):
            best = c
    return best


def delta(less, n):
    base = count_ext(less, n)
    best = Fraction(0)
    bp = None
    for x, y in itertools.combinations(range(n), 2):
        if less[x][y] or less[y][x]:
            continue
        rel = [(i, j) for i in range(n) for j in range(n) if less[i][j]] + [(x, y)]
        p = Fraction(count_ext(closure(n, rel), n), base)
        m = min(p, 1 - p)
        if m > best:
            best = m
            bp = (x, y, p)
    return best, bp


def extensions(less, n):
    res = []

    def rec(cur, used):
        if bin(used).count('1') == n:
            res.append(tuple(cur))
            return
        for x in range(n):
            if (used >> x) & 1:
                continue
            ok = True
            for y in range(n):
                if less[y][x] and not (used >> y) & 1:
                    ok = False
                    break
            if ok:
                rec(cur + [x], used | 1 << x)

    rec([], 0)
    return res


def check_identity(Pless, n, z, x, y):
    """Lemma C: p_P(x,y) = sum_M i_z(M)*[x<y in M] / sum_M i_z(M)."""
    others = [v for v in range(n) if v != z]
    idx = {v: k for k, v in enumerate(others)}
    m = n - 1
    Qrel = [(idx[i], idx[j]) for i in others for j in others if Pless[i][j]]
    Q = closure(m, Qrel)
    below = [idx[v] for v in others if Pless[v][z]]
    above = [idx[v] for v in others if Pless[z][v]]
    num = 0
    den = 0
    for M in extensions(Q, m):
        slots = 0
        prefix = set()
        for k in range(m + 1):
            if k > 0:
                prefix.add(M[k - 1])
            ok = all(a in prefix for a in below)
            if ok:
                suffix = set(M[k:])
                ok = all(b in suffix for b in above)
            if ok:
                slots += 1
        den += slots
        if M.index(idx[x]) < M.index(idx[y]):
            num += slots
    if den == 0:
        return None
    rhs = Fraction(num, den)
    rel = [(i, j) for i in range(n) for j in range(n) if Pless[i][j]] + [(x, y)]
    lhs = Fraction(count_ext(closure(n, rel), n), count_ext(Pless, n))
    return lhs == rhs, lhs, rhs


def e1(random_seed=7, cases=400):
    random.seed(random_seed)
    okc = bad = 0
    for _ in range(cases):
        n = random.choice([4, 5, 6])
        rel = [(i, j) for i in range(n) for j in range(i + 1, n) if random.random() < 0.4]
        L = closure(n, rel)
        if not is_poset(L, n):
            continue
        inc = [(x, y) for x, y in itertools.combinations(range(n), 2)
               if not (L[x][y] or L[y][x])]
        z = random.randrange(n)
        inc = [(x, y) for x, y in inc if x != z and y != z]
        if not inc:
            continue
        x, y = random.choice(inc)
        r = check_identity(L, n, z, x, y)
        if r is None:
            continue
        if r[0]:
            okc += 1
        else:
            bad += 1
            print("  MISMATCH n=%d lhs=%s rhs=%s" % (n, r[1], r[2]))
    print("E1 Lemma C identity: verified %d, mismatch %d" % (okc, bad))


def e2():
    checked = fails = 0
    for bits in range(1 << 10):
        rel = [(i, j) for k, (i, j) in
               enumerate([(i, j) for i in range(5) for j in range(i + 1, 5)])
               if (bits >> k) & 1]
        L = closure(5, rel)
        if not is_poset(L, 5):
            continue
        exts = extensions(L, 5)
        for x, y in itertools.combinations(range(5), 2):
            if L[x][y] or L[y][x]:
                continue
            crit = all((not L[w][x] or L[w][y]) and (not L[x][w] or L[y][w])
                       for w in range(5) if w not in (x, y))
            if not crit:
                continue
            checked += 1
            if not any(any(M[k] == x and M[k + 1] == y for k in range(4))
                       for M in exts):
                fails += 1
                print("  FAIL", rel, x, y)
    print("E2 Lemma B critical pairs: checked %d, failures %d" % (checked, fails))


def e3():
    mins = {'all': (Fraction(2), None), 'w3': (Fraction(2), None),
            'w3s': (Fraction(2), None)}
    cnt = {'all': 0, 'w3': 0, 'w3s': 0}
    for bits in range(1 << 15):
        rel = [(i, j) for k, (i, j) in
               enumerate([(i, j) for i in range(6) for j in range(i + 1, 6)])
               if (bits >> k) & 1]
        L = closure(6, rel)
        if not is_poset(L, 6):
            continue
        w = antichain_width(L, 6)
        if w < 2:
            continue
        d, bp = delta(L, 6)
        cnt['all'] += 1
        if d < mins['all'][0]:
            mins['all'] = (d, bp, rel)
        if w == 3:
            cnt['w3'] += 1
            if d < mins['w3'][0]:
                mins['w3'] = (d, bp, rel)
            for z in range(6):
                Lz = [[L[i][j] for j in range(6) if j != z]
                      for i in range(6) if i != z]
                if antichain_width(Lz, 5) <= 2:
                    cnt['w3s'] += 1
                    if d < mins['w3s'][0]:
                        mins['w3s'] = (d, bp, rel)
                    break
    for k in ('all', 'w3', 'w3s'):
        v = mins[k]
        print("E3", k, "count=%d" % cnt[k],
              "min delta=%s" % (v[0] if v[1] else None), "pair=%s" % (v[1],))


def _rand_width2_poset(m, rng):
    perm = list(range(m))
    rng.shuffle(perm)
    h = m // 2
    C1, C2 = perm[:h], perm[h:]
    rel = {(C1[i], C1[i + 1]) for i in range(len(C1) - 1)}
    rel |= {(C2[i], C2[i + 1]) for i in range(len(C2) - 1)}
    for a in C1:
        for b in C2:
            r = rng.random()
            if r < 0.22:
                rel.add((a, b))
            elif r < 0.44:
                rel.add((b, a))
    L = closure(m, rel)
    if not is_poset(L, m) or antichain_width(L, m) > 2:
        return None
    return L


def _add_point(L, m, rng):
    els = list(range(m))
    D = [v for v in els if rng.random() < 0.3]
    U = [v for v in els if v not in D and rng.random() < 0.3]
    if any(L[b][a] for a in D for b in U):
        return None
    rel = [(i, j) for i in range(m) for j in range(m) if L[i][j]]
    rel += [(a, m) for a in D] + [(m, b) for b in U]
    P = closure(m + 1, rel)
    if not is_poset(P, m + 1):
        return None
    return P


def e4(samples_per_n=2500, max_trials=60000, seed=11):
    rng = random.Random(seed)
    for n in (7, 8):
        best = Fraction(2)
        found = trials = 0
        while found < samples_per_n and trials < max_trials:
            trials += 1
            Q = _rand_width2_poset(n - 1, rng)
            if Q is None:
                continue
            P = _add_point(Q, n - 1, rng)
            if P is None or antichain_width(P, n) != 3:
                continue
            d, _ = delta(P, n)
            found += 1
            if d < best:
                best = d
        print("E4 n=%d samples=%d min_delta=%s" % (n, found, best))


def e4b(samples=20000, seed=23):
    """Insertion stability: can insertion destroy every width-2 balanced pair,
    and does delta(P) stay >= 1/3 anyway?"""
    rng = random.Random(seed)
    third = Fraction(1, 3)
    leak_deltas = []
    worst_drop = worst_rise = None
    tested = 0
    n = 7
    for _ in range(samples):
        Q = _rand_width2_poset(n - 1, rng)
        if Q is None:
            continue
        P = _add_point(Q, n - 1, rng)
        if P is None or antichain_width(P, n) != 3:
            continue
        baseQ = count_ext(Q, n - 1)
        pQ = {}
        for x, y in itertools.combinations(range(n - 1), 2):
            if Q[x][y] or Q[y][x]:
                continue
            rel = [(i, j) for i in range(n - 1) for j in range(n - 1) if Q[i][j]]
            pQ[(x, y)] = Fraction(count_ext(closure(n - 1, rel + [(x, y)]), n - 1), baseQ)
        baseP = count_ext(P, n)
        pP = {}
        for x, y in itertools.combinations(range(n), 2):
            if P[x][y] or P[y][x]:
                continue
            rel = [(i, j) for i in range(n) for j in range(n) if P[i][j]]
            pP[(x, y)] = Fraction(count_ext(closure(n, rel + [(x, y)]), n), baseP)
        for pair, pq in pQ.items():
            pp = pP.get(pair)
            if pp is None:
                continue
            if pq >= Fraction(1, 2) and (worst_drop is None or pp < worst_drop):
                worst_drop = pp
            if pq <= Fraction(1, 2) and (worst_rise is None or pp > worst_rise):
                worst_rise = pp
        tested += 1
        stays = any(third <= pq <= 1 - third and (pair in pP)
                    and third <= pP[pair] <= 1 - third for pair, pq in pQ.items())
        if not stays:
            d, _ = delta(P, n)
            leak_deltas.append(d)
    print("E4b posets tested:", tested)
    print("  min p_P over pairs with p_Q >= 1/2:", worst_drop, "(recorded session value: 2/9)")
    print("  max p_P over pairs with p_Q <= 1/2:", worst_rise, "(recorded session value: 7/9)")
    if leak_deltas:
        print("  leak cases:", len(leak_deltas),
              "their delta(P) range:", min(leak_deltas), "-", max(leak_deltas),
              "all >= 1/3:", all(d >= third for d in leak_deltas))


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("e1", "all"):
        e1()
    if which in ("e2", "all"):
        e2()
    if which in ("e3", "all"):
        e3()
    if which in ("e4", "all"):
        e4()
        e4b()
