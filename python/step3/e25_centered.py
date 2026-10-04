"""E25: centered-z central symmetry.
Claim: if z's D/U sets are mirror-symmetric (iD = s - iU, jD = t - jU), then
gap(Phi(M)) = gap(M) under the reversal+relabel bijection, hence
p_P(k,l) + p_P(s+1-k, t+1-l) = 1 EXACTLY (central symmetry survives in P).
Verify on constructed centered-z corridors; also check odd-odd center cell = 1/2.
"""
import random
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e11_corner.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e22b_s3.py").resolve(), encoding="utf-8").read().split("random.seed(314)")[0].replace("def pmatrix", "def pmatrixOLD"))
exec(open((Path(__file__).parent / "e22b_s3.py").resolve(), encoding="utf-8").read().split("def pmatrix(less")[1].split("random.seed")[0].join(["def pmatrix(less", ""]))

third = Fraction(1, 3)
random.seed(555)

ok_sym = bad_sym = ok_center = na_center = 0
tested = 0
for trial in range(400):
    # random width-2 corridor with cross relations + centered z
    s = random.choice([2, 3, 4])
    t = random.choice([2, 3, 4])
    n = s + t + 1
    z = n - 1
    # chains A = 0..s-1, B = s..s+t-1, random cross relations
    rel = [(i, j) for i in range(s) for j in range(i + 1, s)]
    rel += [(s + i, s + j) for i in range(t) for j in range(i + 1, t)]
    for i in range(s):
        for j in range(t):
            if random.random() < 0.3:
                rel.append((i, s + j))
            elif random.random() < 0.3:
                rel.append((s + j, i))
    # centered z: choose iD in [0..s], iU = s - iD, jD, jU = t - jD; window nonempty
    iD = random.randrange(0, s // 2 + 1)
    jD = random.randrange(0, t // 2 + 1)
    iU = s - iD
    jU = t - jD
    if iD >= iU and jD >= jU: continue  # no window at all
    # z below U-elements... wait: D = elements BELOW z: A[0..iD-1], B[0..jD-1];
    # U = elements ABOVE z: A[iU..s-1], B[jU..t-1]
    for i in range(iD): rel.append((i, z))
    for i in range(iU, s): rel.append((z, i))
    for j in range(jD): rel.append((s + j, z))
    for j in range(jU, t): rel.append((z, s + j))
    L = closure(n, rel)
    if not is_poset(L, n): continue
    if antichain_width(L, n) < 2: continue
    Qz = [v for v in range(n) if v != z]
    Ql, _ = induced(L, n, Qz)
    if antichain_width(Ql, n - 1) > 2: continue
    parts = chain_partitions(Ql, n - 1)
    if not parts: continue
    tested += 1
    for part in parts:
        A, B = part
        s_, t_ = len(A), len(B)
        if s_ != s or t_ != t: continue  # need the full partition matching our construction
        pm = pmatrix(L, n, z, part, None) if False else None
        # build matrix manually (anatomy not needed)
        Q = [v for v in range(n) if v != z]
        mat = {}
        for a1 in range(s_):
            for b1 in range(t_):
                gx, gy = Q[A[a1]], Q[B[b1]]
                if L[gx][gy]: mat[(a1 + 1, b1 + 1)] = Fraction(1)
                elif L[gy][gx]: mat[(a1 + 1, b1 + 1)] = Fraction(0)
                else: mat[(a1 + 1, b1 + 1)] = pval(L, n, gx, gy)
        sym_ok = True
        for k in range(1, s_ + 1):
            for l in range(1, t_ + 1):
                if mat[(k, l)] + mat[(s_ + 1 - k, t_ + 1 - l)] != 1:
                    sym_ok = False
        if sym_ok: ok_sym += 1
        else: bad_sym += 1
        if s_ % 2 == 1 and t_ % 2 == 1:
            c = mat[((s_ + 1) // 2, (t_ + 1) // 2)]
            if c == Fraction(1, 2): ok_center += 1
            else: na_center += 1  # comparable center cell
print("centered-z cases tested:", tested)
print("central symmetry in P: ok=%d bad=%d | odd-odd center=1/2: ok=%d comparable=%d"
      % (ok_sym, bad_sym, ok_center, na_center))
