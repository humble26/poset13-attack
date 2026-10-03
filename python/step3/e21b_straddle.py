"""E21b: do TRUE over-jumps exist at n=6? (adjacent incomparable band cells
with p1 > 2/3 and p2 < 1/3). If yes, record where the balanced cell actually
is relative to the jumping column (same column? another row?).
"""
import random
from fractions import Fraction
from pathlib import Path
exec(open((Path(__file__).parent / "../step2/e5_ideal.py").resolve(), encoding="utf-8").read().split("random.seed(5)")[0])
exec(open((Path(__file__).parent / "e9_anatomy.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e12_c1.py").resolve(), encoding="utf-8").read().split("def main():")[0])
exec(open((Path(__file__).parent / "e20_walk.py").resolve(), encoding="utf-8").read().split("def partA()")[0].split("if __name__")[0].replace("def main", "def _unused"))

third = Fraction(1, 3); twothird = Fraction(2, 3)
random.seed(99)

n = 6
npairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
over = 0
gentle = 0
tried = 0
examples = []
while tried < 6000 and len(examples) < 6:
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
        pm = pmatrix(L, n, z, part, an)
        A, B = part; s, t = an['s'], an['t']
        # over-jumps: adjacent band cells in a column straddling with p1>2/3, p2<1/3
        for l in range(1, t + 1):
            col = [(k, pm[(k, l)][1]) for k in range(1, s + 1)
                   if (k, l) in pm and pm[(k, l)][0] == 'inc']
            for idx in range(len(col) - 1):
                (k1, p1), (k2, p2) = col[idx], col[idx + 1]
                if p1 > twothird and p2 < third:
                    over += 1
                    # where is the balanced cell? same column? row k1? row k2?
                    bal_in_col = [(kk, str(pp)) for (kk, pp) in col
                                  if third <= pp <= twothird]
                    bal_in_row = [(ll, str(pm[(k1, ll)][1])) for ll in range(1, t + 1)
                                  if (k1, ll) in pm and pm[(k1, ll)][0] == 'inc'
                                  and third <= pm[(k1, ll)][1] <= twothird]
                    if len(examples) < 6:
                        examples.append((rel, z, part, ('col', l, k1, k2, str(p1), str(p2)),
                                         bal_in_col, bal_in_row))
                elif p1 >= twothird and p2 <= twothird and p2 >= third:
                    gentle += 1
print("over-jumps:", over, "gentle-landings:", gentle, "tried:", tried)
for e in examples:
    print(e)
