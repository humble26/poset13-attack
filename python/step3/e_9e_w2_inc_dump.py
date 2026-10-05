"""dump 一个 |W|=2 不可比对（w1∥w2）案例的完整结构，理解 m={5/18,13/18} 的来源。"""
import sys, itertools, collections
from fractions import Fraction
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from my_poset import my_close, my_is_poset, my_induced, my_count_ext


def width_le2(Q, mm):
    for a, b, c in itertools.combinations(range(mm), 3):
        if not (Q[a][b] or Q[b][a] or Q[a][c] or Q[c][a] or Q[b][c] or Q[c][b]): return False
    return True


def chain_sort(Q, els):
    return sorted(els, key=lambda v: sum(1 for u in els if Q[u][v]))


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    mm = n - 1
    PAIRS = list(itertools.combinations(range(mm), 2))
    full = (1 << mm) - 1
    third, twothird = Fraction(1, 3), Fraction(2, 3)
    shown = 0
    for mask in range(1 << len(PAIRS)):
        rel = [PAIRS[k] for k in range(len(PAIRS)) if (mask >> k) & 1]
        Q = my_close(mm, rel)
        if not my_is_poset(Q, mm) or not width_le2(Q, mm): continue
        ids = [S for S in range(1 << mm)
               if all(not ((S >> i) & 1) or all((S >> j) & 1 for j in range(mm) if Q[j][i]) for i in range(mm))]
        e_all = [0] * (1 << mm); e_all[0] = 1
        for S in range(1, 1 << mm):
            tot = 0
            for x in range(mm):
                if not (S >> x) & 1: continue
                if all(not Q[y][x] or not (S >> y) & 1 for y in range(mm)):
                    tot += e_all[S ^ (1 << x)]
            e_all[S] = tot
        wgt = {I: e_all[I] * e_all[full ^ I] for I in ids}
        USM = [S for S in range(1 << mm)
               if all(not ((S >> i) & 1) or all((S >> j) & 1 for j in range(mm) if Q[i][j]) for i in range(mm))]
        for D in ids:
            for U in USM:
                if D & U: continue
                Wl = [v for v in range(mm) if not ((D >> v) & 1) and not ((U >> v) & 1)]
                if len(Wl) != 2: continue
                w1, w2 = Wl
                if Q[w1][w2] or Q[w2][w1]: continue
                Dl = [v for v in range(mm) if (D >> v) & 1]
                Ul = [v for v in range(mm) if (U >> v) & 1]
                if any(not (Q[a][b] or Q[b][a]) for i, a in enumerate(Dl) for b in Dl[i+1:]): continue
                if any(not (Q[a][b] or Q[b][a]) for i, a in enumerate(Ul) for b in Ul[i+1:]): continue
                box = [I for I in ids if (I & D) == D and (I & U) == 0]
                eP = sum(wgt[I] for I in box)
                if eP == 0 or len(box) < 2: continue
                fail = True
                m = {}
                for v in Wl:
                    ev = sum(wgt[I] for I in box if not (I >> v) & 1)
                    m[v] = Fraction(ev, eP)
                    if third <= m[v] <= twothird: fail = False
                if not fail: continue
                if shown >= 3: continue
                shown += 1
                print("=" * 70)
                print(f"mask={mask} D={Dl} U={Ul} W={Wl}  eP={eP}")
                print(f"  m(w1)={m[w1]}  m(w2)={m[w2]}")
                # 四个盒点及其 e(I), e(Q∖I)
                print("  盒点 (w1w2 ∈ I):")
                for I in sorted(box):
                    sub = tuple(w for w in Wl if (I >> w) & 1)
                    eI = e_all[I]; eC = e_all[full ^ I]
                    print(f"    {sub or '∅'}: e(I)={eI} e(Q∖I)={eC} μ={Fraction(eI*eC, eP)}")
                # D、U 链长与 w 的关系
                for w in Wl:
                    belowD = [d for d in Dl if Q[d][w]]
                    aboveU = [u for u in Ul if Q[w][u]]
                    print(f"  w={w}: D中低于w {belowD}；U中高于w {aboveU}")
                # Q 的整体关系（D×U, w1w2 之外）
                print(f"  Q 关系: {sorted((a,b) for (a,b) in PAIRS if Q[a][b])}")


if __name__ == '__main__':
    main()
