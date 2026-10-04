"""统计 n=8（Q 有 7 个元素）全穷举规模。"""
import sys, itertools, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from my_poset import my_close, my_is_poset


def width_le2(Q, mm):
    for a, b, c in itertools.combinations(range(mm), 3):
        if not (Q[a][b] or Q[b][a] or Q[a][c] or Q[c][a] or Q[b][c] or Q[c][b]):
            return False
    return True


n = int(sys.argv[1]); mm = n - 1
PAIRS = list(itertools.combinations(range(mm), 2))
t0 = time.time()
nQ = 0; nInst = 0; nQmulti = 0
for mask in range(1 << len(PAIRS)):
    rel = [PAIRS[k] for k in range(len(PAIRS)) if (mask >> k) & 1]
    Q = my_close(mm, rel)
    if not my_is_poset(Q, mm) or not width_le2(Q, mm): continue
    nQ += 1
    DS = [S for S in range(1 << mm)
          if all(not ((S >> i) & 1) or all((S >> j) & 1 for j in range(mm) if Q[j][i]) for i in range(mm))]
    US = [S for S in range(1 << mm)
          if all(not ((S >> i) & 1) or all((S >> j) & 1 for j in range(mm) if Q[i][j]) for i in range(mm))]
    c = 0
    for D in DS:
        for U in US:
            if D & U: continue
            c += 1
    nInst += c
    if c >= 2: nQmulti += 1
print(f"n={n}  Q（width<=2, {mm} 元）: {nQ}")
print(f"      (Q,D,U) 实例: {nInst}")
print(f"      盒可多切割的 Q: {nQmulti}")
print(f"耗时 {time.time()-t0:.0f}s")
