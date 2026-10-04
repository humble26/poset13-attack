"""估算 n=7 全穷举的规模：
Q 有 6 个元素，枚举 2^15 个 mask（固定"恒等排列是线性扩展"的标号，WLOG），
筛 poset + width<=2，再枚举 (D,U)。统计实例总数。"""
import itertools, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from my_poset import my_close, my_is_poset, my_width


def main():
    m = 6
    npairs = [(i, j) for i in range(m) for j in range(i + 1, m)]
    t0 = time.time()
    nQ = 0; nInst = 0
    nQ_multi = 0
    for mask in range(1 << len(npairs)):
        rel = [npairs[k] for k in range(len(npairs)) if (mask >> k) & 1]
        Q = my_close(m, rel)
        if not my_is_poset(Q, m): continue
        if my_width(Q, m) > 2: continue
        nQ += 1
        DS = [S for S in range(1 << m)
              if all(not ((S >> i) & 1) or all((S >> j) & 1 for j in range(m) if Q[j][i]) for i in range(m))]
        US = [S for S in range(1 << m)
              if all(not ((S >> i) & 1) or all((S >> j) & 1 for j in range(m) if Q[i][j]) for i in range(m))]
        cnt = 0
        for D in DS:
            for U in US:
                if D & U: continue
                cnt += 1
        nInst += cnt
        if cnt >= 2: nQ_multi += 1
    print(f"Q 数（width<=2, 6 元）: {nQ}")
    print(f"(Q,D,U) 实例总数: {nInst}")
    print(f"其中盒非单点（|box|>=2）的 Q 占比: {nQ_multi}/{nQ}")
    print(f"耗时 {time.time()-t0:.1f}s")


if __name__ == '__main__':
    main()
