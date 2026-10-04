#!/bin/bash
# n=8 全穷举：把 [0, 2^21) 切成 32 片并行跑，汇总结果。
cd "$(dirname "$0")" || exit 1
PY="python"
OUT=/tmp/exh_n8
rm -rf "$OUT"; mkdir -p "$OUT"
TOTAL=$((1<<21))
SHARDS=32
STEP=$((TOTAL / SHARDS))
pids=()
for i in $(seq 0 $((SHARDS-1))); do
  lo=$((i * STEP))
  hi=$(( lo + STEP ))
  ( "$PY" exhaustive_n.py 8 "$lo" "$hi" > "$OUT/s$i.txt" 2>&1 ) &
  pids+=($!)
done
echo "已启动 $SHARDS 个分片，等待..."
for p in "${pids[@]}"; do wait "$p"; done
echo "=== 汇总 ==="
python - << 'PY'
import glob, ast, collections
tot = collections.Counter()
bad = []
for f in sorted(glob.glob("/tmp/exh_n8/s*.txt")):
    for line in open(f, encoding="utf-8"):
        if line.startswith("SUMMARY|"):
            _, lo, hi, d = line.strip().split("|", 3)
            for k, v in ast.literal_eval(d).items():
                tot[k] += v
        if "反例" in line and "零反例" not in line:
            bad.append((f, line.strip()))
print("n=8 全穷举汇总:")
for k in sorted(tot): print(f"   {k:18s} {tot[k]}")
print("反例:", bad if bad else "零反例 ✅")
PY
