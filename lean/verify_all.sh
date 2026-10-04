#!/bin/bash
# ============================================================================
# poset13-attack / Lean 4 验证驱动器
#
# 用法：
#   bash verify_all.sh light     # 快速组：核内锚点 + n=4 全查（约 10 秒）
#   bash verify_all.sh core      # 常规组：light + 引理 D + n=5 全查（约 20-40 分钟）
#   bash verify_all.sh full      # 完整组：core + n=6 全查 + n=6 AB（32 核并行，约 1 小时）
#   bash verify_all.sh all
#
# 关键点（README 旧版漏了这一步）：
#   编译模块与运行验证文件**都必须带 LEAN_PATH=.**，
#   否则报 unknown module prefix 'Poset13Defs'。
#
# 每个验证文件独立日志 r_<name>.log，末尾记录 **退出码**。
# 旧日志缺退出码，无法区分"跑过"与"跑挂了"—— 本脚本补上这一信息。
# ============================================================================
set -u
cd "$(dirname "$0")" || exit 1

# --- 定位 Lean 工具链 -------------------------------------------------------
# 优先级：本地已装 elan/lake > 常见的解压式安装目录 > PATH
find_lean() {
  if command -v lean >/dev/null 2>&1; then echo "lean"; return 0; fi
  for d in "$HOME/lean4/bin" "$HOME/.elan/bin" "/c/lean4/bin"; do
    [ -x "$d/lean.exe" ] && { echo "$d"; return 0; }
    [ -x "$d/lean" ]     && { echo "$d"; return 0; }
  done
  return 1
}
LEAN_BIN="$(find_lean)" || {
  echo "!! 找不到 lean 可执行文件。安装说明见 README 第 5 节。" >&2
  exit 127
}
export PATH="$LEAN_BIN:$PATH"
echo "工具链：$(lean --version)"
echo "日志目录：$(pwd)"
echo

# --- 构建定义模块（带 LEAN_PATH=.）-----------------------------------------
echo "=== 构建模块 ==="
for m in Poset13Defs Poset13Step3; do
  if LEAN_PATH=. lean -o "$m.olean" "$m.lean" > "/tmp/${m}_build.log" 2>&1; then
    echo "  $m.olean  构建 OK"
  else
    echo "  $m.olean  构建失败："; sed 's/^/    /' "/tmp/${m}_build.log" | head -10
    exit 1
  fi
done
echo

# --- 分组定义 ---------------------------------------------------------------
# 每项：<日志名> <验证文件> <并行度组>
LIGHT=(kernel s4 checks4 fires4)
CORE=(lemmad checks5 "${LIGHT[@]}")
# s5 分片：8 片，可并行
S5_SHARDS=("s5_0" "s5_1" "s5_2" "s5_3" "s5_4" "s5_5" "s5_6" "s5_7")
# s6 分片：16 片，可并行（单片 2.5-4.5 分钟）
S6_SHARDS=("s6_0" "s6_1" "s6_2" "s6_3" "s6_4" "s6_5" "s6_6" "s6_7" \
           "s6_8" "s6_9" "s6_10" "s6_11" "s6_12" "s6_13" "s6_14" "s6_15")

# --- 运行器 -----------------------------------------------------------------
FAILED=0
run_one() {
  local name="$1" file="t_${1}.lean" log="r_${1}.log"
  local t0 t1 dt
  t0=$(date +%s)
  LEAN_PATH=. lean "$file" > "/tmp/${name}_run.out" 2>&1
  local rc=$?
  t1=$(date +%s); dt=$((t1 - t0))
  {
    cat "/tmp/${name}_run.out"
    echo "exit=$rc"
    echo "elapsed=${dt}s"
  } > "$log"
  if [ $rc -ne 0 ]; then
    echo "  [FAIL] $name  exit=$rc  (${dt}s)  -> $log"
    FAILED=1
  else
    echo "  [ok]   $name  exit=0  (${dt}s)"
  fi
}
export -f run_one
export PATH LEAN_BIN

run_serial() {
  for n in "$@"; do run_one "$n"; done
}

# 并行跑一组分片（每个分片一个进程，自身是独立的 native_decide）
run_parallel() {
  local pids=() names=()
  for n in "$@"; do
    ( run_one "$n" ) &
    pids+=($!); names+=("$n")
  done
  for i in "${!pids[@]}"; do
    wait "${pids[$i]}" || true
  done
}

MODE="${1:-all}"
echo "=== 验证模式：$MODE ==="
echo

case "$MODE" in
  light)
    echo "--- 快速组 ---"
    run_serial "${LIGHT[@]}"
    ;;
  core)
    echo "--- 快速组 ---"; run_serial "${LIGHT[@]}"
    echo "--- 引理 D ---"; run_one lemmad
    echo "--- n=5 五类检查 ---"; run_one checks5
    echo "--- n=5 全查（8 分片并行）---"
    run_parallel "${S5_SHARDS[@]}"
    ;;
  full)
    echo "--- 快速组 ---"; run_serial "${LIGHT[@]}"
    echo "--- 引理 D ---"; run_one lemmad
    echo "--- n=5 五类检查 ---"; run_one checks5
    echo "--- n=5 全查（8 分片并行）---"; run_parallel "${S5_SHARDS[@]}"
    echo "--- n=6 全查（16 分片并行）---";  run_parallel "${S6_SHARDS[@]}"
    echo "--- n=6 引理 A+B ---"; run_one checks6AB
    ;;
  all) "$0" full ;;
  *) echo "未知模式：$MODE（可选 light / core / full / all）" >&2; exit 2 ;;
esac

echo
echo "=== 汇总 ==="
printf "%-14s %-8s %s\n" "验证文件" "退出码" "耗时"
for f in r_*.log; do
  rc=$(grep -o 'exit=[0-9]*' "$f" 2>/dev/null | tail -1 | cut -d= -f2)
  el=$(grep -o 'elapsed=[0-9]*s' "$f" 2>/dev/null | tail -1 | cut -d= -f2)
  [ -z "$rc" ] && rc="(无记录)"
  [ -z "$el" ] && el="-"
  printf "%-14s %-8s %s\n" "${f#r_}" "$rc" "$el"
done
echo
[ $FAILED -eq 0 ] && echo "全部通过 ✅" || echo "存在失败 ❌"
exit $FAILED
