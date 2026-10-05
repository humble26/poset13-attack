"""verify_all.py — 全量验证主运行器（三层）

Layer L: Lean 形式化全套（并行，独立日志目录 verify_lean/）
Layer P: Python 验证电池（既有 e-脚本逐个子进程运行 + 期望标记检查）
Layer X: 交叉一致性（关键数字与 paper.md 声明对照）

用法: python verify_all.py [lean|py|all] （默认 all）
每项输出 PASS/FAIL 与关键计数；最终汇总表。退出码 = 失败数。
"""
import subprocess, sys, os, re, shutil, time, glob
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

ROOT = Path(__file__).resolve().parent
LEAN = ROOT / "lean"
PY3 = sys.executable
ENV = {**os.environ, "PATH": os.path.expanduser("~/lean4/bin") + os.pathsep + os.environ.get("PATH", ""),
       "LEAN_PATH": str(LEAN)}

results = []


def record(layer, name, ok, detail=""):
    results.append((layer, name, ok, detail))
    print(("[PASS] " if ok else "[FAIL] ") + f"{layer}/{name}  {detail}")


def run_cmd(cmd, cwd, timeout, env=None):
    try:
        p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True,
                           timeout=timeout, env=env)
        return p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired:
        return -1, "", "TIMEOUT"


# ---------------- Layer L: Lean ----------------
LEAN_JOBS = (["t_kernel.lean", "t_checks4.lean", "t_fires4.lean", "t_lemmad.lean",
              "t_checks5.lean", "t_checks6AB.lean", "t_s4.lean"] +
             [f"t_s5_{k}.lean" for k in range(8)] +
             [f"t_s6_{k}.lean" for k in range(16)])


def lean_one(fname):
    out = LEAN / "verify_lean" / (fname.replace(".lean", ".log"))
    rc, so, se = run_cmd(["lean", fname], str(LEAN), 5400, ENV)
    out.write_text(so + se, encoding="utf-8")
    text = so + se
    ok = rc == 0 and "error" not in text.lower() and ("depends on axioms" in text or "#print" not in
                                                     open(LEAN / fname, encoding="utf-8").read())
    return fname, ok, f"rc={rc}"


def run_lean():
    vdir = LEAN / "verify_lean"
    vdir.mkdir(exist_ok=True)
    # 先确认 olean 存在
    for o in ("Poset13Defs.olean", "Poset13Step3.olean"):
        if not (LEAN / o).exists():
            record("L", f"build {o}", False, "olean 缺失")
            return
    record("L", "oleans", True, "Poset13Defs/Poset13Step3 已构建")
    with ThreadPoolExecutor(max_workers=30) as ex:
        for fname, ok, detail in ex.map(lean_one, LEAN_JOBS):
            record("L", fname, ok, detail)


# ---------------- Layer P: Python 电池 ----------------
PY_CHECKS = [
    # (name, script, args, must_contain[list], must_not_contain[list], timeout)
    ("E19 恒等式", "e19_identities.py", [],
     ["delta identity: ok=", "m-frame reweighting: ok="], ["bad=1", "bad=2", "MISMATCH"], 1200),
    ("E23 自由走廊三事实", "e23_corner.py", [],
     ["corner ok/bad: 16/0", "vandermonde ok/bad: 72/0"], ["FAIL"], 900),
    ("E26 定理7.4 构造", "e26_freelinial.py", [],
     ["E26: constructed-proof verification ok=64 bad=0"], ["FAIL"], 900),
    ("E35 z-孤立归约", "e35_ziso.py", [],
     ["z-isolated p_P=p_Q: ok=98 bad=0"], ["ZISO FAIL"], 900),
    ("E29 边-顶点恒等式", "e29_edge.py", [],
     ["edge/vertex identities: ok=100 bad=0"], ["FAIL"], 900),
    ("E43 链盒闭式", "e43_chainbox_closed.py", [],
     ["closed-form checks: ok="], ["FAIL: ('e("], 1800),
    ("E41 扫过LP不可行", "e41_sweeplp.py", [],
     ["ALL INFEASIBLE"], ["FEASIBLE sweep-failure"], 2400),
    ("E36c 同侧存在性", "e36c_dexist.py", [],
     ["'d_bad': 0", "'u_bad': 0"], [], 2400),
    ("e_9e_chainbox 计数公式", "e_9e_chainbox.py", [],
     ["全部通过 ✅"], ["❌"], 600),
    ("crosscheck_9de 独立对拍", "crosscheck_9de.py", [],
     ["全部一致 ✅"], ["❌"], 1200),
    ("E47 (lo,hi)三项+着陆", "e47_lohi.py", [],
     ["'landed': 1823", "'race_bad': 0", "'notlanded': 0"], ["NOT-LANDED"], 3400),
]


def run_py():
    for (name, script, args, must, mustnot, timeout) in PY_CHECKS:
        path = ROOT / "python" / "step3" / script
        if not path.exists():
            record("P", name, False, "脚本缺失")
            continue
        rc, so, se = run_cmd([PY3, script] + args, str(path.parent), timeout)
        text = so + se
        ok = rc == 0 and all(m in text for m in must) and not any(m in text for m in mustnot)
        detail = ""
        for m in must:
            mm = re.search(re.escape(m) + r"[^|\n]*", text)
            if mm: detail += mm.group(0)[:60] + " ; "
        record("P", name, ok, detail[:120])


# ---------------- Layer X: 交叉一致性 ----------------
def run_x():
    # X1: paper.md 存在且包含关键编号
    p = ROOT / "paper.md"
    ok = p.exists()
    txt = p.read_text(encoding="utf-8") if ok else ""
    keys = ["定理 7.4", "命题 9.A", "命题 9.B", "引理 9.C", "引理 9.D", "引理 9.E",
            "引理 9.F", "猜想 6.1", "定理 6.2", "命题 5.5", "e(Q∖D) = (β₁+1)(β₂+1)"]
    missing = [k for k in keys if k not in txt]
    record("X", "paper.md 编号完整", ok and not missing, "缺失: " + str(missing) if missing else "全部在案")
    # X2: 修正标注在位（E43 纠错不被遗忘）
    record("X", "9.9 修正标注在位", "E43 机检纠错" in txt, "")
    # X3: 死路清单条数
    dd = len(re.findall(r"否证|反例", txt))
    record("X", "死路/否证记录条数(>=13)", dd >= 13, f"count={dd}")


def main():
    phase = sys.argv[1] if len(sys.argv) > 1 else "all"
    t0 = time.time()
    if phase in ("lean", "all"): run_lean()
    if phase in ("py", "all"): run_py()
    if phase in ("x", "all"): run_x()
    print("=" * 60)
    npass = sum(1 for r in results if r[2])
    print(f"汇总: {npass}/{len(results)} 通过  (耗时 {time.time()-t0:.0f}s)")
    for (layer, name, ok, detail) in results:
        if not ok: print("  FAILED:", layer, name, detail)
    sys.exit(0 if npass == len(results) else 1)


if __name__ == '__main__':
    main()
