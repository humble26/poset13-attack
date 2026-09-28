# 偏序集 1/3–2/3 猜想 width-3 情形攻击：项目档案

会话日期：2026-09-13 至 2026-09-27。本文件夹保存该次攻击的全部程序、数据与文档。

## 1. 问题

对有限偏序集 P，e(P) 为线性扩展数，p_P(x,y) = e(P+x<y)/e(P)（x∥y），δ(P) = max_{x∥y} min(p, 1−p)。
**猜想（Kislitsyn 1968）**：P 非全序 ⟹ δ(P) ≥ 1/3。
最强一般结果：Kahn–Saks 1984（δ ≥ 3/11）、Brightwell–Felsner–Trotter 1995（δ ≥ 1/2−√5/10，无穷情形不可再改进）。
本会话主攻 **width-3 情形**：该类无任何一般性结果（height-2 已由 Trotter–Gehrlein–Fishburn 1992 解决；
width-2 由 Linial 1984 解决、Sah 2021 加强），且 Saks 1985 的极值例给出 14/39 > 1/3 的余量。

## 2. 数学结果（完整证明见第 4 节）

| 编号 | 内容 | 状态 |
|---|---|---|
| 引理 0.1–0.4 | 扩展存在性 / 相邻交换 / 前缀下集 / p 的组合刻画 | 已证（初等） |
| 引理 A | e(P) = Σ_{z∈max(P)} e(P∖z)，及 p 的顶部递归（推论 A′） | 已证 |
| 引理 B | 临界对 (x,y)（D(x)⊆D(y), U(x)⊆U(y)）必被某线性扩展紧贴覆盖 | 已证 |
| 引理 C | 插入权重恒等式：p_P(x,y) = Σ_M i_z(M)·[x<_M y] / Σ_M i_z(M) | 已证 |
| 命题 1 | z 全可比单点扩张：p 值逐点不变，继承 1/3–2/3；特例 width(P∖z)≤2 + z 全可比情形 | 已证 |
| 目标定理 | width(P∖z) ≤ 2（z 任意）⟹ δ(P) ≥ 1/3 | **开放（本会话前沿）** |

关键计算发现（n ≤ 8 数值侦察）：
- 目标类不比一般 width-3 容易：n=6 全查 min δ = 4/11；n=7 采样 4/11；n=8 采样 **14/39**（Saks 极值常数）。
- "平衡稳定性"路线被否定：插入 z 可把 p ≥ 1/2 的对压到 2/9。
- 但出现**平衡守恒**现象：879 例中 δ(P) ≥ 4/11 从不失守；19 个"全部宽度 2 平衡对被毁"的例子里 δ(P) ∈ [4/11, 18/37]。
- 下一步：把守恒写成计数恒等式并证明（引理 C 已把问题化归为权重 i_z(M) 的重新分配；z 极大时 i_z(M) = |M|+1−t，t 为 D(z) 在 M 中的最大位次）。

## 3. 形式化检验（Lean 4.33.1）——最终结果：全部通过 ✅

文件：`lean/Poset13.lean`（自足单文件）+ `lean/Poset13Defs.lean`（定义模块）+ `lean/t_*.lean`（分定理验证文件）。
不依赖 mathlib。

### 验证记分板（最终，见 r_*.log）

| 定理 | 范围 | 结果 | 公理披露 |
|---|---|---|---|
| `chain_isPoset` | 链实例 | ✅ | **无任何公理**（纯核内 decide） |
| `chain_e` / `anti_e` | 实例 e=1 / e=24 | ✅ | [propext] |
| `checks4_true` | 全部 2^6 = 64 个 n=4 候选上 A+B+猜想+C+P1 全查 | ✅ | [propext, Lean.ofReduceBool] |
| `cFires4` / `pFires4` | 非空性：C/P1 防卫条件确实在真实偏序集上触发 | ✅ | [propext, Lean.ofReduceBool] |
| `checks5_true` | 全部 2^10 = 1024 个 n=5 候选上五类全查 | ✅ | [propext, Lean.ofReduceBool] |
| `checks6AB_true` | 全部 2^15 = 32768 个 n=6 候选上 A+B 全查 | ✅ | [propext, Lean.ofReduceBool] |

（`Lean.ofReduceBool` 是 `native_decide` 的标准信任公理：声明"编译器求值的 Bool 结果可信"。）

### 被形式化过程抓到的实现错误（人眼复查全部漏过，共 4 + 1 个）

1. `critPair` 未排除 x = y（Python 版 combinations 天然不触发）；
2. `critPair` 下/上界布尔公式笔误；
3. **`skipIdx`/`downIdx` 重编号方向写反**（原标签→Q 标签应为 i−1 而非 i+1），引理 C 检验在反链上即失败——用 Python 镜像交叉定位；
4. **`checkConj` 缺少"非链"前提**：链没有不可比对，`any` 返回 false，而猜想对链空洞真。n=4 的 8 个前向标号链表示全部触发假失败，Python 镜像与 Lean 报出完全相同的失败 mask 集合，互相印证；
5. （性能）`addRel` 的嵌套闭包链无记忆化，单次 `lt` 探测 O((n+1)^fuel) 指数展开——改为物化查找表后从小时级降到分钟级，语义不变。

### 定位方法论

`native_decide` 编译执行快但失败只报真值不报位置；`#eval` 报位置但解释器慢 3 个数量级。
有效组合：Python 精确镜像 Lean 检验器语义快速定位 → Lean `native_decide` 快速确认 →
轻量子集文件先行验证修复 → 分定理独立文件并行验证（结果落各自日志，抗中断）。

## 4. 证明全文

### 记号
(P, ≺) 有限严格偏序；x∥y 不可比；线性扩展、e(P)、P+x≺y（加入关系取传递闭包）、
p_P(x,y) = e(P+x≺y)/e(P)、δ(P) = max_{x∥y} min(p, 1−p)。

### 引理 0.1（扩展存在性）
有限偏序必有极小元（无穷递降链与有限性矛盾）；取极小元 m，归纳排 P∖{m}，m 置首。∎

### 引理 0.2（相邻交换）
扩展中相邻的不可比元素可交换仍为扩展：逐一验证涉及 a,b 的五种关系位置情形。∎

### 引理 0.3（前缀是下集）
u 在前缀、v ≺ u ⟹ v 在前缀（v 必须排在 u 前）。∎

### 引理 0.4（p 的组合刻画）
ext(P+x≺y) = {L ∈ ext(P) : x 在 y 前}；P+x≺y 的每个关系要么是旧关系、要么是 x≺y、
要么由传递性导出且 L 满足之。∎

### 引理 A
扩展的最后元素必为极大元；按最后元素分划 ext(P)，"末尾追加 z"给出 ext(P∖z) ≅ E_z 的双射。∎

### 推论 A′（顶部递归）
p_P(x,y) = q_y + Σ_{w∈max(P)∖{x,y}} (e(P∖w)/e(P))·p_{P∖w}(x,y)，q_y = e(P∖y)/e(P)。
（按最后元素 w = y / w = x / 其它三种情形分块计数。）∎

### 引理 B（临界对覆盖）
临界对 (x,y)：x∥y 且 D(x)⊆D(y)、U(x)⊆U(y)。在 P′=P+x≺y 中取扩展 L′（x 在 y 前）。
夹层元素 s：s≺x 排除（s 须在 x 前）；x≺s 排除（U(x)⊆U(y) ⟹ y≺s ⟹ s 在 y 后）；y≺s 排除。
故 s∥x（P′ 中），用引理 0.2 逐步上移 x 贴住 y。每步均为 P（从而 P′）的扩展。∎

### 引理 C（插入权重恒等式）
i_z(M) = #{k : D(z) ⊆ prefix_k(M) 且 U(z) ⊆ suffix_k(M)}。
L ↦ (L∖z, z 的位置) 给出 ext(P) 与 {(M,k) : M ∈ ext(Q), k 合法空隙} 的双射
（合法空隙恰是使插入结果满足全部关系者）；删 z 不改 x,y 相对次序。
故 e(P) = Σ_M i_z(M)，e(P+x≺y) = Σ_M i_z(M)·[x<_M y]，相除即得。∎

注记：z 极大时 U(z)=∅，i_z(M) = |M|+1−t（t = D(z) 在 M 中的最大位次）。

### 命题 1
z 与其余一切元素可比 ⟹ 合法空隙唯一（k = |D(z)|），i_z(M) = 1；
代入引理 C 得 p_P(x,y) = p_{P∖z}(x,y)；width ≤2 母图情形由 Linial 1984 / Sah 2021 承接。∎

### 目标定理（开放）
width(P∖z) ≤ 2 ⟹ δ(P) ≥ 1/3。剩余步骤：证明"平衡守恒"（见第 2 节末条）。

## 5. 复现指南

Python（无需依赖）：
```
cd python && python poset_experiments.py all     # 或 e1 / e2 / e3 / e4
```

Lean（工具链 v4.33.1，GitHub 直连不可用时用镜像
`https://ghproxy.net/https://github.com/leanprover/lean4/releases/download/v4.33.1/lean-4.33.1-windows.zip`，
解压后 `PATH=$HOME/lean4/bin:$PATH`）：
```
cd lean
lean -o Poset13Defs.olean Poset13Defs.lean          # 构建定义模块
LEAN_PATH=. lean t_kernel.lean                      # 核内锚点（秒级）
LEAN_PATH=. lean t_checks4.lean                     # n=4 全查（约 1 分钟）
LEAN_PATH=. lean t_fires4.lean                      # 非空性（约 1 分钟）
LEAN_PATH=. lean t_checks5.lean                     # n=5 全查（约 20 分钟）
LEAN_PATH=. lean t_checks6AB.lean                   # n=6 A+B（约 30–60 分钟）
# 每个文件打印 #print axioms 披露信任等级；预期全部 exit=0
```
单文件版 `Poset13.lean` 含全部定义与定理，可 `lean Poset13.lean` 直跑（较慢）。

## 7. 第二阶段：守恒机制与理想分解（2026-09-27 会话）

### 已证新结果

**引理 D（理想分解）**：e(P) = Σ_{I} e(I)·e(Q∖I)，I 取遍 Q = P∖z 的理想且 D(z) ⊆ I、I ∩ U(z) = ∅。
证明：L ↦ (before-z 集 I, 前缀扩展, 后缀扩展) 是双射（理想的下封闭性排除 u∉I<v 型矛盾）。
Python 验证 300/300；Lean 验证至 n=5（`t_lemmad.lean`，ld4/ld5 全过，axioms 干净）。

**条带公式**：对 u ∥ z，p_P(z,u) = Σ_{I ∈ box, u ∉ I} e(I)e(Q∖I) / e(P) —— z-对的 p 值恰是
"盒子内避开 u 的理想质量占比"。Python 验证 267/267；Lean 同上验证。

**引理 E（多数传递性 / 等价形式）**：p(x,z) ≥ p(x,y) + p(y,z) − 1 恒成立（并集界）。
故若所有不可比对都不平衡，则"多数关系" χ(x,y): p(x,y) > 2/3 可传——
**目标定理 ⟺ P = Q+z（width Q ≤ 2）不存在"2/3 多数线性扩张"**。

**角点对下界（并集界推论）**：x "低"（x ∈ D 或 x ∥ z 且 p(z,x) < 1/3）且 y "高"
（y ∈ U 或 y ∥ z 且 p(z,y) > 2/3）⟹ p_P(x,y) ≥ P(x∈I, y∉I) ≥ 1/3。

### 数值发现（E5–E8）

- E5（n=5,6 穷举 126,159 例）：目标定理零反例（min δ = 1/3，由 width-2 的 P 达到，与 Aigner 猜想相容）；
  **守恒二择一零例外**：z-对全失败（35,753 例，28%）时必有 Q-对平衡（both_fail = 0）。
- E6（n=7,8,9 各 1500 例采样）：target_fail = 0，both_fail = 0。（注意：P 为链时猜想空洞，
  分类须排除——第一版采样把链误计为失败。）
- E7：z-对全失败 ⟹ 存在质量 > 1/3 的"重层"（沿链的层状累积跳过 [1/3,2/3]）；重层质量最小观测值 4/11。
- E8：在 245 个可评估案例中，**角点对 (a_{i*+1}, b_{j*+1}) 全部平衡（0 反例）**；其余案例因切割
  落在链端点等结构原因需推广定义（失败族 2 显示救援对可含 D 中元素，如 (最后 D 元素, 孤立元) 对 p=2/5）。

### 当前缺口（下一步）

并集界只给单侧界：角点对 p ≥ 1/3 已证，但 p ≤ 2/3 未证——上界需要 e(I)e(Q∖I) 权重的精细结构
（前缀/后缀扩展的内部序统计）。候选路线：(a) 多数序 ≡* 的 z-切割 + 单调 p-矩阵阶梯 + 重层量化矛盾；
(b) 把 1/3 放宽为 14/39 寻找余量；(c) 对 corridor 上 e(i,j) 的 Pascal 递推做精确恒等式。

### 形式化新抓到的 bug（第 6 个）

`t_lemmad` 的 lemmaDAll 守卫方向写反（`!(可比) || strip` 变成"可比才检查"）——
用 #eval 逐层定位：countSub ✓ → lemmaDCheck ✓ → stripCheck（u==z 自对泄露）✓ → 守卫方向 ✗。

### 新增文件

`lean/t_lemmad.lean`、`lean/r_lemmad.log`；`python/step2/`（e5_ideal / e5_classify / e6_sample /
e7_heavy / e7_debug / e8_corner）。

## 8. 文件清单

```
lean/Poset13.lean        Lean 4 完整单文件（定义 + 检验器 + 穷举定理）
lean/Poset13Defs.lean    定义模块（供分定理验证文件 import）
lean/t_kernel.lean       核内锚点验证（chain/antichain 实例）
lean/t_checks4.lean      n=4 全查（五类检查 × 64 候选）
lean/t_fires4.lean       C/P1 防卫条件非空性
lean/t_checks5.lean      n=5 全查（五类检查 × 1024 候选）
lean/t_checks6AB.lean    n=6 引理 A+B（32768 候选）
lean/t_lemmad.lean       引理 D + 条带公式（n=4,5 全查）
lean/r_*.log             每个验证文件的最终运行日志（含 #print axioms）
lean/lakefile.lean       Lake 构建配置
lean/lean-toolchain      工具链版本钉
python/poset_experiments.py  四个数值实验（会话输出数值已注明）
python/mirror4.py        Lean 检验器的 Python 精确镜像（bug 定位工具）
python/step2/            第二阶段实验（E5–E8）
README.md                本文档
```
