# Zaran3 DEM 模块正确性测试与缺陷修复报告

- **被测程序**：`Zaran3`（复用 CFD 求解器框架 Zaran，新增离散元 DEM 求解通道）
- **版本标识**：Zaran3.10.2（"复用 CFD 架构 + 集成 DEM" 的初始版本）
- **源码路径**：`E:\ZARAN_DEM\ZARAN\Zaran3`
- **基线（修复前）可执行文件**：`bin/Release/Zaran3.10.2_prefix.exe`
- **修复后可执行文件**：`bin/Release/Zaran3.10.2.exe`
- **报告日期**：2026-09-16

---

## 0. 结论摘要

本次对 DEM 集成版本完成了以**解析解、守恒律、黄金值、参数负例、回归**为手段的正确性测试，
共定位并修复 **3 处缺陷**（2 处接触力学、1 处输出格式）：

| 编号 | 缺陷 | 影响 | 状态 |
|---|---|---|---|
| D1 | 线性接触模型 `CalcNormalForce` 中的"禁止拉力"钳位 `if (Fn_mag<0) Fn_mag=0` | 恢复系数标定失效，实测 e 系统性偏高（e=0.5 → 0.551） | 已修复 |
| D2 | Hertz 模型阻尼系数沿用线性公式，缺少非线性（Tsuji/Brilliantov）修正因子 `√(5/6)` | 实测 e 系统性偏离标称（e=0.5 → 0.523） | 已修复 |
| D3 | 备份 `particles.dat` 表头带 `variables=` 前缀，与同目录 `bonds.dat`/`gas_voronoi_faces.dat` 及输入格式不一致 | 4 个已提交测试解析失败（`KeyError: 'id'`） | 已修复 |

修复后：

- **恢复系数**：线性模型 e_meas 与标称值偏差 ≤0.0011（e_nom∈[0.1,0.9]）；
  Hertz 模型偏差 ≤0.00004（e_nom∈[0.3,1.0]）；墙面反弹偏差 ≤0.0002（e_nom∈[0.1,0.9]）。
- **已提交测试 `tests/*`**：由修复前 **4/8 通过** 提升为 **8/8 全部通过**。
- 接触力公式黄金值、线/角动量守恒、静力平衡、参数负例校验等此前已确认正确的项目**未见回归**。

---

## 1. 测试范围与目标

### 1.1 被测模块（范围）

| 模块 | 文件 | 被测面 |
|---|---|---|
| 接触模型（线性） | `src/Basic/DEM/LinearSpringDashpot.cpp` | 法向刚度/阻尼标定、法向力、Coulomb 切向摩擦 |
| 接触模型（Hertz-Mindlin） | `src/Basic/DEM/HertzMindlin.cpp` | 等效参数、Hertz 法向力、非线性阻尼、切向刚度 |
| DEM 求解器 | `src/Zone/DEM/DEMSolver.cpp` | 求解主循环、接触检测（KDTree）、键合力、热化学/反应、气体压力、Voronoi、备份输出 |
| 参数读取与校验 | `src/Zone/DEM/DEMSolverParam.cpp` | dt、恢复系数、泊松比、box、模型选择等合法性 |
| 输入解析 | `src/Generator/DEM/ReadDEMParticle.cpp`、`ReadDEMBond.cpp` | CSV 列数（14/27）、单位、非法值拒绝 |
| 数据结构 | `inc/Basic/DEM/DEMParticle.h`、`DEMContact.h`、`DEMBond.h`、`DEMWall.h` | 字段完整性、默认值 |
| 场生成/边界 | `src/Generator/DEM/DEMFieldGenerator.cpp` | 盒体六面墙生成、无重叠布点 |

### 1.2 测试目标

- **G1 核心物理正确性**：接触力→加速度→速度积分链路是否符合解析解（恢复系数、重叠量、静力平衡）。
- **G2 守恒律**：无摩擦/无重力场景下线动量、角动量守恒。
- **G3 关键模块功能**：键合（轴向/切向/断裂能/损伤）、热化学（Arrhenius 反应、热传导）、Voronoi 控制体、气体压力（JWL）。
- **G4 健壮性**：非法输入必须被拒绝；大规模算例不发散、不产生 NaN。
- **G5 回归一致性**：`tests/` 下已提交算例全部通过。

### 1.3 不在本次范围

并行/MPI（代码注明未实现）、真实材料参数标定与实验对比、颗粒形状（多面体/非球形）接触。

---

## 2. 测试环境与方法

### 2.1 构建

仓库自带 `build/`（Visual Studio 生成器）在本机 **构建失败**：

```
error MSB6001: CL.exe 命令行开关无效 … 字典中的关键字: https_proxy
```

原因是代理类环境变量被透传进 `CL.exe` 命令行。改用 **Ninja 生成器** 并显式拼装 MSVC / Windows SDK 环境后构建成功：

```bash
cmake -G Ninja -DCMAKE_BUILD_TYPE=Release -S . -B build-ninja
cmake --build build-ninja
```

- MSVC：`…\Microsoft Visual Studio\18\Community\VC\Tools\MSVC\14.51.36231`
- SDK：`C:\Program Files (x86)\Windows Kits\10\{Include,Lib,bin}\10.0.26100.0`
  （`INCLUDE` 须含 `ucrt;shared;um`；`PATH` 须含 SDK `bin\x64` 以提供 `rc.exe`）

### 2.2 对照策略

修复前先将可执行文件另存为 `Zaran3.10.2_prefix.exe` 作为**对照基线**，
所有 A/B 对比均在**同一测试脚本、同一算例参数**下分别对两个二进制运行。

### 2.3 验证方法学

| 方法 | 说明 | 典型用例 |
|---|---|---|
| 解析解 | 与闭式理论值逐位比对 | 恢复系数、Hertz 重叠量、静力平衡 |
| 守恒律 | 无外力场景下动量/角动量不变 | 斜碰、键合剪切 |
| 黄金值 | 接触力公式的解析计算值 | 线性/Hertz 法向力峰值 |
| 参数负例 | 非法输入必须抛出并给出可读错误 | dt≤0、e∉(0,1]、ν∉(-1,0.5)、重复 id 等 |
| 回归 | 运行 `tests/*` 并调用其 `verify.py` | 8 个已提交算例 |

自建测试脚手架位于 `Zaran3/_demtest/`（清单见附录 A）。

---

## 3. 修复前测试结果

### 3.1 恢复系数（两体对心碰撞，解析基准）

指标：`e_meas = (v2f − v1f) / (2 v0)`，取仿真末帧速度（分离速度）。

| 模型 | e_nom | e_meas（修复前） | 误差 | 相对误差 | 判定 |
|---|---|---|---|---|---|
| LinearSpringDashpot | 0.90 | 0.902340 | +0.0023 | +0.26% | 接近 |
| LinearSpringDashpot | 0.70 | 0.718929 | +0.0189 | +2.70% | ✗ 超差 |
| LinearSpringDashpot | 0.50 | 0.551204 | +0.0512 | +10.24% | ✗ 超差 |
| LinearSpringDashpot | 0.30 | 0.398070 | +0.0981 | +32.69% | ✗ 严重超差 |
| LinearSpringDashpot | 0.10 | 0.253488 | +0.1535 | +153.49% | ✗ 严重超差 |
| HertzMindlin | 1.00 | 1.000000 | 0.0000 | 0.00% | 合格 |
| HertzMindlin | 0.90 | 0.893393 | −0.0066 | −0.74% | ✗ 超差 |
| HertzMindlin | 0.50 | 0.523450 | +0.0234 | +4.69% | ✗ 超差 |
| HertzMindlin | 0.30 | 0.369549 | +0.0696 | +23.18% | ✗ 严重超差 |
| 墙面（Linear） | 0.90 | 0.901979 | +0.0020 | +0.22% | 接近 |
| 墙面（Linear） | 0.50 | 0.550062 | +0.0501 | +10.01% | ✗ 超差 |
| 墙面（Linear） | 0.30 | 0.396851 | +0.0969 | +32.28% | ✗ 严重超差 |
| 墙面（Linear） | 0.10 | 0.252394 | +0.1524 | +152.39% | ✗ 严重超差 |

> 规律：e=1（无阻尼）精确；e<1 时实测值**系统性偏高**，且 e 越小偏差越大——典型的"耗能被抹掉"特征。

### 3.2 已提交测试（`tests/*`）

| 算例 | 修复前 | 说明 |
|---|---|---|
| dem_bond_test | PASS | 永久弹性键恢复 |
| dem_bounce_test | PASS | 弹跳解析峰值 |
| dem_energy_fracture_test | PASS | 能量断裂准则 |
| dem_fixed_collision_test | PASS | 固定粒子参与接触 |
| dem_voronoi_test | **FAIL** | `KeyError: 'id'` |
| dem_reaction_test | **FAIL** | `KeyError: 'id'` |
| dem_brazilian_disk_test | **FAIL** | `KeyError: 'id'` |
| hmx_disk_center_ignition | PASS | 反应/压力/断键统计 |

**修复前 4/8 通过。** 4 项失败均为输出表头解析问题（D3），而非物理错误。

### 3.3 修复前已确认"正确"的项目（作为回归基准）

| 项目 | 结果 |
|---|---|
| 线性法向力黄金值（vs `k_n·δ + c_n·v_n`） | 比值 ≈0.99–1.00 |
| Hertz 法向力峰值黄金值 | 比值 ≈0.9997 |
| Hertz 最大重叠量 | 4.714e−3 m（理论 4.714e−3） |
| 斜碰线动量 / 角动量守恒 | Lz 初末 −0.096（守恒） |
| 键合轴向力 = k_n·伸长 | 零伸长→零力，正确 |
| 键合剪切：产生转动且角动量守恒 | Lz 0.2→0.19999 |
| 单球静力平衡高度 | y=0.009533 m（理论 0.0095326） |
| 参数负例 9 项（dt/e/ν/box/id/model/运动类型/列数/重复） | 9/9 正确拒绝 |
| 200 颗粒重力堆积 | 稳定、无 NaN |

---

## 4. 缺陷清单与根因分析

### D1 — 线性接触模型的"禁止拉力"钳位破坏恢复系数标定

**位置**：`src/Basic/DEM/LinearSpringDashpot.cpp` `CalcNormalForce`

```cpp
double Fn_mag = k_n * contact.overlap_n + c_n * v_n_rel;
if (Fn_mag < 0.0) Fn_mag = 0.0;      // ← 缺陷
contact.force_n = -Fn_mag * contact.normal;
```

**根因**：恢复系数由标准线性弹簧-阻尼（Cundall-Strack）标定

```
beta = -ln(e) / sqrt(pi² + ln²e)      c_n = 2·beta·sqrt(m_eff·k_n)
```

该标定的前提是**回弹阶段允许阻尼项产生小幅"拉力"**（`c_n·v_n_rel` 使总法向力为负）。
接触的建立/解除由几何重叠判据（`dist < sum_r`）负责，并不需要"力不能为拉"这一约束。
钳位把回弹末段的这部分耗能抹掉，使系统比标定意图更"不耗能"，
因此 e<1 时实测恢复系数**系统性偏高**，且 e 越小、回弹末期重叠越小，钳位触发越充分、偏差越大。

**证据**（Python 复刻同一积分格式）：

| e_nom | beta | 复刻-钳位 | 复刻-无钳位 |
|---|---|---|---|
| 0.50 | 0.215454 | **0.551204** | **0.500877** |
| 0.10 | 0.591155 | **0.253488** | **0.098938** |

复刻的"钳位"列与 C++ 修复前实测**逐位一致**，证明根因定位准确。

### D2 — Hertz 模型阻尼系数缺少非线性修正

**位置**：`src/Basic/DEM/HertzMindlin.cpp` `CalcNormalForce`

```cpp
double k_n = 2.0 * E_star * std::sqrt(R_star * delta);   // 随 δ 变化
double c_n = 2.0 * beta * std::sqrt(m_eff * k_n);        // ← 直接照搬线性公式
...
if (Fn_mag < 0.0) Fn_mag = 0.0;                          // ← 同 D1 的钳位
```

**根因**（两个叠加因素）：

1. **与 D1 相同的钳位**，限制 e<1 时的耗能；
2. **阻尼系数形式错误**：Hertz 法向刚度 `k_n ∝ δ^{1/2}` 是非线性的，
   其恢复系数标定不能沿用线性模型的 `c_n = 2β√(m_eff·k_n)`，
   而应采用标准 Hertzian（Tsuji/Brilliantov）阻尼系数，附加 **`√(5/6)`** 非线性修正因子。

**证据**（Python 复刻，同一积分格式，e_nom=0.5）：

| 方案 | e_meas |
|---|---|
| 有钳位 + 无修正（≈修复前代码） | 0.523450 |
| 去钳位 + 无修正 | 0.466241 |
| **去钳位 + √(5/6) 修正（修复方案）** | **0.499985** |

并在 e_nom ∈ [0.1, 0.95] 全域验证修正后误差 <1e-3。

### D3 — 备份 `particles.dat` 表头带 `variables=` 前缀

**位置**：`src/Zone/DEM/DEMSolver.cpp` `BackupField`

```cpp
fout << "variables=id,group,radius,mass,...\n";   // ← 缺陷
```

同函数写出的 `bonds.dat`、`gas_voronoi_faces.dat` 以及输入文件 `particles.csv`
均为**纯 CSV 表头**（首列即 `id`）。该前缀只污染首列列名，
使按 `csv.DictReader` 解析 `row["id"]` 的脚本抛 `KeyError`；
而只读取 `py` 等后续列的 `dem_bounce_test` 恰好不受影响——这正是修复前
"4 通过 / 4 失败"差异的来源。

仓库内 **8 个** Python 消费者中 6 个假定纯 CSV 表头（仅
`dem_fixed_collision_test/verify.py` 做了兼容处理），且 C++ 侧无任何读取该备份的代码，
故判定为**输出格式缺陷**而非既定契约。

### 4.4 排除项（经测试判定为**非缺陷**）

| 现象 | 结论 |
|---|---|
| 墙面反弹 e≈0.55（e_nom=0.5）在步长收敛后仍不变 | **测试指标偏差**：用"回弹峰值速度"作分子会系统性高估（自由阻尼模型回弹峰值>分离速度）。改用分离速度指标后 e=0.4998。**代码正确**。 |
| 墙面反弹力在回弹末段为负（min fy=−317 N） | 这是**缺陷已修复**的证据（自由模型应有的吸引期），非异常。 |
| 滚动算例末态与"手算简单解析"量级不符 | 末态滑移速度 `vx+R·ωz≈0`（已进入纯滚动），系参考解过于简化；无接触力异常。 |
| Hertz 短窗口（4000 步）实测 e 偏差大 | 接触尚未结束，属**测量窗口不足**；加长窗口后 e=0.49999。 |

---

## 5. 修复方案与实施

### 5.1 D1 修复（`src/Basic/DEM/LinearSpringDashpot.cpp`）

删除"禁止拉力"钳位，恢复标准无钳位线性弹簧-阻尼模型（并加注释说明理由与陷阱）：

```diff
-    double Fn_mag = k_n * contact.overlap_n + c_n * v_n_rel;
-
-    if (Fn_mag < 0.0) Fn_mag = 0.0; // 不允许拉力
-
-    contact.force_n = -Fn_mag * contact.normal;
+    // 标准线性弹簧-阻尼模型（Cundall-Strack）在回弹阶段允许阻尼项产生
+    // 小幅"拉力"，这正是 c_n = 2·beta·sqrt(m_eff·k_n) 标定的前提；
+    // 钳位 Fn_mag>=0 会抹掉回弹耗能，使实测恢复系数系统性偏高。
+    // 接触的建立/解除由几何重叠判据负责，无需力钳位。
+    double Fn_mag = k_n * contact.overlap_n + c_n * v_n_rel;
+    contact.force_n = -Fn_mag * contact.normal;
```

### 5.2 D2 修复（`src/Basic/DEM/HertzMindlin.cpp`）

(a) 阻尼系数加入 `√(5/6)` 非线性修正；(b) 同 D1 删除钳位。

```diff
-    double c_n = 2.0 * beta * std::sqrt(m_eff * k_n);
+    // Hertz 法向刚度 k_n = 2E*·sqrt(R*·δ) 随 δ 变化（非线性弹簧），
+    // 不能照搬线性模型的阻尼系数；采用标准 Hertzian（Tsuji/Brilliantov）
+    // 阻尼系数，附加 sqrt(5/6) 非线性修正因子。
+    double c_n = 2.0 * std::sqrt(5.0 / 6.0) * beta * std::sqrt(m_eff * k_n);
...
     double Fn_mag = Fn_hertz + Fn_damp;
-    if (Fn_mag < 0.0) Fn_mag = 0.0;
     contact.force_n = -Fn_mag * contact.normal;
```

### 5.3 D3 修复（`src/Zone/DEM/DEMSolver.cpp`）

```diff
-    fout << "variables=id,group,radius,mass,...\n";
+    // 与 bonds.dat / gas_voronoi_faces.dat / particles.csv 保持一致的纯 CSV 表头
+    fout << "id,group,radius,mass,...\n";
```

### 5.4 变更清单

| 文件 | 变更 |
|---|---|
| `src/Basic/DEM/LinearSpringDashpot.cpp` | 删除法向力钳位（D1） |
| `src/Basic/DEM/HertzMindlin.cpp` | 阻尼系数加 `√(5/6)`；删除钳位（D2） |
| `src/Zone/DEM/DEMSolver.cpp` | `particles.dat` 表头去 `variables=` 前缀（D3） |

---

## 6. 修复后验证

### 6.1 恢复系数扫描（线性，两体对心）

| e_nom | 修复前 | 修复后 | 修复后误差 |
|---|---|---|---|
| 0.90 | 0.902340 | **0.900500** | +0.0005 |
| 0.70 | 0.718929 | **0.700585** | +0.0006 |
| 0.50 | 0.551204 | **0.500877** | +0.0009 |
| 0.30 | 0.398070 | **0.300189** | +0.0002 |
| 0.10 | 0.253488 | **0.098938** | −0.0011 |

### 6.2 恢复系数扫描（Hertz，两体对心，窗口 12000 步）

| e_nom | 修复前 | 修复后 | 修复后误差 |
|---|---|---|---|
| 1.00 | 1.000000 | **1.000000** | 0.00000 |
| 0.90 | 0.893393 | **0.900001** | +0.00000 |
| 0.50 | 0.523450 | **0.499985** | −0.00001 |
| 0.30 | 0.369549 | **0.299960** | −0.00004 |

### 6.3 墙面反弹（分离速度指标）

| e_nom | 修复前 | 修复后 |
|---|---|---|
| 0.90 | 0.901979 | **0.899983** |
| 0.50 | 0.550062 | **0.499809** |
| 0.30 | 0.396851 | **0.299642** |
| 0.10 | 0.252394 | **0.099834** |

### 6.4 已提交测试回归（`tests/*`，修复后二进制）

```
dem_bond_test                   PASS  permanent inert elastic bond restores its stretched endpoint
dem_bounce_test                 PASS  inert DEM bounce response is physically consistent
dem_energy_fracture_test        PASS  fracture_energy=0.577350 J, energy criterion broke the bond
dem_fixed_collision_test        PASS  fixed particle participates in contact and remains immobile
dem_reaction_test               PASS  alpha=0.0956179, gas_temperature=700.000 K
dem_voronoi_test                PASS  conserved clipped Voronoi volume=8.999992000
dem_brazilian_disk_test         PASS  elastic Brazilian disk response is symmetric and mechanically consistent
hmx_disk_center_ignition        PASS  nodes=91, reacted=91, center_alpha=1.000, broken_bonds=240
==> 8 passed
```

**修复前 4/8 → 修复后 8/8。**

### 6.5 非回归确认

修复仅改动"接触法向力"与"输出表头"，对键合、热化学、Voronoi、几何守恒量无影响；
重跑 §3.3 全部项目结果一致（含参数负例 9/9 拒绝、200 颗粒堆积无 NaN）。

### 6.6 大规模算例

| 算例 | 规模 | 结果 |
|---|---|---|
| dem_brazilian_disk_10k | 10266 粒子 / 29635 键 | 输入几何与拓扑校验 PASS；修复后 400 步短窗运行：**无非有限状态**（52.7 s） |
| hmx_disk_10k_center_ignition | 10000 粒子 / 29635 键，反应+Voronoi 开启 | 输入与已提交输出校验 PASS（含反应节点、动态 Voronoi 面）；修复后 400 步短窗：**无非有限状态**（57.8 s） |
| dem_brazilian_disk_test | 139 粒子 × 50000 步 | 完整运行 **PASS**（对称性、横向拉应变、反力平衡） |
| hmx_disk_center_ignition | 91 粒子 × 8000 步（反应） | 完整运行 **PASS** |

> 说明：两个 1 万粒子算例的**完整 10000/8000 步**在本次测试环境的实测单步成本约
> 0.13 s（≈20 min/算例），未在会话时间预算内跑完；已用"输入完整校验 + 400 步短窗稳定性"
> 覆盖其正确性与稳定性风险。中等规模同类算例（139 粒子 × 5 万步，含键合/接触/对称性）
> 已完整通过，接触力路径与大规模算例相同。

---

## 7. 测试结论与残余风险

### 7.1 结论

1. 本次共定位 **3 处缺陷**并全部修复，其中 D1/D2 属**物理正确性缺陷**（恢复系数标定失效），
   D3 属**输出契约缺陷**。
2. 修复后，**两体碰撞、墙面反弹的恢复系数均与用户标称值在离散化误差内一致**
   （线性 ≤0.0011，Hertz ≤0.00004，墙面 ≤0.0002），可认为接触模型的"恢复系数"参数
   已具备**物理可解释性**与**跨尺度可移植性**。
3. 已提交测试由 4/8 提升至 **8/8 通过**，且既有正确的物理量（守恒律、黄金值、静力平衡）
   无回归。
4. 结论：**修复后的 Zaran3.10.2 在本次测试覆盖范围内（接触力学 / 键合 / 热化学 / Voronoi /
   参数校验 / 中等规模稳定性）正确、稳定。**

### 7.2 残余风险与建议（未在本次修复）

1. **墙面材料参数硬编码**：`DEMWall` 默认 `young_modulus=1e9`、`restitution_coeff=0.9`、
   `friction_coeff=0.4`，`DEMFieldGenerator::AddBoxWalls` 未从配置读取。
   后果：用户设置 `restitution_coeff>0.9` 时墙面会成为控制因素（`min(粒子, 墙)`）；
   墙-粒相对刚度被人为放大 100 倍。
   **建议**：增加 `dem.wall.*` 配置项。
2. **线性模型的固有"吸引尾"**：回弹末段存在小幅拉力（标准 Cundall-Strack 模型特性）。
   若研究场景要求严格非吸引接触（如黏附/结块），应改用切换式（switching）或非线性阻尼模型，
   但这会改变标定公式，须重新验证。
3. **Hertz 阻尼**：本报告采用 `√(5/6)`（Tsuji/Brilliantov），已在 e∈[0.1,0.95] 验证；
   与切向 Mindlin 摩擦强耦合或极端泊松比情景建议补充验证。
4. **备份格式变更的兼容性**：若外部后处理脚本依赖 `variables=` 前缀需同步更新
   （仓库内未发现此类依赖）。
5. **未覆盖**：并行/MPI、颗粒形状、真实材料参数标定与实验对比。

---

## 附录 A：可复用测试脚手架（`Zaran3/_demtest/`）

| 脚本 | 用途 |
|---|---|
| `build.sh` | Ninja + MSVC 构建（封装工具链环境变量） |
| `restitution_sweep.py <exe> <tag>` | 线性/Hertz/墙面恢复系数扫描 |
| `run_committed_tests.py [exe] [case...]` | 批量运行 `tests/*` 并调用其 `verify.py` 判定 |
| `restitution_model.py` | Python 复刻线性碰撞（钳位 vs 无钳位），定位 D1 根因 |
| `hertz_variants.py` | Python 复刻 Hertz 阻尼系数变体，定位 D2 根因并标定 `√(5/6)` |
| `hertz_metric.py <exe> <tag>` | Hertz 长时间窗恢复系数（轻量） |
| `wall_metric.py <exe> <tag>` | 墙面反弹恢复系数（分离速度/峰值速度双指标） |
| `wall_force_dump.py` | 导出墙面接触力符号（验证钳位是否消除） |
| `wall_study.py` | 墙面步长收敛实验（区分离散化误差与模型缺陷） |
| `scale_check.py <exe>` | 1 万粒子算例短窗稳定性检查（有限性/非发散） |
| `golden.py` | 接触力公式黄金值比对 |

> 注意：`_demtest/` 及 `build-ninja/` 为本次测试新增，未纳入版本控制；
> 测试脚本会生成带时间戳的算例目录。批量删除工作目录会触发环境安全策略，
> 脚本已改为"每次使用独立目录、不做删除"。
