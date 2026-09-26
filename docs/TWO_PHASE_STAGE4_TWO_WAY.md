# DEM-CFD 双向耦合 + 方形域二维球形黎曼问题（阶段 4）—— 实现与验证

> **阶段目标**：把阶段 3 的**单向**耦合升级为**双向**耦合：
> ① 颗粒对流场有反作用（动量 + 做功）；② 颗粒自身在水动力下运动。
> 同时**换掉基础算例**：不再用"x 方向 1 m、y/z 只有 2–3 格"的长条薄板，
> 改用**方形计算域**，并以**二维球形黎曼问题**（内部球形静止高压、外部静止低压）
> 作为主要考察对象。
>
> 上游：`docs/TWO_PHASE_STAGE3_DEM_CFD.md`（单向耦合）、
> `docs/TWO_PHASE_EULER_STAGE1.md`（ε 加权 Euler）、
> `docs/TWO_PHASE_STAGE2_MAPPING.md`（粒子 → 网格映射）。

---

## 1. 为什么换掉旧激波管算例

| 问题 | 旧 `tests/demcfd_shock_tube` | 新设计 |
|---|---|---|
| 网格各向异性 | x 方向 1 m、`dy=dz=0.02` 而 `dx=0.0025` | **dx = dy = dz**，三个方向的 CFL 贡献相同 |
| y/z 分辨率 | 只有 1 个单元厚且与 x 不成比例 | 方形 x-y 截面（81×81），z 方向 1 单元厚（二维化） |
| 粒子可摆布性 | y/z 里只能占 1 格 | 能在 x-y 平面里按对称群布置 |
| 流场结构 | 1 维波系 | 2 维球（柱）形黎曼：外传激波 + 内传稀疏 + 接触间断 |

"二维问题放进三维求解器"仍然用薄板：`nz = 3`、`dz = dx`，
初值在 z 方向**严格不变**（`init.explosion.dim = 2` 只用 (x,y) 距离 ⇒ 柱形），
六面 wall 使解保持 z 无关。

---

## 2. 反作用力：加什么、不加什么（最容易双计的地方）

气相的离散方程（与阶段 1 一致，`cons = ε·cons_intrinsic`）：

```
∂(ε ρ u)/∂t + ∇·(ε ρ u u + ε p I) = p ∇ε + S
```

连续形式的分工（Model A）：

| 力 | 在气相方程里的对应项 | 要不要写进 S |
|---|---|---|
| 颗粒受的**压力梯度力** `−V_p∇p` | 已经是 `−ε_g∇p` 的**对偶项**：颗粒侧的 `Σ_p(−V_p∇p)/V = −ε_s∇p`，与气相的 `−ε_g∇p` 相加恰好是混合物的 `−∇p` | ❌ **不加**（再加一次就是双计） |
| 颗粒受的**曳力** | 气相方程里**没有**任何对应项 | ✅ 加 `−Σ_p F_drag,p / V_cell` |

把 `p∇ε` 这一项拿掉再看更清楚：`−∇·(εpI) + p∇ε ≡ −ε∇p`。
所以 `S` 里**只放曳力**。开关 `coupling.react_pressure_gradient_force` 有 `auto/on/off` 三档：
`auto` 会按 `two_phase.porosity_gradient_force` 自动选择 ——
取 `ε∇p` 形式（=1）时是 `off`，取完全守恒的 `∇(εp)` 形式（=0）时才是 `on`（那时浮力必须自己补）。
两者配错会告警。

**相间功**：相间力对气相要做功，能量方程里应当有 `S·u`：

```
∂(ε E)/∂t + ∇·(ε (E+p) u) = S·u
```

少了这一项，气相就是"白拿动量"（拿到动量却不多出一份动能）。
默认打开（`coupling.react_interphase_work = true`），且因为它在整个 RK 阶段是常量
（`S` 滞后、`u` 取步初值），`Σb_i = 1` 使一步的总功**精确**等于 `Σ_cells S·u·V·Δt` ——
于是能量账本可以逐位核算。

---

## 3. 实现

### 3.1 求解器：一个源项钩子

`EulerTwoPhaseStructUniform`：

```cpp
void SetInterphaseMomentumSource(const double* src3);  // 3*N 节点数组，nullptr = 关
void SetInterphaseWork(bool enabled);                  // 是否同时加 S·u 进能量方程
void CalcSourceResidual() override;                    // 基类是空实现
```

`CalcSourceResidual()` 把 `src[3*idx+dir]` 加到第 `dir` 个动量残差，
并（若开启相间功）把 `src[3*idx+dir]*u_dir[idx]` 加到能量残差。
源项在整个 RK 阶段是常量，所以**一步的总贡献 = S·Δt**，可以被精确核算。

### 3.2 耦合器：散列 + 推进颗粒

```cpp
void DEMCFDCoupler::ScatterReaction(results, include_pressure_gradient_force,
                                    src_dimensionless, ReactionSummary&);
int  DEMCFDCoupler::AdvanceParticles(dt, MotionOptions, results);
```

`ScatterReaction` 用**与体积分配/插值完全相同**的权重把 `−F_p` 摊回网格节点，
因此 `Σ_cells src·V_cell = −Σ_p F_p` 是构造保证的（Σw = 1 时逐位成立）。
转换关系：`S* = F* / V*_cell`，`F* = F_SI/F_ref`，`V*_cell = V_cell/L_ref³`
⇒ `factor = L_ref³/(F_ref·V_cell)`，**只能对 `F_ref` 除一次**。

`AdvanceParticles` 支持两种积分器：

| 积分器 | 更新式 | 特点 |
|---|---|---|
| `exponential`（默认） | `u ← u_∞ + (u − u_∞)·exp(−k·Δt)`，`k = V_p(β/ε_s)/m`，`u_∞ = u_g + a_ext/k` | 对**冻结线性曳力是解析解**，无条件稳定（Stokes 下实测逐步偏差 0.000e+00） |
| `rk2` | 显式中点法（对冻结气态） | 需要子步（自动按 `k·Δt ≤ 0.25` 定） |

### 3.3 驱动器：交错（lagged）耦合

```
CouplingPreStep(iter)                       ← 每个时间步**之前**
   ① 颗粒动过 ⇒ 用当前位置重建 ε_g 并注入（ghost 随后由 BoundaryCondition 补）
   ② 用**步初**气态算受力 → m_results_pre
   ③ ScatterReaction → 交给求解器（这一步的源项）

   … PreSolve / SolveOneStep(含 S 的 RK) / PostSolve …

CouplingPostStep(iter)                      ← 每个时间步**之后**
   ① 用**同一份**力 m_results_pre、按气体真实走过的 Δt 推进颗粒
   ② 累计动量/能量账本
   ③ 在**步末**气态上重算受力 → m_results（供输出与下游校验）
   ④ 输出 CSV / Tecplot 点云
```

②③ 分开保存是刻意的：散射进流体的力与推进颗粒的力必须是**同一份**（能量/动量账本才闭合），
而输出必须对应当前帧的气态（这样 `particle_force.csv` 的每一行都能与 `result/<iter>.dat` 逐值对照，
阶段 3 的 V4/V5/V6/V8 才继续成立）。

### 3.4 参数（`[coupling]` 段新增）

```toml
particle_motion   = "fixed" | "prescribed" | "dynamic"
particle_integrator = "exponential" | "rk2"
particle_substeps   = 0            # 0 = 自动
particle_max_step_parameter = 0.25 # rk2 自动子步的上限 k·Δt
body_acceleration_x/y/z = 0        # 体积力/质量 (m/s²)，如重力
react_on_fluid      = true|false   # 反作用力是否进流体
react_interphase_work = true       # 相间功 S·u 是否进能量方程
react_pressure_gradient_force = "auto" | "on" | "off"
remap_volume_fraction = true       # 颗粒在动 ⇒ 每步重建 ε
```

### 3.5 粒子输出：ZONE 行

按用户要求，**FEPOINT 点云的 ZONE 行默认只写 `T=` 与 `SOLUTIONTIME=`**：

```
ZONE T="Particles", SOLUTIONTIME=0.0
```

之前写的 `, N=8, E=0, ZONETYPE=FEPOINT, DATAPACKING=POINT` 在某些 Tecplot 里不被识别。
`[output] particle_zone_keywords = none | modern | classic` 三档可切
（`modern` = `ZONETYPE=FEPOINT` + `DATAPACKING=POINT`；`classic` = `F=FEPOINT`）。

⚠️ **键合的 `FELINESEG` zone 必须保留** `ZONETYPE=FELINESEG, DATAPACKING=BLOCK, VARLOCATION=…`
—— 线元 zone 带着单元中心量，去掉这些关键字就不可能被正确读取。
`benchmarks/tecplot_ascii/check_tecplot_dat.py` 已同步：遇到"既没有 ZONETYPE 也没有 I/J/K"
的 zone 就按"一行一个点"的 FEPOINT 处理（并用 backup CSV 的行数交叉核对）。

---

## 4. 算例 A：颗粒在静止流场中运动（`tests/demcfd_still_fluid`）

**目的**：把"颗粒运动对流体产生的反作用"单独拎出来验证。

* **方形计算域**：x,y ∈ [0,1]（41×41 节点，`dx=dy=dz=0.025`），z ∈ [0,0.05]。
* **静止流场**：全场 ρ=1、p=1、u=0，六面 wall。初值用同一个状态填 `[init.riemann]` 左右两段
  ⇒ 没有间断，什么都不会发生。
* **颗粒**：9 个（3×3，间距 0.05），r=4 mm（Δ/d_p = 3.125），ρ_p=1000，**匀速 0.05 无量纲**
  （`particle_motion = "prescribed"`，SI 速度 16.56 m/s）；颗粒会平移 ⇒ 每步重建 ε。
* 跑到 t = 0.2（57 步），`write_interval = 1` 便于拿相邻两帧做单步检验。

```bash
python tests/demcfd_still_fluid/generate_case.py
bin/Release/Zaran3.10.2.exe tests/demcfd_still_fluid
python tests/demcfd_still_fluid/verify.py
```

### 4.1 实测结果（15 项断言全过）

| 检验 | 判据 | 实测 |
|---|---|---|
| **W1 牛顿第三定律**：源冲量 == −Σ_p F* | 相对 L1 尺度 ≤ 1e-12 | x/z **0.000e+00**，y 2.4e-22；C++ 侧逐步最大散列残差 **2.04e-16** |
| **W2 单步逐格注入** | 相对残差 ∝ Δt | CFL=0.5：8.41e-2；**CFL/4：2.10e-2 ⇒ 降了正好 4.00 倍** |
| W2 第一步能量 | 无虚假相间功 | 相对 9.3e-8（O(Δt²) 的通量效应） |
| **W3 全局动量账本**（第一步） | ≤ 1e-10 | **2.01e-14** |
| W3 全程 | 仅报告 | P* = 5.203e-8 vs 冲量 5.211e-8（差 0.16%，壁面压强回推） |
| **W4 退化**：关反作用 | 流场逐位不变 | **58 帧 × 5043 节点，max\|Δ\| = 0.000e+00；max\|u\| ≡ 0** |
| **W5 颗粒积分器** | 逐步递推 ≤ 1e-12 | **0.000e+00**（指数积分器对冻结线性曳力是解析解） |

**W2 是这一阶段最有信息量的一条**：均匀静止态的**所有界面通量都相同**（含壁面），
所以第一步的动量变化**逐格**必须等于 `S·Δt`（`S` 由 Python 用同一套 SubCell 权重从
`particle_force.csv` 独立复算）。实测残差不是 0 而是 `S·Δt` 的 8.4% ——
它不是实现误差，而是**通量散度**：反应力让第一步末的流场不再均匀，
第 2 个 RK 子步的界面通量差就不再为零，量级 `∝ c·Δt/Δx`。
于是判据写成**收敛性**：Δt 减 4 倍 ⇒ 相对残差降到 1/4（实测 4.00 倍）。

> 一个值得记住的物理结论：**单颗粒在静止气体里运动时，反作用是"声学筛除"的**。
> 本算例里气体只被加速到 ~1e-6 的无量纲速度（马赫 1e-6），
> 因为低马赫下压强偶极子几乎立刻把施加的力顶回去。
> 想看到可观的反作用，要么把颗粒铺满很大一部分体积（堆积带），要么让气体能逃逸（开口边界）。
> 这也解释了为什么"反作用验证"必须走**精确性**路线（牛顿第三定律、逐格注入、
> 退化检验），而不是靠"看得见的变化"。

---

## 5. 算例 B：方形域 + 二维球形黎曼问题（`tests/demcfd_ball_riemann`）

**初始条件**：半径 `R0 = 0.15` 的圆内 (ρ=1, p=1, 静止)，圆外 (ρ=0.125, p=0.1, 静止)。
用 `field_type = "Explosion"` + **新增的** `init.explosion.dim = 2` 实现 ——
`dim=2` 时距离只用 (x,y)，薄板里初值严格与 z 无关（即柱形黎曼问题）。

```bash
python tests/demcfd_ball_riemann/generate_case.py --single-phase --out tests/demcfd_ball_riemann/ref_ball
python tests/demcfd_ball_riemann/generate_case.py
python tests/demcfd_ball_riemann/verify.py       # 29 项断言
```

* 网格 81×81×3，`dx=dy=dz=1/80`；六面 wall；t_end = 0.15（激波传到内切圆边界约需 0.1998）。
* **粒子按 D4 群对称布置**（全部精确落在网格节点上）：
  圆心 1 个 + 坐标轴 r=0.2 处 4 个 + 对角线 r=0.3005 处 4 个；
  r=2 mm（Δx/d_p = 3.125），ρ_p = 20 kg/m³（取得轻，冲击波才能明显推动它）。

### 5.1 为什么这个算例的主判据是"对称性"而不是"动量守恒"

D4 对称布置（粒子分组互为镜像）⇒ 流体给颗粒的合力**矢量恒为 0**
⇒ 反作用注入的总动量为 0（实测 ΣF* ≈ 1e-22，即机器零）。
所以"动量账本"在这个算例里是**平凡**的，动量与能量的账本检验留给算例 A。
反过来，这份对称性给了一组很硬的判据。

### 5.2 实测结果

| 检验 | 判据 | 实测 |
|---|---|---|
| **B1 初值自检** | 内外状态、静止 | 内 1263 节点 / 外 18216 节点全部合规，max\|u\| = 0 |
| **B2 D4 对称性**（流场不变量 ρ, p, \|u\|, u_r） | ≤ 1e-9·尺度 | 参考 **≤ 2.44e-15**；耦合 **≤ 1.89e-15**（x↔y 恰好 0.000e+00） |
| **B3 守恒** | wall 下 ∭ρ、∭E | 单相 **1.52e-13 / 1.73e-13**；耦合 ∭ερ **1.51e-13** |
| B3 能量账本（耦合） | Δ∭εE vs −Σ∫F_drag·u_g dt | Δ∭εE = −3.029e-7 vs 功 = −3.050e-7，**相对偏差 0.683%** |
| **B4 体积配分（端到端）** | ≤ 1e-9 | **3.51e-16** |
| **B5 颗粒群对称性** | 互为镜像者受力/速度逐位一致 | 轴上 4 个散布 **8.3e-16**；对角 4 个 **1.9e-15** |
| B5 圆心颗粒 | 受力 ≡ 0、一动不动 | \|F\| = 3.4e-20 N（相对 4.6e-19）；\|u_p\| = 3.0e-15 m/s（相对 8.4e-17） |
| **B6 冲击波扫过顺序** | 由内到外 | 内圈 r=0.200 起始 t=0.01261 < 外圈 r=0.3005 起始 t=0.07526 |
| **B7 颗粒动力学自洽** | ≤ 5% | r=0.200：36.165 vs Σ(F/m)Δt=36.407（**0.67%**）；r=0.3005：10.068 vs 10.212（**1.43%**） |
| **B8 ε≡1 + 无反作用 == 单相参考** | 逐位 | **0.000e+00**（16 帧 × 19683 点 × 8 列） |
| B8 耦合算例 ≠ 参考 | 必须不同 | 末帧 max\|Δρ\| = **6.25e-3** |
| **B9 牛顿第三定律** | 对称布置下 ΣF* ≡ 0 | 源冲量 1.176e-22 / 颗粒冲量 1.157e-22（尺度 4.41e-7）；散列残差 **1.51e-16** |

颗粒响应：轴上 4 个被激波加速到 **36.17 m/s**（≈11% 声速），
对角 4 个到 **10.07 m/s**（r 更大、激波已被几何发散削弱），圆心颗粒一动不动。
参考算例末帧激波前沿 R = 0.3750（R0=0.15）⇒ 平均速度 1.50（平面 Sod 是 1.7522，
柱形发散确实把它拖慢了）。

---

## 6. 过程中踩到的四个新坑

### 6.1 无量纲 `dt` 被直接用在 SI 的颗粒上

`AdvanceParticles` 拿到的是求解器的无量纲步长（`t* = t·a_ref/L_ref`），
而颗粒状态（位置 m、速度 m/s）是 SI。
直接相乘 ⇒ 位移被放大 `a_ref/L_ref ≈ 331` 倍 ——
实测颗粒在 0.15 的无量纲时间里"跑"了 0.34 m，**直接飞出计算域**，
连带把体积配分残差顶到 0.444（4/9 的粒子出了域，权重为空）。
修正：`dt_si = dt·L_ref/a_ref`，并在 `AdvanceParticles` 末尾检查粒子是否跑出控制体，
只告警一次（不点名的话很容易把"初值/时长选错"误读成"耦合 bug"）。

### 6.2 散列源项的单位换算只能除一次 `F_ref`

```cpp
// 错：ForceToDimensionless 已经除过 F_ref 了
src -= w * ForceToDimensionless(f_si) * (L³ / (F_ref * V_cell));
// 对：直接乘原始 SI 力
src -= w * f_si * (L³ / (F_ref * V_cell));
```
错版本的等效系数是 `1/F_ref²`，反作用被缩小 **5 个量级**
（实测源冲量 4.8e-13 vs 颗粒冲量 5.2e-08）。
**判据**：`Σ_cells src·V_cell` 必须等于 `−Σ_p F*` —— 修复后实测 2.04e-16。

### 6.3 `iteration.dt` 与实际推进的时间不一致

`NSFieldSimulation::CalcTimeStep()` 先把**未削减**的 `dt` 写进 `GlobalData`，
之后才为了"正好落在 `end_time`"把它削减并**单独**更新 `current_time`。
于是最后一步里 `iteration.dt`（3.52e-3）与实际推进的时间（2.80e-3）不一致，
颗粒被多推了 26% 的一步，轨迹末点偏差 0.36%。
修正：颗粒推进改用 **`current_time` 的增量** —— 那才是气体真正走过的物理时间。
（这也让"用 CSV 的 time 列重构积分"变成严格可比的检验。）

### 6.4 静止流场会被"收敛判据"误杀；对称布置会让残差分母为零

* 均匀静止态的残差恒为 0 ⇒ `log10(0) = -inf < -min_residual_order` **立刻满足收敛条件**
  ⇒ 算例只跑 1 步就"算完了"（实测 `demcfd_still_fluid_settling` 只有 2 行输出）。
  静止类算例要把 `residual_interval` 取得极大，让 `CalcResidual` 根本不跑。
* 散列残差的分母**不能取 ΣF\***：D4 对称布置下它恒为 0，
  残差会虚报 **277%**。要取 `Σ_p Σ_d |F*_d|`（L1 尺度）。修好后 ≤ 1e-10。

---

## 7. 局限与下一步

**本阶段的已知局限**

1. **单向耦合仍然是"点源"**：曳力用插值到粒子中心的气态；Δx/d_p ≳ 3 时成立。
   分辨率不足时插值气态代表不了真实绕流。
2. **颗粒自身仍被算进自己的 ε_g**（`ε_g` 里含自己贡献的 α_s），曳力偏大 ~0.07%；
   `particle_force.csv` 里显式输出 `ε_g` 便于识别。
3. **时间耦合是一阶（lagged）**：气相加颗粒都用步初的力。要二阶需要
   `CouplingPreStep` 里预测半步 + `CouplingPostStep` 里校正。
4. **无湍流、无对流换热、无相变/燃烧**；粘性项未启用。
5. **单位制仍是"求解器无量纲 + 颗粒 SI"**，两个方向都要手工换算一次
   （本阶段两个 bug 都出在这里）。彻底的做法是让求解器也支持 SI 输入。
6. **球形黎曼算例缺一个定量的流场参考解**：目前只做了对称性/守恒/拓扑（内外顺序）。
   要做定量对照，最直接的是写一个**一维柱对称 HLLC 参考解**（同一初值、细径向网格），
   把 2-D 的径向剖面与它比较 —— 这是下一步最容易见效的一件事。
7. **二维化靠薄板 + wall**：z 方向 1 个单元；若要真正的 2-D，需要 z 方向的
   对称/周期边界或一个二维求解器。

**下一步的自然阶梯**

| 级别 | 内容 |
|---|---|
| 1 | 一维柱对称参考解 → 2-D 径向剖面的定量对照（补 B2/B3 的定量缺口） |
| 2 | 二阶时间耦合（半步预测 + 校正），看能量账本残差能否再降一个量级 |
| 3 | 单球定常绕流的 `C_D–Re`、固定床压降 vs Ergun、最小流化速度（真正的双向耦合验收） |
| 4 | 颗粒床被激波穿透（算例 B 的堆积带 + 双向），激波衰减/反射系数 |

---

## 8. 文件清单

```
新增
  tests/demcfd_still_fluid/{generate_case.py, verify.py, zaran.toml, mesh.*, particles.csv}
  tests/demcfd_ball_riemann/{generate_case.py, verify.py, zaran.toml, mesh.*, particles.csv}
  tests/demcfd_ball_riemann/ref_ball/（单相 Euler 参考解，由生成器产出）
  docs/TWO_PHASE_STAGE4_TWO_WAY.md（本文）

改动
  inc/Zone/Solver/EulerTwoPhaseStructUniform.h   + SetInterphaseMomentumSource / SetInterphaseWork
  src/Zone/Solver/EulerTwoPhaseStructUniform.cpp + CalcSourceResidual（动量 + 相间功）
  inc/Basic/TwoPhase/DEMCFDCoupler.h             + MotionMode/Integrator/MotionOptions/
                                                   ReactionSummary/ScatterReaction/AdvanceParticles
  src/Basic/TwoPhase/DEMCFDCoupler.cpp           同上（含单位换算与出域告警）
  inc/Main/DEMCFDSimulation.h  src/Main/DEMCFDSimulation.cpp
                                                 + [coupling] 双向参数、动量/能量账本、
                                                   CSV 由 36 列扩到 47 列（末尾追加，旧列名不变）
  src/Zone/Solver/EulerSolverStructUniform.cpp   + init.explosion 的
                                                   dim / inner|outer_density / pressure / velocity
  src/Basic/PostProcess/Visual.cpp  inc/.../Visual.h
                                                 + output.particle_zone_keywords（none|modern|classic）
  benchmarks/tecplot_ascii/check_tecplot_dat.py  + 无 ZONETYPE 时推断为 FEPOINT
```
