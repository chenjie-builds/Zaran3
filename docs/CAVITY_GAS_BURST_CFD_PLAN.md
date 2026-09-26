# 把 HMX 高压气腔撑碎算例从"纯 DEM"改为"DEM-CFD 耦合"：工作清单

> 目标：现在的 `tests/hmx_cavity_burst*` 走 `solver = "DEM"`，
> 气体是 `CavityGasModel` 提供的**0 维集中参数压力**（整个空腔一个 p，
> 作用在空腔边界粒子环上）。希望换成真正的 CFD 场：
> 压力有空间分布、气体能沿着裂纹流动、能把压力作用到**裂纹面**上。
>
> 本文只做盘点与规划，**不改代码**。

---

## 一、一句话结论

**这不是"加一个 CFD 模块"，而是把两条目前互斥的链路合并。**

| 链路 | 入口 | 有什么 | 缺什么 |
|---|---|---|---|
| A | `solver = "DEM"` | `DEMSolver`：键 / 接触 / 断裂 / Weibull / 热化学 / JWL / 0 维气腔 | **没有 CFD 场**（无网格、无流场） |
| B | `solver = "EulerDEM"` | 两相 Euler（ε 加权、可压缩、HLLC）+ `DEMCFDCoupler`（双向耦合，阶段 3/4 已完成） | 颗粒是**点粒子 + 曳力 ODE**：无键、无接触、无断裂 |

证据：
- 分派在 `src/Main/Application.cpp:187-245`（按 `task.solver` 二选一）。
- 链路 B 的颗粒推进只有平动积分：`DEMCFDCoupler::AdvanceParticles()`（`src/Basic/TwoPhase/DEMCFDCoupler.cpp:360-458`），
  **不实例化 `DEMSolver`**。
- 链路 A 完全不碰 `EulerSolverStructUniform` / `DataManagerNSTwoPhase`。

⇒ 工作量的重心在**时间循环与数据所有权的合并**，外加两块物理/数值的新东西：
**自由表面压力**与**高压场景下的数值稳健性**。

---

## 二、现有能力的准确边界（先把"能复用"和"缺"分清）

### 2.1 CFD 侧（可复用）

| 项 | 现状 | 位置 |
|---|---|---|
| 方程 | 5 方程可压缩无粘 Euler，含能量；MUSCL + HLLC + SSP-RK2 | `inc/Zone/Solver/EulerSolverStructUniform.h:48` |
| 两相 | `cons = ε·cons_intrinsic`，ε 加权的单气相（不含颗粒相动量方程） | `EulerTwoPhaseStructUniform.h:15-22` |
| 网格 | **只支持均匀笛卡尔结构网格**（`CalcCoordTransCoef` 强制 `max_dev < 1e-6`） | `EulerSolverStructUniform.cpp:433-486` |
| 维度 | 2D / 3D；1D 用薄板模拟（`nk=1 → 2D`） | `GridGeneratorStructGridgen.cpp:148-155` |
| 边界 | **只有 `inlet` / `outlet` / `wall`**，且为常值、**不随时间变化** | `EulerSolverStructUniform.cpp:1051-1131` |
| 初始条件 | `init.explosion.*`，`dim=2` 为二维柱形、`dim=3` 为球 | `EulerSolverStructUniform.cpp:232-309` |
| EOS | **只有理想气体 `PerfectGas`**；`RealGas` 枚举存在但无实现 | `Gas.h:16-20`、`FlowSolverPara.cpp:38` |
| 粘性 | CFD 方程无粘；Sutherland 只被耦合器用来算 Re | `EulerSolverStructUniform.h:86`、`DEMCFDSimulation.cpp:568-594` |

### 2.2 耦合侧（可复用）

| 项 | 现状 |
|---|---|
| 双向耦合 | 曳力 + 压力梯度力 → 颗粒；反作用力散列 + 相间功 → 流体（阶段 4） |
| 曳力模型 | 5 种（stokes / schiller_naumann / wen_yu / ergun / gidaspow，默认 gidaspow） |
| 粒子→网格 | `ParticleGridMapper`：`subcell` / `exact`，**构造保证体积守恒**（实测残差 3.5e-16） |
| 颗粒积分器 | `exponential`（对冻结线性曳力**无条件稳定**）/ `rk2` |
| ε 注入 | `VolumeFractionType::External` + `coupling.remap_volume_fraction` |
| 输出 | 流场 Tecplot、`particle_force.csv`（36 列）、粒子点云、`coupling_report.csv` |

### 2.3 明确不存在的（必须新写）

1. **完整 DEM 与 CFD 的组合驱动**（两条链路没有共同的 `Simulation`）。
2. **CFD 与 DEM 的时间步协调**：链路 B 用 CFD 的 dt，链路 A 用 `dem.time_step`，两者无 sub-cycling。
3. **自由表面压力（面力）**：耦合器只给 `−V_p∇p` 体积力，没有"空腔壁/裂纹面上的压力面积分"。
4. **非理想 EOS**：CFD 侧没有 JWL，也没有 `p(ρ,e)` 形式的 EOS 接口。
5. **随时间变化的边界**（活塞、随时间变化的背压）。
6. **ε_g→0 的刚性问题**：`two_phase.volume_fraction_min`（默认 1e-3）**只用来告警**，求解器内部不做任何钳位；
   耦合器侧 `coupling.eps_min_clip`（默认 1e-3）是唯一的地板。Ergun 支的 `β/ε_s ∝ ε_s/ε_g²`
   在 `DragForceModel.cpp:113-115` 只把 ε 钳到 1e-6，**限制但不消除**爆炸。

---

## 三、工作块 1：架构合并（工作量最大）

### 1.1 新的组合驱动
- **内容**：新建一个 `Simulation`（暂名 `DEMCFCTwoWaySimulation`），同时持有
  `DEMField`（粒子 + 键 + 接触模型）与 `NSFieldStructTwoPhase`（两相 Euler 场），
  并复用 `FieldSimulation` 的 `PrepareCoupling / CouplingPreStep / CouplingPostStep` 三段式。
- **现状**：`DEMFieldSimulation` 与 `DEMCFDSimulation` 各写一套；`DEMCFDSimulation` 的
  「颗粒」是自造的 `DEMCFDCoupler::Particle` 结构，不是 `DEMParticle`。
- **要点**：把 `DEMCFDCoupler::Particle` 换成/桥接 `DEMParticle`（含 `mass/radius/young_modulus/
  friction_coeff/inertia/omega/...`），否则键、接触、转动都无从谈起。
- **验收**：DEM 场与 CFD 场同时初始化、同时逐帧输出、时间严格一致。

### 1.2 `DEMSolver` 接受外部时间步
- **现状**：`DEMSolver::Solve()` 内部十余处直接用 `GetDEMParam()->GetTimeStep()`。
- **内容**：把 dt 提升为可覆盖量（`SetTimeStepOverride(double)`，或把 dt 变成 `SolveOneStep(dt)` 的参数）。
- **验收**：覆盖值与直接写 `dem.time_step` **逐位一致**（否则会引入不可见的行为漂移）。

### 1.3 时间步协调（sub-cycling）——必须有
- **量级**：CFD 的 CFL dt（12 GPa 气体、dx = 0.1 mm 时 ~1e-8–1e-9 s）
  与 DEM 的稳定 dt（`dem.time_step = 5e-10 s`，由接触/键频率与断裂事件定）**差 2–100 倍**。
- **方案**：外层 CFD 主步 `dt_cfd`；内层跑 `N = ceil(dt_cfd / dt_dem)` 个 DEM 子步；
  气动力在子步内**冻结**（交错一阶）。气相比热/压力在子步内不变。
- **关键约束（阶段 4 的教训）**：**散列进流体的力与推进颗粒的力必须是同一份**，
  否则动量/能量账本不闭合（`docs/TWO_PHASE_STAGE4_TWO_WAY.md` §3.3 的 `m_results_pre` 设计）。
- **验收**：
  - `N = 1` 时与"等步长"逐位一致；
  - `N > 1` 时能量/动量账本仍闭合；
  - `N` 加倍结果收敛（一阶，误差 ∝ N·dt_dem）。

### 1.4 ε_g 的重建与所有权
- **内容**：DEM 粒子动过之后每步重建 ε_g 并注入；`CavityGasModel` 的"边界环"概念在新路径里取消
  （至少不作为压力来源）。
- **验收**：`Σ V_p = Σ (1 − ε_g)·V_cell`，残差 ≤ 1e-12（映射器本身能到 3.5e-16）。

---

## 四、工作块 2：气固相互作用力（物理最关键）

### 2.1 自由表面压力（面力）——**必须新写**
- **为什么不能只用 `−V_p∇p`**：`−V_p∇p` 是"p 在粒子尺度上线性"的一阶近似。
  空腔壁与裂纹面是**自由表面**，压力在这里从 p 跳到环境值；
  激波/Hugoniot 面前 `λ_p ~ d_p`，一阶近似直接失效。
  现有的 0 维气腔（边界环面力 `p·ℓ·t`）在这一点上**反而比 `−V_p∇p` 更准**。
- **正确做法**：对与气体直接相邻的粒子做压力面积分
  `F_i = Σ_faces p_f · A_f · n̂_f`。
- **实现抓手（好消息）**：`ParticleGridMapper` 已经能算出
  「粒子被切成的小块各自落在哪个单元、占多大体积」——
  这套权重**同时就是"粒子暴露给气体的面积与方向"**。
  `Exact` 方法的球–盒交集体积、`SubCell` 的八分体都能直接转成面元。
- **参照**：`inc/Basic/DEM/CavityGasModel.h`（边界环面力 + 散度定理算体积）、
  LSM `gas_phase_mod.f90` 的 `force_between_gas_and_spring_wall`（气楔入裂纹）。
- **验收**：
  1. 静止均匀压力场中，任意粒子的面力合力 ≡ 0（散度定理，机器零）；
  2. 单粒子在已知线性梯度场中，面力与 `−V_p∇p` 的相对差 = O((d_p/λ_p)²)；
  3. 静弹性变体上对 Lamé 解复核（见 4.3 V-A）。

### 2.2 曳力（已有，但要处理刚性）
- 5 个模型可直接复用；`β` 与 `β/ε_s` 差一个 ε_s 的坑已在阶段 3 踩过并有校验。
- 需要新增：ε_g 很小时的**上界/亚松弛**（见 3.4）。

### 2.3 压力梯度力的对偶性复核
- `coupling.react_pressure_gradient_force`（`auto/on/off`）与
  `two_phase.porosity_gradient_force` 必须配套，配错会双计或漏项（现有代码会告警）。
- 在新链路里这条恒等式要重新验一遍：`F_drag + F_∇p ≡ (−dp/dz)·V_p/ε_s`。

### 2.4 相间功（含面力做功）
- 已有 `coupling.react_interphase_work`（`S·u` 进能量方程）。
- ⚠ **面力也要做功**（`F_surf·u_p`），不能只算曳力的功，否则能量账本会缺项。

---

## 五、工作块 3：CFD 侧为高压/爆炸场景补齐的能力

### 3.1 气体状态方程：JWL 或多方
- **现状**：只有 `PerfectGas`。12 GPa 下理想气体完全失真（温度会算出 1e5 K 量级）。
- **内容**：给 `Gas` 抽象基类加三个 EOS 钩子
  `p(ρ,e)` / `e(ρ,p)` / `c(ρ,e)`，然后实现 JWL：
  `p = A(1 − ωρ/(R₁ρ₀))·exp(−R₁ρ₀/ρ) + B(1 − ωρ/(R₂ρ₀))·exp(−R₂ρ₀/ρ) + ωρe`
- **参数**：可复用 `dem.jwl_a/b/r1/r2/omega`（HMX 基 PBX 的一组值已在 `DEMSolverParam.h:179-183`），
  或新开 `[gas.jwl]` 段避免与 DEM 侧耦合。
- **兜底**：若只需"定性撑碎"，多方气体 `p = (γ−1)ρe` 就够（理想气体已支持可配 γ），
  但要注意 `γ` 与真实气体的对应关系。
- **验收**：
  1. 静止均匀场解析解；
  2. 与 0 维 `CavityGasModel` 的绝热膨胀 `pV^γ = const` 对照（同一初值、同一膨胀历史）。

### 3.2 初始条件：空腔高压区
- **已有**：`init.explosion.{dim, center_x/y/z, radius, inner/outer_density, inner/outer_pressure, inner/outer_velocity_x/y/z}`
  （`EulerSolverStructUniform.cpp:232-309`）。`dim = 2` 正是"只用 (x,y) 距离"的柱形，
  与二维圆盘 DEM 晶格天然匹配。
- **要确认**：文件里取值都是**无量纲**（`p* = p/(ρ_ref a_ref²)`），
  12 GPa / 0.40 GPa 需要先除以 `p_ref`。`coupling_report.csv` 会写参考量，可直接取。
- **验收**：t = 0 帧上，空腔节点 p 精确等于目标值、外部等于环境值（用 `result/0.dat` 逐点核）。

### 3.3 边界条件：远场 / 低压环境 / 计算域尺寸
- **现状**：只有 `outlet`（零梯度外推）/ `inlet`（常值来流）/ `wall`（镜像反射），**无随时间变化边界**。
- **需要**：
  - 气体喷出方向用 `outlet`；
  - 计算域要留净空（现有经验：碎片速度 ~1.7 km/s、`t_end ~1e-5 s` ⇒ 至少 2R 净空，D = 20 mm 时盒 ±30 mm）；
  - ⚠ **不要让 `outlet` 外推到真空**：理想 Euler 的零梯度外推在"喷入真空"时会产生负压/负密度。
    外部给一个不低的低压（0.1–1 MPa 或大气压），并加"p、ρ 必须为正"的检查。
- **验收**：无气体对照（ε ≡ 1、p 均匀）下全场 `max|v|` 只剩舍入漂移（阶段 4 的口径：严格为 0 或 ~1e-14）。

### 3.4 固相密排区的数值稳健性 —— **风险最高**
- **量级**：二维三角密排 `ε_g ≈ 1 − π/(2√3) ≈ 0.093`；离散化后单元内可能低到 0.05。
- **两个刚性来源**：
  1. **曳力**：`β/ε_s` 在 Ergun 支 `∝ ε_s/ε_g²`，ε_g = 0.09 时比稀疏区大 3 个量级；
   `k = V_p(β/ε_s)/m` 随之爆炸。好消息是 `AdvanceParticles` 的 `exponential` 积分器
   **对冻结线性曳力无条件稳定**，颗粒侧不会炸；坏消息是**散列进流体的源项 `S` 是显式的**，
   会在 ε_g 小的单元里把动量方程打死。
  2. **`p∇ε` 源项**：`ε_g` 从 0.09 跳到 1.0 的界面上，`∇ε` 很大，源项是刚性的。
- **候选手段（按代价递增，建议先 (a) 再 (c)）**：
  - (a) 对 `β/ε_s` 加 ε_g 相关的**上界**（例如 ε_g < 0.2 时改用 Ergun 的层流线性项并限幅）；
  - (b) 对源项 `S` 做**半隐式**（至少把曳力的线性部分隐式化，与 `exponential` 积分器配套）；
  - (c) **固相区（ε_g < ε_c）退出 CFD**：在 ε_g 极小处把气固作用退化成"面力 + 封闭腔"，
    即实质上是"CFD 只在气体连通的区域求解，固体内部用面力"。
- **验收**：1e4 步无 NaN；能量账本闭合；全程 `p*`、`ρ*` 为正；`ε_g` 无越界告警。

### 3.5 网格分辨率 —— 有一个**必须正视的结构性问题**
- **需求估算**：裂纹宽度 ~ `L0 = 0.25 mm`。要"看见裂纹"至少 `dx ≲ L0/2 = 0.125 mm`。
  D = 20 mm + 外围（±30 mm）⇒ 60 mm / 0.125 mm = 480 格 → `480×480×3 ≈ 6.9e5` 节点。
  规模本身可接受（均匀网格、无 AMR，内存几十 MB/场），但要评估运行时间。
- ⚠ **真正的麻烦是 `Δ/d_p` 落在最坏区间**：
  - 粒子直径 `d_p = 2r = L0 = 0.25 mm`；若 `dx = 0.125 mm` ⇒ `Δ/d_p = 0.5`；
    若 `dx = 0.25 mm` ⇒ `Δ/d_p = 1`。
  - 工程里的经验是"体积分数法要求 `Δ/d_p ≳ 3`"（只在 `DEMCFDCoupler.cpp:111-124` **告警**，不报错）。
  - 也就是说：要么 `Δ ≫ d_p`（体积分数法中值化，`Δ/d_p ≳ 3` ⇒ dx ≥ 0.75 mm，
    **根本分辨不了裂纹**），要么 `Δ ≪ d_p`（用 `exact` 映射 + 面力，退化成"一个粒子跨多格"）。
    **现在恰好卡在中间两头不讨好。**
  - ⇒ 这需要一次**专门的分辨率研究**（在 `Δ/d_p = 0.25 / 0.5 / 1 / 2 / 4` 各跑一次，
    看碎块形态与能量账本是否收敛）。若收敛性不成立，应转向 §七的备选路线。
- **验收**：`dx` 减半、`dt` 减半，碎块统计（最大碎块质量占比、裂纹带大小、径向断键率）收敛。

---

## 六、工作块 4：算例、输出、验证

### 4.1 算例生成器
- `tests/hmx_cavity_burst/generate_case.py` 增加 CFD 网格生成：
  写 `mesh.dat`（节点坐标）与 `mesh.inp`（块尺寸 + 六个面的边界类型），
  并把 `[freestream]`、`init.explosion.*`、`[coupling]` 一起写进 `zaran.toml`。
- 注意二维化：`nk = 3` + `z` 方向 wall，或 `dim=2` 的 `init.explosion` —— 与现有
  `tests/demcfd_ball_riemann`（`81×81×3`）同一套路。

### 4.2 输出
- 流场：`result/<iter>.dat`（Tecplot ORDERED，已有 `X,Y,Z,Density,Velocity_*,Pressure,Volume_fraction`）
  —— **这正是你一直想看而没有的东西**（当时纯 DEM 路径根本没有流场）。
- 粒子点云 + 键：已有。
- 建议新增：`gas_energy.csv`（气相内能、气固功交换、压力极值）—— 能量账本要靠它。

### 4.3 验证组（新增 6 组，每组都要有判据）

| 编号 | 内容 | 判据 |
|---|---|---|
| V-A | 静弹性场对 Lamé 解（压力由 CFD 提供，取代 0 维气腔） | 径向/环向应力与 Lamé 解误差 < 5%（沿用现有口径：抬阈值 + 慢升压 + 加粘性） |
| V-B | 能量账本 | `U_gas0 + W_ext = U_gas + KE + ΣE_bond + ΣE_断裂 + ΣE_耗散 + ΣE_接触`，残差 < 1% |
| V-C | 与 0 维气腔的**早期时刻**对照 | 裂纹尚未形成（`t ≲ 起裂时刻`）时，两者空腔压力历史相对差 < 5% |
| V-D | 退化对照 | `ε ≡ 1` / 无气体 ⇒ 流场冻结；关闭耦合 ⇒ 与纯 DEM 逐位一致 |
| V-E | 收敛性 | `dx` 减半、`dt` 减半，碎块统计收敛 |
| V-F | 回归 | `_demtest/run_committed_tests.py` 的 **22 项不能坏** |

---

## 七、需要你定的 6 件事

| # | 问题 | 建议 |
|---|---|---|
| 1 | 气体模型：多方/理想（γ 可配）够，还是必须 JWL？ | 先多方跑通流程，再上 JWL |
| 2 | 维度：2D 薄板还是真 3D？ | **2D 薄板**（与 DEM 晶格一致、省 ~3 倍） |
| 3 | 外边界：开放自由空间还是封闭腔？外部环境压力？ | 开放 + 环境压力 0.1–1 MPa（避免喷入真空） |
| 4 | 目标压力档：0.40 GPa（6.8×p_crit）还是 12 GPa？ | **先 0.40 GPa**（压力比 4e3 而非 1.2e5） |
| 5 | 网格预算：粗网格先跑通（`dx = L0`，~1.7e5 节点）还是直接上细网格？ | 先粗网格跑通，再做 V-E 收敛 |
| 6 | `Δ ≈ d_p` 最坏区间怎么处理？ | 先做专项分辨率研究；不收敛就转 §八 备选路线 |

---

## 八、备选路线：Voronoi 气相场（LSM 的做法，值得同时评估）

`docs/CAVITY_GAS_BURST_GAPS.md` 的 B 组已经列过：

| 步骤 | 内容 | 工程现状 |
|---|---|---|
| B1 | 解除 `UpdateGasVoronoiMesh` 的 `reaction_enabled && voronoi_enabled` 门控 | 代码**已存在**，只是被门控 |
| B2 | 支持"事先存在的高压气体"初值（`p₀, T₀` 逐格点，不依赖 `reaction_progress`） | `JwlPressure` 要求 `reaction_progress > 0` |
| B3 | 固相破碎 → 气相单元的**守恒转移**（体积/内能/质量） | 无 |
| B4 | 气相 → 固壁（**含裂纹面**）作用力 = 气楔入 | 无 |
| B5 | 气相格点间压力平衡 / 导热 | 无 |
| B6 | 气相场输出 | 无 |

**优点**：网格自动跟随粒子，不受 "Δ vs d_p" 困扰；裂纹里的气体**自动**有单元；
与 LSM 可以直接逐项对照（LSM 就是这么做的）。
**缺点**：不是教科书 CFD，流动精度低于有限体积。

**建议**：工作块 1（架构合并）与工作块 2（面力）对**两条路线都通用**，先做；
把 Voronoi 路线保留为"无法在均匀网格上分辨裂纹"时的退路。

---

## 九、建议的实施顺序（4 步，每步都有可验收产物）

| 步 | 内容 | 产物 |
|---|---|---|
| **Step 1** | 合链路（最小可跑）：新驱动 + dt 覆盖 + 等步长 + 体积力 + 多方气体 | 静弹性变体对 Lamé 解；两个场同时输出 |
| **Step 2** | 时间与力的正确性：sub-cycling + 自由表面面力 + 相间功（含面力做功） | `N=1` 逐位一致；`N` 收敛；面力散度定理；Lamé 仍成立 |
| **Step 3** | 高压场景：JWL/多方 + 初始高压区 + 远场边界 + ε_g→0 刚性处理 | 0.40 GPa 档跑通；能量账本闭合；无 NaN；有流场云图 |
| **Step 4** | 验证与工程化：V-A…V-F + 分辨率研究 + 回归 | 完整验证报告；22 项回归不坏 |

**风险最高的两项（如果这两项过不去，整个路线要重估）**：
1. **ε_g→0 的刚性 + `Δ/d_p` 落在最坏区间**（§3.4、§3.5）；
2. **面力的正确性**（§2.1）—— 如果继续只用 `−V_p∇p`，等于白加 CFD。
