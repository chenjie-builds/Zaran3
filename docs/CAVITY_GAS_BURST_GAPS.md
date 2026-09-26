# 高压气腔撑碎算例：力学模型缺口清单与实施建议

> **状态更新（同日）：A 组核心已实施并验证完毕。**
> 详细设计与验证见 **`docs/COHESIVE_FRACTURE_MODEL.md`**。已完成的项：
>
> | 项 | 内容 | 状态 |
> |---|---|---|
> | A1/A2 | γ 控制的双线性内聚律（升-峰-软化-干净断裂，δ_f = 2·Efract/F_p 由面积反解） | ✅ 已实现 |
> | A3 | 能量阈值口径统一（不再"储能 ≥ Efract"瞬断；表面能单独记账） | ✅ 已实现 |
> | A3b | γ ↔ σ_t ↔ 网格自洽性检查（γ_min、l_cz/L0，不满足即报错） | ✅ 已实现 |
> | A6 | 断键能量分配 | ⛔ **经核实不需要** —— 见下方更正 |
> | A10 | 内聚律参数输出 + 碎片统计 + 拉格朗日径向损伤剖面 | ✅ 已实现 |
> | D1/D2 | 载荷量级自检（起裂阈值扫描）+ 时间步收敛判据 | ✅ 已实测 |
> | A7/A8/A9 | 卸载规则 / 单元级破损 / 初始缺陷场 | ⏳ 仍未做 |
> | B / C2 | Voronoi 气相场、气相输出、碎块黏附/再键合 | ⏳ 仍未做（本轮限定纯 DEM） |
> | D3 | 载荷速率可控（升压时间）与准静态/动态两条验证路径 | ⏳ 仍未做（`ramp_time = 0`） |
>
> **状态更新（第二轮）：A4/A5 + C1/C3 已实施并验证完毕。**
> 详细设计与验证见 **`docs/MULTICHANNEL_FRACTURE_AND_CONTACT.md`**。
>
> | 项 | 内容 | 状态 |
> |---|---|---|
> | A4 | 切向独立断裂 + 部分失效（`Et ≥ UtIII` ⇒ 切向先断，按 LSM `GBratio` 折减法向**能量预算**） | ✅ 已实现 |
> | A5 | 压剪通道（压缩侧 `Et ≥ UtIII` 独立成通道，与"压缩压溃"分开） | ✅ 已实现 |
> | A5b | 拉伸混合通道（`En + Et ≥ Un_eff`） | ✅ 已实现 |
> | C1 | 滚动阻力矩（`M_r = −μ_r·R*·|F_n|·ω̂_rel`，含过冲截断） | ✅ 已实现 |
> | C3 | 切向接触阻尼（`c_t = λ·c_n`，与弹簧力合并后统一受库仑上限约束） | ✅ 已实现 |
> | C4 | 碎片统计输出 | ✅ 第一轮已实现（`fragment_stats.csv`） |
>
> 四项新功能**全部默认关闭**（`dem.bond_shear_energy_ratio = 0`、
> `dem.rolling_friction = 0`、`dem.tangential_damping_scale = 0`），
> 所以 22 项（19 历史 + 3 新增）回归逐位不变。
> 在本文档的**标准撑碎档**（0.40 GPa）上做开关对照后：
> 最大裂纹带 220 → 1086 粒子（×4.9）、径向断键率从"中段回升"变成**单调外减**、
> 断键率 4.63% → 10.18% 但最大完整碎块质量占比仍从 30.9% 升到 33.3%
> ⇒ **裂纹真正连通、碎得"有结构"而不是全盘粉化**（详见 `MULTICHANNEL_FRACTURE_AND_CONTACT.md` §6.4）。
>
> ### ⚠️ 一处更正（原先判断有误）
> 我在上一版把"断键不产生飞溅"列为根因之一（A6）。读了 LSM 的实现后更正：
> `break_spring_energy_dissipation.f90` 里 `Ubreak = UnIII − Efract`，
> 而 `UnIII = Efract + Estore`、`Estore = UnI − Ecr`、`Ecr = UnI·(r00/r0)`。
> **未反应、未压实（r00 = r0）时 Estore ≡ 0**，于是 `Ubreak = 0` ⇒ `Ekbreak = 0、Qbreak = 0`，
> **LSM 在标准配置下也完全不产生速度增量**。那个 kick 是绑在"气核收缩"机制
> （`magnify = (r0/r00)^Di`）上的，不是通用断裂机制。
> 因此**不需要**人造 kick：正确的做法是让断裂能 γ 控制"造新表面 vs 碎片动能"的分配，
> 飞溅是**涌现**的（实测标准档碎片速度峰值 2984 m/s）。

> 下文为原始清单（保留作为完整对照），其中"现状"列可能已被上面的状态表覆盖。
>
> 本文回答两个问题：
> (1) 当前圆盘"撑碎"到底是怎么算的、为什么没有流场网格输出；
> (2) 为什么现在的结果没有裂纹扩展、碎得很散，以及**还需要补哪些力学模型**（对照 LSM）。
>
> 参照对象：`LSM/code_utf8/`（Fortran）。经核实，`calculation.f90` 里
> `force_between_bonded_lattices` / `force_between_discrete_lattices` 的调用**已被注释掉**，
> 当前 LSM 实际使用的是**弹簧阵列模型**：
> `covalent_ionic_spring_mod.f90`、`elastic_plastic_spring_mod.f90`、
> `general_ineteraction_spring_mod.f90`，气相为 `gas_phase_mod.f90`。

---

## 一、当前算例算的是什么（为什么没有流场网格输出）

`tests/hmx_cavity_burst/zaran.toml` 里是：

```toml
[task]
simulation = "SOLVE_FIELD"
solver     = "DEM"
```

即走 `Application` 的 `DEM` 分支 → `DEMFieldSimulation`。这条路径**只做 DEM（键合弹簧网络）**，
**根本不创建任何流体/网格求解器**（单相 `Euler`、两相 `EulerTwoPhase`、耦合 `EulerDEM` 都没有被实例化）。
因此 `result/` 里只有三个东西：

| 文件 | 内容 | 是否"网格" |
|---|---|---|
| `particles_<iter>.dat` | 5760 个节点的点云（FEPOINT），逐点带 `velocity_*`、`temperature`、`gas_pressure`、`gas_volume`、`damage` 等 | 点云，**无 I/J/K** |
| `bonds_<iter>.dat` | 16974 条键的线元（FELINESEG），带 `active`、`damage`、`strain`、`force_*`、`elastic_energy` | 线元，**无 I/J/K** |
| `cavity_gas.csv` | 本算例新增的 0 维气腔诊断（p、V、U、能量账本） | 表格 |

### 气体是怎么处理的

没有气相网格。`CavityGasModel`（本次新增）是一个 **0 维集中参数模型**：

- **加载面**：把空腔壁离散成"边界粒子环"（30 个粒子），每个粒子受力
  `F_i = p · ℓ_i · t · n̂_i`（`ℓ_i` 是角向分片宽度、`t` 是面外厚度）；
- **体积**：由同一套边界离散用散度定理给出，`V = ½|Σ (r_i × r_{i+1})|`；
  受力与体积共用同一次几何刷新 ⇒ `ΣF_i·Δu_i ≡ p·ΔV` 逐位恒等；
- **状态**：`U ← U − pΔV`、`p = (γ−1)U/V`（绝热，`γ=1.4`），U₀ = 21.97 J；
- **假设**：**整个空腔只有一个压力**，且这个压力均匀作用在边界环上。

### 为什么没有用框架里已有的 Voronoi 气相

Zaran3 **确实**有一套节点中心的气相场（`DEMSolver::UpdateGasVoronoiMesh()` 建 Delaunay/Voronoi
单元与内部面，逐粒子 `gas_volume / gas_pressure / gas_temperature / gas_internal_energy`，
压强用 JWL 状态方程，还有气相间导热）。但它被**硬门控**在反应路径上：

```cpp
void DEMSolver::UpdateGasVoronoiMesh(bool force) {
    if (!GetDEMParam()->GetReactionEnabled() || !GetDEMParam()->GetVoronoiEnabled())
        return;                       // ← 惰性算例直接退出，一个单元都不建
    ...
}
double DEMSolver::JwlPressure(...) const {
    if (p.reaction_progress <= 0.0 || ...) return 0.0;   // ← 未反应的粒子压力恒为 0
    ...
}
```

本算例 `reaction_enabled = false`、且气体是**事先存在的惰性高压气体**（不是炸药反应产物），
所以这套气相场用不上——**这就是"没有流场输出"的根本原因，也是必须补的第二大块能力**。

### LSM 的做法（对照）

LSM 的 `gas_phase_mod.f90`（2978 行）就是把气相建在格点上：

| LSM 子程序 | 作用 |
|---|---|
| `record_new_gas_lattice` | 固相单元破碎后，**把它转成新的气相格点**（体积、内能守恒转移） |
| `remove_solid_lattice` | 移除固相格点、更新气相占据 |
| `build_gas_cell_list` / `bulid_Voronoi_polygon` | 在格点上建 Voronoi 元胞（气相体积 = Voronoi 面积 × 厚度） |
| `calculate_Voronoi_polygon_state` | 逐格点状态：`P`、`T`、`E`、比容、膨胀功；初始气相用 `gasC(1,n)=P`、`gasC(2,n)=T` |
| `force_between_gas_lattices` | 气相格点之间的压力/黏性/导热交换 |
| `force_between_gas_and_spring_wall` | 气相 → 弹簧壁面的作用力（**裂纹里的气体楔入**） |
| `outputData_with_gasInformation_mod.f90` | 输出气相信息（就是用户想看到的"流场"） |

⇒ **LSM 会把气体"灌进"裂纹并推动裂纹继续扩展；我们的 0 维气腔做不到这件事。**

---

## 二、诊断：为什么没有裂纹扩展、而且碎得很散

### 硬数据（用 `_demtest/diag_crack.py` 从 `bonds_<iter>.dat` 反算）

| iter | t | 断键率 | 最大**完整**碎块（粒子数） | ≥20 粒子碎块 | 孤立单粒子 | 处于软化段(0<d<1)的键 |
|---|---|---|---|---|---|---|
| 0 | 0 | 0.00% | 5760 | 1 | 0 | 0 |
| 100 | 1e-7 s | 1.32% | 5694 | 1 | 60 | 1 |
| 500 | 5e-7 s | 7.87% | 5404 | 1 | 331 | 118 |
| 2000 | 2e-6 s | **56.09%** | 3645 | 1 | 1979 | 1191 |
| 10000 | 1e-5 s | **94.04%** | **26** | **4**（仅占 1.7% 粒子） | **4309** | 537 |

**断键率的径向分布**（12 个等宽环带，从圆心到外缘）：

```
iter=100   : 0.0%  76.2%  0.0%  0.0% ...            ← 裂纹在孔壁起裂 ✔
iter=500   : 0.0%   0.0% 99.1%  56.6%  2.9%  0.0%... ← 一个"环状破碎带"
iter=2000  : 0.0%   0.0%100.0% 100.0%  99.7%  99.7%  98.8%  71.2%  30.3%...
iter=10000 : 100%  100%  100%  100%  100%  100%  99.0%  92.9%  89.1%  86.4%...
```

**读法**：裂纹**确实从孔壁起裂并向外走**（iter=100 时第一环 76%），但它不是一条"细裂纹"，
而是一个**越来越厚的破碎环带**——到 iter=2000，从圆心到 r/R≈0.58 已经 100% 断光。
最终 94% 的键断掉、最大完整碎块只剩 **26/5760 = 0.45%** 粒子、**4309 个孤立单粒子**，
即"**全盘粉化**"，而不是"几大块 + 飞溅细粒"。

### 五条根因

**(1) 载荷严重超驱动 —— 整个圆盘从 t=0 就全部在破坏包线之外。**
Lamé 解孔壁环向应力 `σ_θ(a) = p(R²+a²)/(R²−a²) = 1.0202·p₀`；按 W2 反推的起裂压力
`p_crit ≈ 118–154 MPa`，而算例给了 `p₀ = 12 GPa`，即 **p₀/p_crit ≈ 78–101 倍**。
12 GPa 下孔壁 12.24 GPa、外缘仍有 0.245 GPa，全部 ≫ HMX 抗拉强度(~0.1 GPa)。
**裂纹扩展需要"裂尖超限、远场弹性"的应力梯度**；现在是全场上限，所以只能同时、全域、弥散地粉碎。

**(2) 断裂能没有被真正表达 —— 碎块尺寸不可控。**
本算例用的是**应变**口径：`bond_peak_strain = 1.5e-3`、`bond_break_strain = 3e-3`，
软化段下面积由 `kn·(ε_break−ε_peak)` 决定，**与材料的表面能 γ 无关**。

框架里其实有两种口径，**都不等于 LSM 的完整内聚律**（`src/Zone/DEM/DEMSolver.cpp:956/971`）：

| 口径 | 触发条件 | 行为 | 问题 |
|---|---|---|---|
| γ 口径 | `dem.surface_energy > 0` | 阈值 `½k·ext² ≥ Efract·A·s²`，**达阈值即瞬时断裂并 `continue`** | **完全没有软化段**；且"储能 ≥ 表面能"比 LSM 的 `≥ Efract + Estore` 少一项；断键后不扣表面能 |
| 应变口径（本算例） | `fracture_energy = 0` | 双线性软化 peak→break | 曲线面积与 γ 无关 ⇒ 碎块尺寸不可控 |

LSM 是反过来的：`Efract = gamma2 · surface`（`surface = √3·r₀/3`），
`UnIII = Efract + Estore`，并由 `Efract` **反解出软化终点 `rnmax`**，保证力-分离曲线下的面积
= 表面能 + 储能。**γ 决定"多少能量用于造新表面、多少变成碎片动能"，直接决定碎块尺寸分布。**
现在两种口径都缺这条闭环 ⇒ 再调网格也不会收敛到真实碎块数。

**(3) 断裂判据只有"法向应变"一条通道，且没有渐进/部分失效。**
LSM 每个键有四条通道（`covalent_ionic_spring_mod.f90:489-527`）：

| LSM `breakmod` | 判据 | Zaran3 是否有 |
|---|---|---|
| 1 拉伸（+剪切混合） | `En + Et ≥ UnIII` | 只有应变口径的近似 |
| 2 **压剪** | `rn ≤ r0` 且 `Et ≥ UtIII` | ✘ 无 |
| 3 **压缩压溃** | `rn < rnmax·0.6` | 有（ratio=8 的脆性阈值），但口径不同 |
| 4 **切向弹簧部分断裂** | `Et ≥ UtIII` → 切向失效，并把剩余法向阈值按 `GBratio` **折减** | ✘ 无 |

`breakmod=4` 的"部分失效 + 阈值折减"是裂纹能"走走停停、沿路径择优"的关键；
现在没有它，破坏一旦发生就是同时、全域、不可逆。

**(4) 断键不产生飞溅 —— 这直接对应"没有飞溅的颗粒"。**
LSM `break_spring_energy_dissipation.f90` 的断键能量分配：

```
Ubreak = UnIII − Efract          ! 扣除表面能后的剩余能量
Eklimit = 2·(0.5·m·vlimit²)      ! vlimit = max(10, 0.1·vmax)，动能上限
if Ubreak ≤ Eklimit: Ekbreak = 0.5·Ubreak,  Qbreak = 0
else:                Ekbreak = 0.5·Eklimit, Qbreak = 0.5·(Ubreak − Eklimit)
vi = +√(Ekbreak/mi), vj = −√(Ekbreak/mj)   → 给两端加减速度（速度增量）
```

即**一半残余能量变成两端节点的速度增量**（沿法向/切向按 `rnmax`、`rtmax` 配比），
一半记热。Zaran3 现在是**把全部断键弹性能记成热**（`dissipated_fracture_energy`），
两端速度**一点都不改** ⇒ 碎片只有几何排斥带来的低速散开，没有"爆开"的颗粒。

**(5) 碎后接触太"滑"、且没有气体回压。**
碎片之间只有 `LinearSpringDashpot`：法向（含重叠保护）+ 库仑切向（μ=0.3），
**没有滚动阻力矩**，没有黏附/再键合，没有碎块锁固。更关键的是：
0 维气腔模型**无法把气体压力作用到裂纹面上**，而 LSM 的
`force_between_gas_and_spring_wall` 正是靠"气楔入裂纹"把外侧碎块撑住/复位。

**(6) 缺初始缺陷/损伤场 —— 没有裂纹起始点与取向。**
LSM 有 `set_damage_initial_condition.f90`、`set_distribution_accordingTo_ignitionSite.f90`
（`damage ∈ [0,2]`，`>1` 即预破损，参与建链与刚度）来给出裂纹萌生位置。
Zaran3 只有逐键 Weibull 强度散布（m=8），**空间上不相关**，只会把破坏弥散化。

**(7) 没有单元级破损记录。** LSM 一旦断键就把**整个三角形单元**标记破损
（`break_node=1`、`element(4,j)=0`、`break_mod` 记 1/2），并把该单元转成气相格点。
Zaran3 只有逐键 `active`，没有"单元"和 `break_mod`，所以既无法保证裂纹路径连通，
也没有裂纹/碎片单元这种可导出、可统计的载体。

---

## 三、待实施清单（按优先级；**本次不实施**）

### A 组：断裂本构（决定"有没有裂纹"与"碎成几块"）——最高优先级

> **先更正一处**（`src/Zone/DEM/DEMSolver.cpp:178-182`、`1252`、`1317`）：γ 口径**不是**死代码——
> `dem.surface_energy` 与 `dem.grain_boundary_energy` 都实现了，`bond.fracture_energy`
> `= γ·L0·t/√3·(1 − E_gb/γ)`（与 LSM 的 `Efract = gamma2·√3·r₀/3` 在 `r₀=L0` 时完全一致），
> `hmx_disk_10k_center_ignition` 就在用（γ=1680 J/m²）。
> **但本算例没有给 `surface_energy`，于是回退到应变双线性口径。**
> 更要紧的是这两种口径**都不等于 LSM 的完整内聚律**，见 A1–A3。

| # | 要实施的内容 | 现状 | LSM 参照 | 必要性 |
|---|---|---|---|---|
| A1 | **完整牵引-分离（内聚）律**：非线性硬化段 → 峰值 `(rncr, fncr)` → 线性软化 → 干净断裂；参数 `γ`、`fncr`、`rncr`、`rnmax`、`nn`、`Kn1/Kn2/Kn3` | 只有两种**互斥**的简化：<br>(a) `fracture_energy>0`（γ 口径）⇒ **瞬时脆断、无软化**；<br>(b) `fracture_energy=0`（应变口径）⇒ 双线性软化、**无断裂能控制**（本算例走的这条） | `covalent_ionic_spring_mod.f90:186-262`；`read_material_properties.f90:435-473` | 让"键断裂能 = γ·A + 储能"与软化段**同时**成立，碎块尺寸才可标定/可收敛 |
| A2 | **按 γ 反解软化终点**，使力-分离曲线下面积严格等于 `Efract + Estore` | 无（`rnmax` 无对应物） | 同上：`dU = Efract − Ecr` → `rnmax`，并以 `r0max = 1.7·r0` 封顶 | 这是"断裂能守恒"的唯一正确写法 |
| A3 | **能量阈值补 `Estore` 项**。现在判据是 `½k·ext² ≥ Efract·A·s²`，即"**储能 ≥ 表面能**"；LSM 是 `En ≥ Efract + Estore`（硬化段累积面积）。且断键后 Zaran3 把**全部**储能记热、**不扣表面能** | 阈值口径差一个 `Estore`；表面能没有单独记账 | `UnIII = Efract + Estore`；`break_spring_energy_dissipation` 里 `Ubreak = UnIII − Efract` | 否则"造新表面消耗多少能量"根本没有被表达，碎块尺寸不可控 |
| A3b | **γ 与 σ_t 的一致性约束**。文档已记：HMX `γ=8.4 J/m²` 与 `σ_t=60 MPa` **只在 `L0 ≈ 70–130 μm` 时自洽**；本算例 `L0 = 250 μm`，直接用 γ 会得到与 `σ_t` 不匹配的强度 | 需人工试算 | 由 `Efract`、`fncr`、`rncr` 三者关系反解 | 换 L0 或换 γ 时强度会漂移，必须给出选择准则（先定 γ 与 σ_t ⇒ 反解 L0） |
| A4 | **切向独立断裂 + 部分失效**：`UtIII = α·UnIII`；`Et ≥ UtIII` 时切向先断，按 `GBratio` 折减剩余法向阈值 | 无 | `covalent_ionic_spring_mod.f90:489-527` | 渐进损伤；裂纹路径择优的关键 |
| A5 | **压剪通道**（`rn ≤ r0` 且 `Et ≥ UtIII`）与**压缩压溃通道**（`rn < rnmax·0.6`）分开 | 压缩只有单一 ratio 阈值 | 同上 | 孔洞算例里压缩/剪切破坏占主导（现能量账本里压缩失效贡献远大于拉伸） |
| A6 | **断键能量分配**：`Ubreak − Efract` 一半转两端速度增量（带动能上限 `vlimit`），一半记热 | 100% 记热 | `break_spring_energy_dissipation.f90` | **直接决定"飞溅颗粒"是否存在**，用户已观察到此缺失 |
| A7 | **卸载/再加载规则**（软化段卸载线性回原点、裂纹面反复开合耗能正确） | `damage` 单调、刚度按 (1−d) 折减，卸载不严格 | 各段的解析 `En` 表达式 | 碎块碰撞/回弹阶段能量正确 |
| A8 | **单元级破损记录与 `break_mod` 输出**（三节点 1/2/3 分类），并支持把破损单元标记为"裂纹/气相新增体积" | 只有逐键 `active` | `force_between_bonded_lattices.f90:228-274`；`element(4,j)=0` | 裂纹可视化、碎片统计、与气相的接口 |
| A9 | **初始缺陷/损伤场**（从文件读 `damage`/`matID`，或强度场带空间相关长度） | 只有不相关的 Weibull m=8 | `set_damage_initial_condition.f90`；`set_distribution_accordingTo_ignitionSite.f90` | 提供裂纹萌生点与取向；不相关噪声只会弥散化 |
| A10 | **输出 `Damage` 场与主裂纹长度时序** | `bonds_*.dat` 有 `damage`/`active`，但没有裂纹长度、没有单元破损 | `output_sample_reaction_file` | 定量验证"裂纹扩展速度 ~0.3–0.6·C_R" |

### B 组：气相（决定"能不能看到流场"与"裂纹能否被气体继续撑开"）——次高优先级

| # | 要实施的内容 | 现状 | LSM 参照 |
|---|---|---|---|
| B1 | **解除气相 Voronoi 的 `reaction_enabled` 门控**，让惰性算例也能建/输出气相场 | `UpdateGasVoronoiMesh` 直接 return | `gas_phase_mod.f90` |
| B2 | **支持"事先存在的高压气体"初值**：空腔区域按 `p₀, T₀` 初始化逐格点 `P/T/E/V_gas`；不依赖 `reaction_progress`、不用 JWL | `JwlPressure` 要求 `reaction_progress > 0` | `state%P = gasC(1,n)`、`state%T = gasC(2,n)`（`gas_phase_mod.f90:231-252`） |
| B3 | **固相破碎 → 气相单元的守恒转移**（体积、内能、质量） | 无 | `record_new_gas_lattice` / `remove_solid_lattice` |
| B4 | **气相 → 固壁（含裂纹面）作用力 = 气楔入** | 0 维气腔只能作用于"边界环" | `force_between_gas_and_spring_wall` |
| B5 | **气相格点间压力平衡/导热/黏性**（气相动量与能量在裂纹网络里重新分配） | 无 | `force_between_gas_lattices` |
| B6 | **气相场输出**（Voronoi 元胞的独立 zone，或增强点云字段） | 点云带 `gas_pressure/gas_volume` 但恒为 0 | `outputData_with_gasInformation_mod.f90` |
| B7 | 泄压/泄漏判据（孔洞与外界连通后压力释放的物理化） | 现有 `vent_area_ratio` 是几何粗判据 | `crevice_deepness` / `crevice_connection` |

### C 组：碎后接触（决定"看起来散不散"）

| # | 要实施的内容 | 现状 |
|---|---|---|
| C1 | **滚动阻力矩**（碎块的形状不规则性 / 角锁） | 无（已知缺口） |
| C2 | 碎块的黏附 / 再键合（可选，用于压实区） | 无 |
| C3 | 接触的率相关阻尼按恢复系数标定 + 切向阻尼 | 无切向阻尼 `c_t` |
| C4 | 碎片统计输出（碎片 id、尺寸、质量分布） | 需从 `bonds` 后处理，无原生字段 |

### D 组：驱动与数值

| # | 要实施的内容 | 说明 |
|---|---|---|
| D1 | **载荷量级自检**：`p₀` 应按 `p_crit` 的倍数给。撑碎/裂纹扩展用 `1.5–3×p_crit ≈ 0.2–0.4 GPa`；`12 GPa` 属于爆轰级，不是"撑碎级" | 本算例最直接的"结果不像"的原因之一 |
| D2 | 多判据时间步（`dt_vmax / dt_cl / dt_gas / dt_visco`）+ 断裂事件判据 `dt ≲ ε_t·L₀/v_rel` | 参照 `calculation.f90:2465-2474` |
| D3 | 载荷速率可控（升压时间）与准静态/动态两条验证路径 | 现在 `--ramp-time 0`（瞬时） |

---

## 四、配套验证方案（每项功能都要有判据）

1. **单键拉伸试验**：力-分离曲线与 LSM 解析分段式**逐点对照**；曲线下面积
   = `γ·A + Estore`（相对误差 < 1e-10）。
2. **断裂能网格无关性**：`L₀` 减半、`γ` 不变，碎块统计应**收敛**（现在的应变口径不该收敛，正好是对照组）。
3. **裂纹形态判据**：内部加压圆盘在 `p₀ ≈ 2·p_crit` 下应出现 **3–8 条主裂纹**，
   主裂纹连通分量长度 ≫ 键长，而不是现在这种"r/R<0.6 全断"的厚破碎环。
4. **裂纹扩展速度**：主裂纹长度(t) 的斜率应落在 `0.3–0.6·C_R`（HMX `C_R ≈ 2.4 km/s`）。
5. **碎块尺寸分布**：Weibull/Mott 分布，最大碎块质量占比与实验/文献对照；
   现在最大完整碎块 26/5760 = 0.45% 明显偏小。
6. **飞溅速度**：断键动能分配后，细粒速度应有明显高频尾（LSM 的 `Ekbreak` 口径可直接复算）。
7. **能量账本**（沿用现有 V2）：新增项 `表面能 + 断键动能分配` 后仍须闭合（现 0.368%）。
8. **气相一致性**：`Σ_cells V_gas = 空腔+裂纹总体积`；气相总内能随膨胀按 `Σ pΔV` 下降；
   气相与固壁的相互作用反作用力守恒（牛顿第三定律，参照阶段 4 的 W1）。
9. **回归**：19 项 DEM 回归必须保持 19/19。

---

## 五、建议的最小实施顺序

1. **D1 参数修正 + A6 断键动能分配 + A10 输出** —— 成本最低、立刻改善"飞溅"与可视化，
   且能马上用现有 25 项校验量出差异（V8 的速度上界判据要同步调整）。
2. **A1–A3 内聚律 + γ 闭环** —— 把 `Estore` 项、表面能记账、按 γ 反解 `rnmax` 补齐，
   并把软化段与 γ 口径统一起来（现在两者互斥）；同时给出 `γ ↔ σ_t ↔ L0` 的选择准则（A3b）。
   做完这一步，碎块尺寸才有物理含义。
3. **A4/A5/A8/A9 多通道失效 + 单元破损 + 初始缺陷** —— 才有真正的裂纹路径。
4. **B1–B6 Voronoi 气相场** —— 才有流场输出，也才能算"气体楔入驱动裂纹继续扩展"。
5. 最后再回到 `p₀ ≈ 2–3×p_crit` 的标准撑碎算例上做裂纹形态/速度/碎片分布验证。
