# 机械耗散生热：摩擦 / 阻尼耗散计入温度

> 用途：把**接触切向摩擦、接触法向阻尼、键合阻尼**不可逆耗散的能量计入粒子温度；
> 并把温度更新从 `reaction_enabled` 开关里解耦，使**惰性算例也会升温**。
> 参考实现：LSM `code_utf8/calculation.f90:2064` 与
> `code_utf8/force_between_discrete_lattices.f90`。

- 开关参数：`dem.mechanical_heating`（默认 `true`；置 `false` 恢复"温度只由键合导热与外部热源驱动"的旧行为）
- 改动文件：
  - `inc/Basic/DEM/DEMContact.h`（新增耗散字段）
  - `src/Basic/DEM/LinearSpringDashpot.cpp`、`src/Basic/DEM/HertzMindlin.cpp`（写入耗散）
  - `src/Zone/DEM/DEMSolver.{h,cpp}`（汇总耗散、温度解耦）
  - `inc/Zone/DEM/DEMSolverParam.h` / `src/Zone/DEM/DEMSolverParam.cpp`（开关）

---

## 1. 改动前的问题

1. **机械耗散完全不生热**：接触路径只把 `force_n + force_t` 累加到粒子并产生力矩，
   全代码没有任何 `F·v` 形式的耗散累加；LSM 的 `dissipation(1/2,·)` 数组**从未移植**。
   因非弹性而"消失"的机械能无声蒸发。
2. **温度更新被反应开关绑死**：`DEMSolver::CalcThermalReaction()` 首行是
   `if (!reaction_enabled) return;`，因此 `reaction_enabled = false` 的算例
   **温度完全冻结** —— 连键合导热与外部体热源都不生效。

## 2. LSM 的做法（参考基准）

`force_between_discrete_lattices.f90` 中把"格点摩擦功率"记为

```fortran
fnv = -visco1*vn/rn ;  ftv = -visco2*vt/rn
fn  = fn + fnv + fn3 ; ft  = ft + ftv          ! 此时 ft 已含库仑项 -sign(vt)*friction*abs(fn)
fw  = 0.5*(ft*vt + (fnv)*vn)                   ! 每节点记一半
dissipation(1,ni) = dissipation(1,ni) - fw
dissipation(1,nj) = dissipation(1,nj) - fw
```

然后在 `calculation.f90:2064` 并入固相热功率并推进温度：

```fortran
if(decomposition(1,nn) < 0.8) then
    heat(1,nn) = (heat(1,nn) + dissipation(1,nn))*dt + dissipation(2,nn)
    heat(2,nn) = heat(1,nn)/(EC(11,nodec(2,nn))*(1.0 - decomposition(1,nn)))
    heat(3,nn) = heat(3,nn) + heat(2,nn)          ! 固相温度
end if
```

要点：**法向只计入粘性项 `fnv`**（弹性项 `k·δ` 是可逆储能，不能算耗散）；
两节点各记一半；`dissipation(2)` 是断裂能 `Qbreak`（本实现未纳入，见 §6）。

## 3. 本实现

### 3.1 耗散记账字段

`DEMContact` 新增（单位 **J**，即"本步该对被耗散的能量"，而非功率）：

```cpp
double dissipation = 0.0;
```

放在 `DEMContact` 上而不是粒子上的原因：接触模型是**成对**量（`v_n`、`delta_t`、
库仑上限都只在这一对里可见），由模型就地写入最自然，也不需要在模型里持有求解器状态。

### 3.2 法向阻尼耗散

两个接触模型的 `CalcNormalForce`（`dt` 参数此前未被使用，现改为使用）：

```cpp
contact.dissipation = c_n * v_n_rel * v_n_rel * dt;      // 功率 c_n·v_n² ≥ 0
```

由动量守恒可推出该对被粘性耗散的功率恰为 `c_n·v_n_rel²`（恒非负）；
**Hertz 模型同理**（其 `δ ≤ 0` 的提前返回处把耗散清零）。

### 3.3 切向摩擦耗散

关键点：LSM 的离散体模块里切向力是**纯库仑滑动**（`ft = -sign(vt)·μ|fn|`），
所以 `ft·vt ≤ 0` 恒成立，`fw` 天然非负；而本工程的接触模型是
**增量切向弹簧 + 库仑截断**（Cundall–Strack），切向力在粘滞-滑移转换过程中
可能瞬时与滑移同向，直接照搬 `ft·vt` 会出现**负的耗散**（等于给粒子"制冷"）。

因此取**库仑截断瞬间切向弹簧储能的跌落**：

```cpp
if (Ft.norm() > Ft_max)
{
    const double stored_before = 0.5 * k_t * contact.delta_t.squaredNorm();
    Ft = Ft.normalized() * Ft_max;
    contact.delta_t = -Ft / k_t;                       // 滑动时重置弹簧位移
    const double stored_after  = 0.5 * k_t * contact.delta_t.squaredNorm();
    contact.dissipation += stored_before - stored_after;
}
```

它**严格非负**（截断条件 `k_t|δ_t| > μ|F_n|` 保证 `|δ_t|` 必然下降），而且可以证明
在持续滑动时等于库仑滑动功：设 `d_new = μ|F_n|/k_t`、`d_pre ≈ d_new + |v_t|·dt`，

```
ΔE = 0.5·k_t·(d_pre² − d_new²) ≈ k_t·d_new·|v_t|·dt = μ|F_n|·|v_t|·dt   ✓
```

即与 LSM 的 `ft·vt` 在滑动状态下的取值一致，但不会在粘滞段出现负值。
弹簧"充能"到库仑上限之前的能量属弹性储存，不计入耗散（这也意味着
下文的解析校验需要让充能段足够短）。

### 3.4 键合阻尼耗散

`DEMSolver::CalcBondForce`（阻尼力为 `-nd·v_n·n - td·v_t`）：

```cpp
const double damping_energy = (bond.normal_damping * velocity_n * velocity_n
    + bond.tangential_damping * velocity_t.squaredNorm()) * dt;
m_tmp_dissipation[bond.idx_a] += 0.5 * damping_energy;   // 两端各半
m_tmp_dissipation[bond.idx_b] += 0.5 * damping_energy;
```

`velocity_n` / `velocity_t` 是该键两端在**力臂处**的相对速度（含 `ω×arm` 贡献），
与阻尼力所用的完全一致。

### 3.5 汇总与分配

`DEMSolver` 新增成员缓冲 `dynamic_array<double> m_tmp_dissipation`（每粒子本步耗散能 J）：

| 来源 | 分配 |
|---|---|
| 粒子-粒子接触 | `0.5·c.dissipation` 给两端（同 LSM 的 `dissipation -= fw` 记两节点） |
| 粒子-墙面接触 | 全部给该粒子（墙不可动、无质量） |
| 键合阻尼 | 两端各半 |

在 `ZeroForce()` 中清零/扩容，与 `force`/`torque` 一起每步重置。
若 `contact_stiffen_ratio` 触发了法向排斥力放大，法向阻尼耗散同步乘同一放大系数
（放大发生在切向计算之前，此时 `dissipation` 只含法向项）。

### 3.6 温度解耦

`CalcThermalReaction()` 改为：

```cpp
const bool reaction = GetDEMParam()->GetReactionEnabled();
if (reaction) UpdateGasVoronoiMesh();
...
if (reaction) { /* 长度比/体积比/气固状态/气相压力汇总/内部换热/核心燃尽 */ }
// 键合导热：无条件
if (reaction) { /* 颗粒间表面燃烧 */ }
for (i : particles) {
    if (reaction) { p.reaction_rate = 0; ... }
    // ==== 以下对所有算例（含惰性）都执行 ====
    energy_delta[i] += p.heat_source * dt;                       // 外部体热源
    if (GetMechanicalHeatingEnabled()) energy_delta[i] += m_tmp_dissipation[i];
    if (reaction) { /* Arrhenius 基体热分解、气相内能 */ }
    /* ΔT = ΣE /(m·c_p·(1−φ))、上下限夹紧 */                      // 无条件
    if (reaction) { /* phase 判定、UpdateGasSolidState */ }
}
```

只把**会改写惰性算例输出**的反应量（Voronoi 网、`volume_ratio`、气固状态、
气相压力、燃烧、Arrhenius、phase）留在 `reaction` 分支内；
键合导热（含 `bond.heat_flow_a` 输出）、外部热源、机械耗散与温度更新无条件执行。
`reaction_enabled = false` 时反应相关分支整段跳过，因此既有惰性算例的
`gas_*`、`volume_ratio` 等字段输出与改动前一致。

## 4. 验证

### 4.1 法向阻尼（解析精确）——`tests/dem_damping_heating_test`

两个等质量球对心相撞（`±1 m/s`、`e = 0.5`、`reaction_enabled = false`）。
线性弹簧-阻尼在接触结束时弹性项完全卸载，故**总耗散 = 动能损失**
`m·u²(1−e²)`，两端各半：

```
kinetic energy       : 1.000000 J -> 0.250000 J (e=0.5)
dissipated total     : 0.750000 J  -> 0.375000 J per particle
temperature rise     : 3.750000e-03 K / 3.750000e-03 K  (解析 3.750000e-03 K, err 0.000%)
bounce velocities    : -0.500106, +0.500106 m/s  (解析 ∓0.500)
final separation     : 0.103661 m > 2r  (接触已解除，故全部损失都是耗散)
```

同时校验两端各得一半、末态对称、接触已解除。

### 4.2 切向摩擦（解析）——`tests/dem_friction_heating_test`

滑块（`r = 0.05`）骑在**远大于它的规定运动基底球**（`r = 5.0`，曲率近似平面）上。
重力调到使接触处于静力平衡（`Fn = m·g = k_n·δ`），于是全程法向力恒定、
法向相对速度为 0，唯一耗散就是库仑摩擦 `μ·Fn`：

```
slider velocity      : 4.978570 m/s  (解析 4.978240, err 0.007%)
friction deceleration: 21.7604 m/s^2 = mu*g
frictional work total: 1.085653e-01 J  -> 5.428267e-02 J 给每一端
temperature rise     : 滑块 5.500000e-04 K, 基底 5.500000e-04 K  (解析 5.428267e-04, err 1.32%)
```

> **为什么基底要用大球**：若两球半径相当，切向滑移会让接触点离开正下方、
> 中心距超过 `r_a+r_b` 而**间歇性丢失接触**，摩擦时断时续（实测只有理论值的 59%）。
> 这是几何问题而非实现问题 —— 用大基底把接触法向的倾角压到 ~1e-3 rad 即可。

> **为什么参数取 δ=2e-4、u0=5 m/s**：切向弹簧需先"充能"到库仑上限
> （`t_charge = 2μδ/u0`），充能段不产生耗散。选这组参数使 `t_charge` 仅为窗口的 1.6%。
> 1.32% 的残差主要来自温度输出只有 6 位有效数字（1.000550 K 的量化步长 1e-5）。

### 4.3 惰性算例的大规模解耦验证——`_demtest/inert_heating_check.py`

惰性外切巴西盘（10266 粒子、`reaction_enabled = false`），在 27 列输入里给试件粒子
加外部体热源 `heat_source = 1 W`，4000 步（`t = 8e-4 s`）：

| 配置 | 试件平均温度 | 说明 |
|---|---|---|
| `heat_source = 1 W`、`mechanical_heating = true` | 300.000 → **300.041 K** | 解析 `heat_source·t/(m·c_p) = 4.1354e-2 K`；实测 4.100e-2 K，**差 0.85%（6 位有效数字量化）** |
| `heat_source = 1 W`、`mechanical_heating = false` | 同上 | 开关只管机械耗散，不影响外部热源 |
| `heat_source = 0` | **300.000 → 300.000 K** | 无耗散源则不变（见下） |

> **注意**：该算例即使开着机械耗散生热也**不会因为压板加载而升温**，因为
> `friction_coeff = 0` 且位移控制压板的推压是准静态的（接触点法向相对速度趋于 0），
> 本身就没有可耗散的机械能。这是正确的物理结果，不是开关失效 ——
> "惰性算例会升温"这一能力由 §4.1 / §4.2 两个解析算例直接验证
> （两者都是 `reaction_enabled = false`）。

> 统计口径：试件粒子 `m·c_p = 1.93e-5 × 1000 = 0.0193 J/K`，而压板粒子
> `m·c_p = 1.0 × 1000 = 1000 J/K`，两者升温差 5 个量级；
> 用全局平均会把升幅稀释约 3%，必须按 group 统计。

### 4.4 回归

**16/16** 全部通过；精确值未变：`alpha=0.0956179`、`broken_bonds=240`、
`volume=8.999992000`、`fracture_energy=0.577350`、
HMX 的 `max_core_dalpha=6.002e-04` / `max_neighbor_dalpha=1.118e-03`。

## 5. 使用

```toml
[dem]
mechanical_heating = true     # 默认；关掉则温度只由键合导热与外部体热源驱动
# reaction_enabled = false    # 惰性算例现在也会升温，不再冻结
```

注意惰性算例要真正看到升温，还需要**存在耗散源**：
`friction_coeff > 0`（切向摩擦）或 `restitution_coeff < 1`（法向阻尼）或
配置了 `normal_viscosity` / `tangential_viscosity`（键合阻尼）。
`dem_damping_heating_test` 与 `dem_friction_heating_test` 都显式把初始温度设为
**1.0 K**，因为备份文件只有 6 位有效数字，`300 K` 上下无法分辨 1e-4 K 量级的温升。

## 6. 未纳入的部分（如需可再补）

1. **断裂释放能**：LSM 的 `dissipation(2)`（`Qbreak`）会并入热；本实现未纳入。
   它会直接改变反应算例的热量收支（裂纹尖端热点正是含能材料点火机制之一），
   需配合重新标定 `hmx_*` 与 `dem_reaction_test` 的期望值。
2. **切向阻尼**：接触模型本身没有 `c_t`（见 §4 的对照），因此可计入的切向耗散
   只有库仑截断那一项。
3. **墙面为非物理近似**：`DEMWall` 的材料参数是硬编码的（`young_modulus = 1e9`、
   `friction = 0.4`、`restitution = 0.9`），因此墙面摩擦生热用的是 `min(粒子 μ, 0.4)`。
4. `delta_t` 因法向转动而做切向平面投影时也会丢弃少量能量，未计入耗散（二阶小量）。
