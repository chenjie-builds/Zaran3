# DEM 初始弹簧连接网络（可选功能）——实现要点

> 用途：在仿真开始时（t=0）按几何邻近关系一次性建立粒子间的**弹簧连接**，
> 用于描述脆性材料的断裂、破坏等物理现象。
> 该机制与**逐步接触检测**相互平行、彼此独立，默认关闭，开启后不影响原有逻辑。

- 入口开关：`dem.spring_network_enabled`
- 代码位置：
  - `src/Zone/DEM/DEMSolver.cpp` → `DEMSolver::BuildSpringNetwork()`
  - `src/Zone/DEM/DEMSolver.cpp` → `DEMSolver::CalcBondForce()`（力的计算，复用）
  - `inc/Zone/DEM/DEMSolverParam.h` / `src/Zone/DEM/DEMSolverParam.cpp`（参数）
  - 数据存储：`inc/Basic/DEM/DEMBond.h` + `inc/Zone/DEM/DEMFieldData.h`

---

## 1. 初始化方式与数据存储

### 1.1 调用时机

在 `DEMSolver::InitField()` 中调用，位置经过刻意安排：

```
粒子材料默认值应用（含 k_n 推导所需的 E、ν）
        ↓
if (dem.spring_network_enabled) BuildSpringNetwork();   ← 仅在 t=0 执行一次
        ↓
键合统一后处理（由表面能换算断裂能、由粘度换算阻尼、合法性校验）
        ↓
UpdateGasVoronoiMesh()
```

放在"后处理"之前，是为了让自动生成的连接与 `bonds.csv` 提供的键合**共用同一套**
断裂能/阻尼/校验逻辑，保证两种来源物理一致。

### 1.2 建链判据与邻域搜索

- 判据：`dist(a,b) ≤ (r_a + r_b) · (1 + dem.spring_network_gap)`
  - `gap = 0` → 仅连接**恰好相切/重叠**的粒子；
  - `gap > 0` → 允许连接存在初始缝隙的近邻（脆性材料通常用小幅 gap）。
- 邻域搜索**复用接触检测的均匀网格**（`BuildUniformGrid` / `GridCellOf`，链表式单元表）：
  取单元边长 `2·r_max·(1+gap)`，则任一候选连接必落在相邻 3×3×3 单元内，
  复杂度 O(N)，只在 t=0 执行一次。
- 每对粒子只建一条连接（`j > i` 过滤）；与 `bonds.csv` 已有连接按 (idx_a,idx_b) 去重，
  因此**两种来源可同时使用而不重复**；`bond.id` 从已有最大 id 之后续编。

### 1.3 存储

复用 `DEMBond`，存于 `DEMFieldData::m_bonds`（`dynamic_array<DEMBond>`），
与接触对列表 `m_contacts` **完全独立**：

| 字段 | 自动生成时的取值 |
|---|---|
| `idx_a` / `idx_b` | 粒子**数组索引**（非 id） |
| `rest_length` | 初始中心距（→ 初始无应力） |
| `normal_stiffness` | 用户给定，或按与线性接触模型一致的 `k_n = 2·E*·R*` 自动推导 |
| `tangential_stiffness` | 用户给定，或 `0.5·k_n` |
| `normal_damping` / `tangential_damping` | 0，随后由 `normal_viscosity`/`tangential_viscosity` 统一换算 |
| `fracture_energy` | 由表面能换算；或由断裂应变换算 `0.5·k_n·(ε_f·L0)²` |
| `active` / `damage` / `delta_t` / `extension` / `force_a` / `elastic_energy` / `dissipated_fracture_energy` | 初始为默认（0 / true） |
| `source` | **1**（= 初始时刻自动生成的弹簧连接网络）；由 `bonds.csv` 读入的为 **0** |

连接状态会随每次输出写入 `backup/iter=*/bonds.dat`（含 `active`、`damage`、
`extension`、`fx/fy/fz`、`elastic_energy`、`fracture_energy`、`dissipated_fracture_energy`、
`source`），因此断裂模式可直接后处理。`source` 列用于区分两种来源，便于脚本筛选。

---

## 2. 每个时间步中与接触检测的协同顺序

`DEMSolver::Solve()` 的固定顺序：

```
ZeroForce()              // 力/力矩清零
CalcBondForce()          // ← 弹簧连接网络（自动生成 + bonds.csv），与接触无关
ContactDetection()       // ← 逐步接触检测（均匀网格 + 接触历史）
CalcContactForce()       // ← 接触力（法向/切向摩擦）
CalcWallForce()          // ← 墙面接触力
AdvanceContactHistory()  // 推进切向接触历史
CalcGravity()
CalcThermalReaction()
CalcGasPressureForce()
Integrate()              // 半隐式 Euler 积分
```

要点：

- 弹簧网络在**每个时间步**先于接触检测求值，但**不做任何邻域搜索**——连接在 t=0 已固定，
  每步只按 `idx_a/idx_b` 取端点状态算力，成本 O(连接数)。
- 弹簧网络与接触检测**互不读取对方数据**，仅通过 `pa.force / pb.force` 独立累加
  （`pa.force += force_a; pb.force -= force_a;`，作用力-反作用力）。
  因此两者可任意叠加，也可各自单独启用。
- 连接不会因为断裂而删除数组元素，故粒子索引始终稳定，`t=0` 建立的连接可长期有效。

---

## 3. 弹簧力与断裂条件的计算（`CalcBondForce`）

对每条 `active` 连接：

1. 几何：`Δx = x_b − x_a`，`length = |Δx|`，`n = Δx/length`，`extension = length − rest_length`
2. 轴向：`F_n = k_n·extension·n`（拉伸为正，压缩为负）
   - 压缩时可启用 LSM LEP/REP 高压硬化支路（`high_pressure_rnn/knn/exponent`）
3. 切向：`delta_t += v_t·dt` 后投影回切向平面，`F_t = −k_t·delta_t`
4. 阻尼：`−normal_damping·v_n` 与 `−tangential_damping·v_t`
5. `force_a = F_n + F_t`，`pa.force += force_a`，`pb.force −= force_a`

**断裂判据（两条，任一满足即断）**

| 准则 | 触发条件 | 说明 |
|---|---|---|
| 能量准则（优先） | `elastic_energy ≥ fracture_energy` | `fracture_energy` 由 `surface_energy`（含晶界能修正）换算；未设表面能时由 `spring_network_fracture_strain` 换算为 `0.5·k_n·(ε_f·L0)²` |
| 应变准则 | `maximum_tensile_strain > bond_break_strain` | 先按 `bond_peak_strain → bond_break_strain` 区间累积 `damage`（不可逆双线性损伤），超过上限即断 |

其中 `elastic_energy = 0.5·(k_n·extension² + k_t·|delta_t|²)`（反应算例再乘固相面积因子）。

另有两条**反应耦合**的失效路径（仅 `reaction_enabled = true` 时生效）：
任一端点燃尽（`reaction_progress ≥ burnout_threshold`）→ 连接失效。

---

## 4. 断裂后的状态更新

断裂当步（`CalcBondForce` 内）：

```cpp
bond.dissipated_fracture_energy += bond.elastic_energy; // 能量账本
bond.damage = 1.0;
bond.active = false;
bond.force_a.setZero();      // 不再传力
continue;                    // 本步不再向粒子累加该力
```

- **连接保留在列表中**（只标记 `active=false`），用于：
  - 输出断裂模式与能量耗散统计（`bonds.dat`）；
  - 保持粒子索引与连接拓扑稳定。
- 后续时间步遇到非 `active` 连接直接跳过，因此断裂后该粒子对的相互作用退化为
  "接触检测（若几何重叠）+ 其它体积力"，与弹簧网络无关。

---

## 5. 通过配置启用/禁用

### 5.1 开关与参数（`zaran.toml` 的 `[dem]` 段）

| 参数 | 默认值 | 含义 |
|---|---|---|
| `spring_network_enabled` | `false` | **总开关**。关闭时完全不建立任何连接，行为与改动前逐位一致 |
| `spring_network_gap` | `0.02` | 建链的相对间隙容差：`dist ≤ (r_a+r_b)·(1+gap)` |
| `spring_network_stiffness` | `0.0` | 法向刚度 (N/m)；`≤0` 表示按 `k_n = 2·E*·R*` 自动推导 |
| `spring_network_tangential_stiffness` | `0.0` | 切向刚度 (N/m)；`≤0` 表示取 `0.5·k_n` |
| `spring_network_fracture_strain` | `0.02` | 断裂应变阈值；`0` 表示回退到 `bond_break_strain`。**当 `surface_energy > 0` 时被能量准则取代** |

与之配套的既有参数（对自动连接同样生效）：
`surface_energy`、`grain_boundary_energy`、`bond_break_strain`、`bond_peak_strain`、
`normal_viscosity`、`tangential_viscosity`、`lattice_thickness`。

参数校验：`gap ≥ 0`、两个刚度 `≥ 0`、`fracture_strain ≥ 0`，否则
`DEMSolverParam::Init()` 抛 `ZaranError("DEM spring network parameters are invalid")`。

### 5.2 启用示例

```toml
[dem]
# 脆性材料：t=0 建立弹簧网络，拉伸 1% 即断裂
spring_network_enabled = true
spring_network_gap = 0.05
spring_network_fracture_strain = 0.01
# 若改用表面能控制断裂，则设置下面这项（优先生效）
# surface_energy = 1680.0
```

### 5.3 禁用

`spring_network_enabled = false`（默认值），或直接省略该键 —— 不生成任何连接，
`CalcBondForce` 仅在存在 `bonds.csv` 键合时才有工作量。

---

## 6. 验证

新增两个已提交算例（`_demtest/run_committed_tests.py` 纳入回归）：

| 算例 | 校验内容 |
|---|---|
| `tests/dem_spring_network_test` | 启用：t=0 自动生成 **1 条**连接且 `rest_length` = 初始中心距（0.11 m > 2r，证明该对**本不接触**）；拉伸阶段 `extension>0`、轴向力为回拉力（+x，实测 +267.97 N）、`elastic_energy = 0.5·k_n·extension²`；达到断裂应变后 `active=0`、`damage=1`、`dissipated_fracture_energy>0`、`fx=0`；断裂后线动量守恒、粒子以恒定速度飞离且相比初值已被弹簧减速 |
| `tests/dem_spring_network_off_test` | 禁用：完全相同的粒子与初速度下连接数为 **0**、粒子受力为 0、速度保持 ±1 m/s 不变，证明开关关闭时不影响原有行为 |

规模验证：在 10266 粒子 / 29635 键的巴西圆盘算例上开启（`gap = 0.10`、`fracture_strain = 0.05`）
运行 400 步，日志实测：

```
DEMSolver: spring network built at t=0 -> 280 connections
           (29635 pre-existing pairs, 29915 bonds total, gap=0.1, fracture_strain=0.05)
```

即 t=0 一次性建成、与 `bonds.csv` 已有连接自动去重，运行无发散、无 NaN。

全部已提交测试（12 项）通过；恢复系数扫描、接触力黄金值、守恒律等结果与开启该功能前**完全一致**（默认关闭）。

---

## 7. 与 `bonds.csv` 的关系

| | `bonds.csv` 给定的键合 | 自动生成的弹簧网络 |
|---|---|---|
| 来源 | 输入文件显式给出拓扑与刚度 | t=0 按几何邻近关系推导 |
| 开关 | `dem.bond_file` 是否存在 | `dem.spring_network_enabled` |
| 存储 | `DEMFieldData::m_bonds`（同一容器） | 同左 |
| 求解 | `CalcBondForce` | 同左 |

两者可**同时启用**：自动生成时会跳过已由文件给出的粒子对，避免重复连接。
输出中由 `source` 区分：`0` = 文件给定，`1` = 自动生成的弹簧连接网络。

---

## 8. VTP 可视化：如何一眼看出是否断裂

`.vtp` 中每条连接是一个**线单元**（`<Lines>`），并随附以下单元数据与点数据。
键合状态数组由 `output.bond_details` 控制（**默认 `true`**；置 `false` 则退回历史输出，
便于超大规模键数时压缩体积）。

### 8.1 单元数据（CellData，每条连接一个值）

> 注意：VTK 的 CellData 需覆盖全部单元，因此数组前 `N` 个值是粒子（Verts）的占位 `0`，
> 第 `N+1 … N+B` 个值才对应第 1…B 条连接。

| 数组名 | 类型 | 含义 |
|---|---|---|
| **`bond_active`** | Int32 | **`1` = 完好，`0` = 已断裂** —— 判断断裂的核心字段 |
| **`bond_source`** | Int32 | `0` = `bonds.csv` 给定，`1` = t=0 自动生成的弹簧连接网络 |
| `bond_stiffness` | Float64 | 法向刚度 `k_n` (N/m) |
| `bond_rest_length` | Float64 | 无应力长度 `L0` (m) |
| `bond_strain` | Float64 | 当前轴向应变 `extension/L0`（正 = 拉伸，负 = 压缩） |
| `bond_max_strain` | Float64 | 历史最大拉伸应变，与 `spring_network_fracture_strain` / `bond_break_strain` 直接可比 |
| `bond_force_magnitude` | Float64 | 键合力大小 `|F_a|` (N)，断裂后恒为 `0` |
| `bond_energy_loss` | Float64 | 断裂耗散的能量 (J)，`> 0` 即真实发生过断裂 |

原有的 `bond_id`、`bond_extension`、`bond_force_a`(3 分量)、`bond_heat_flow_a`、
`bond_damage`、`bond_elastic_energy`、`bond_fracture_energy` 保持不变。

### 8.2 点数据（PointData，每个粒子一个值）

| 数组名 | 类型 | 含义 |
|---|---|---|
| `bonds_incident` | Int32 | 与该粒子相连的连接条数（配位数） |
| `bonds_broken` | Int32 | 其中已断裂的条数 |
| `bonds_broken_ratio` | Float64 | 已断裂比例 `bonds_broken / bonds_incident`，**用它上色可让裂纹带/断裂面直接凸显** |

### 8.3 ParaView 操作建议

1. 打开 `.vtp` → 需要同时显示点与线时，用 **Extract Cells By Type / 直接显示**；
   球体可用 `Sphere` glyph（半径取 `radius`）叠加线框查看。
2. **只看断裂的连接**：对线单元做 `Filters → Threshold`，选 `bond_active`，
   范围 `0 … 0`（正是断裂的那些）；调 `1 … 1` 则只看完好连接。
3. **只看弹簧连接网络**：`Threshold` 选 `bond_source`，范围 `1 … 1`
   （与 `bonds.csv` 给定的晶格键合区分开）。
4. **着色看受力/伸长**：对线单元用 `bond_force_magnitude` 或 `bond_strain` 上色；
   已断裂连接会自然落在大片零值区（已归零），与完好段形成对比。
5. **看断裂程度/逼近断裂**：用 `bond_elastic_energy / bond_fracture_energy` 的
   计算器表达式（或直接比较 `bond_max_strain` 与配置的断裂应变）。
6. **整体损伤云图**：切到点数据用 `bonds_broken_ratio` 上色，
   断裂区会连成清晰的高值带；也可用 `bonds_broken > 0` 做 `Threshold` 提取受损粒子。

### 8.4 验证

`_demtest/vtp_bond_check.py` 与 `tests/dem_spring_network_test/verify.py` 会解析 `.vtp` 并校验：

- `NumberOfLines` == 连接条数；每 个数组长度 == `(N+B) × 分量数`（点数据为 `N`）；
- t=0：`bond_source = 1`、`bond_active = 1`、`bonds_incident` 合计正确；
- 断裂后：`bond_active = 0`、`bond_force_magnitude = 0`、`bond_energy_loss > 0`、
  两个端点的 `bonds_broken` 各计 1、`bonds_broken_ratio = 1`；
- `bond_extension` / `bond_strain` / `bond_force_magnitude` 与 `bonds.dat` 数值一致；
- `output.bond_details = false` 时新数组全部消失，**旧数组逐值不变**（输出可比）。

实测开销（`dem_brazilian_disk_10k`：10266 粒子 / 29635 键，400 步 / 2 帧）：

| `output.bond_details` | VTP 总体积 | 每帧 | 端到端耗时 |
|---|---|---|---|
| `true`（默认） | 27.3 MB | 13.6 MB | 10.8 s |
| `false` | 17.7 MB | 8.9 MB | 10.8 s |

即体积 **+54%**（每帧 +4.8 MB），耗时无可测差异（写盘走整块缓冲，仅多格式化若干列）。
键数极大（≫10⁶）且只关心几何时，可置 `false` 压体积。
