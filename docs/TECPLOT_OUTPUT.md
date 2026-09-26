# Zaran3 的 Tecplot 输出（粒子 / 键合 / 流场）

> 目标：**只用 Tecplot 也能看完 DEM 与流场的一切**。
> 之前 DEM 侧只写 VTK PolyData（`.vtp`），Tecplot 读不了；现在统一为 Tecplot ASCII。

---

## 1. 文件清单

一次输出（每个输出步）最多产生三个文件，全部是 Tecplot ASCII：

| 文件 | 内容 | zone 类型 | 打包 |
|---|---|---|---|
| `result/particles_<iter>.dat` | DEM 粒子（点云）| `FEPOINT` | `DATAPACKING=POINT`（一行一粒子）|
| `result/bonds_<iter>.dat` | 永久键合（线元）| `FELINESEG` | `DATAPACKING=BLOCK` |
| `result/<iter>.dat` | 流场（结构网格）| `ORDERED`（I/J/K）| `DATAPACKING=POINT` |

三者都带 `STRANDID=1` 与 `SOLUTIONTIME=<物理时间>`，因此**可以一起做时间动画**。

---

## 2. 为什么分成三个文件（而不是一个）

Tecplot ASCII 的 **`VARIABLES` 是整个文件共享的**，不是每个 zone 各自一套。
粒子是点云（惰性 25 列、含能 43 列），键合是线元（3 个节点量 + 20 个单元中心量 = 23 列），
流场是网格（8 或 9 列，其中 X/Y/Z 的含义虽然相同但数据条数完全不同）。

硬塞进一个文件只有两条路，代价都不可接受：

- 用**并集变量表**去凑 ⇒ 每个 zone 都要写一大堆空列（粒子文件凭空多 20 列零、键合文件多 18 列），
  体积翻倍且表格难读；
- 让三者共用一套变量名 ⇒ 语义直接错位（键合量不是节点量，网格量不是粒子量）。

所以按"数据性质"分开写。这在 Tecplot 里不是问题：多文件可以同时载入，
Tecplot 360 支持同一帧显示多个数据集（见 §5.5）。

> 如果你只关心粒子（大多数 DEM 场景），载入 `particles_*` 一个文件就够了。

---

## 3. 变量表

### 3.1 粒子文件 `particles_<iter>.dat`

前 3 列恒为坐标，随后是**核心列（恒输出，共 25 列）**：

| # | 变量 | 说明 |
|---|---|---|
| 1–3 | `X`,`Y`,`Z` | 质心坐标 (m) |
| 4–5 | `id`,`group` | 全局编号、分组 |
| 6–7 | `radius`,`mass` | 半径 (m)、质量 (kg) |
| 8–9 | `active`,`kinematic` | 是否参与计算 / 是否按给定速度运动 |
| 10–12 | `velocity_x/y/z` | 平动速度 (m/s) |
| 13–15 | `omega_x/y/z` | 角速度 (rad/s) |
| 16–18 | `force_x/y/z` | 合力 (N) |
| 19–21 | `torque_x/y/z` | 合力矩 (N·m) |
| 22 | `temperature` | 固相温度 (K) |
| 23–25 | `bonds_incident`,`bonds_broken`,`bonds_broken_ratio` | 相邻键数 / 已断键数 / 断裂比例 |

**含能算例**（任一粒子 `energetic=true`，或 `output.particle_fields="full"`）在 `temperature`
之后**再追加 18 列**：`energetic`,`phase`,`rotation_x/y/z`,`reaction_progress`,`reaction_rate`,
`gas_temperature`,`gas_pressure`,`volume_ratio`,`solid_volume`,`gas_volume`,`solid_core_radius`,
`gas_internal_energy`,`internal_heat_transfer`,`body_reaction_increment`,`core_burn_increment`,
`neighbor_burn_increment`（共 43 列）。

> `bonds_broken_ratio` 是看损伤最顺手的着色变量：裂纹带会直接凸显。

### 3.2 键合文件 `bonds_<iter>.dat`

23 列 = **3 个节点变量**（`X`,`Y`,`Z`）+ **20 个单元中心变量**：

`bond_id`, `active`, `source`, `normal_stiffness`, `tangential_stiffness`, `rest_length`,
`extension`, `strain`, `max_strain`, `max_compression`, `strength_scale`,
`force_x`, `force_y`, `force_z`, `force_magnitude`, `damage`,
`elastic_energy`, `fracture_energy`, `energy_loss`, `heat_flow_a`

关键列：

- `max_compression`：历史最大压缩应变。与"压缩断裂阈值 = `bond_compression_strength_ratio` ×
  抗拉阈值 × `strength_scale`"直接可比 ⇒ **判断破坏是压碎还是劈裂**就看它。
- `strength_scale`：逐键 Weibull 强度折减（1 = 均质）。
- `active`：1 = 完好，0 = 已断。用 Threshold 只留 `active = 1` 即可隐藏断键。
- `energy_loss`：断裂时耗散的能量（>0 即真的断过）。

### 3.3 流场文件 `<iter>.dat`

`X`,`Y`,`Z`,`Density`,`Velocity_x`,`Velocity_y`,`Velocity_z`,`Pressure`
（两相算例再追加 `Volume_fraction`）。

---

## 4. 配置键（全部写在 `[output]` 段内）

```toml
[output]
result_folder   = "result"

# ── DEM 粒子/键合 ───────────────────────────────────────────────
particle_format = "tecplot"   # tecplot(默认) | vtp | both | none
particle_fields = "auto"      # auto(默认) | core | full   —— 是否附加含能材料的 18 列
bond_details    = true        # Tecplot 模式下=false 则不写 bonds_*.dat
                              # ParaView 模式下=false 只裁掉部分数组

# ── 流场 ───────────────────────────────────────────────────────
tecplot_format  = "ascii"     # ascii(默认，.dat) | binary(.plt)

# ── 通用 ───────────────────────────────────────────────────────
tecplot_precision = 12        # 有效数字位数，6..17，默认 12
```

开关实测（`tests/dem_bond_test` / `tests/euler_uniform_sod`）：

| 配置 | 产生的文件 |
|---|---|
| `particle_format = "tecplot"` | `particles_*.dat` + `bonds_*.dat` |
| `particle_format = "both"` | 上面两个 + `particles_*.vtp` |
| `particle_format = "vtp"` | 只有 `particles_*.vtp` |
| `particle_format = "none"` | 不写任何可视化文件（只留 `backup/` 的 CSV）|
| `tecplot_format = "ascii"` | `result/<iter>.dat` |
| `tecplot_format = "binary"` | `result/<iter>.plt` |

兼容性：旧的 `output.tecplot_ascii`（bool）仍被识别，但**仅在未写新键时生效**
（新键优先），并在日志里提示改用新键。

算例生成器也已跟进：
`make_euler_case.py --tecplot-format ascii|binary`（默认 ascii）；
旧的 `--tecplot-ascii` 保留为别名，新增 `--tecplot-binary`。
生成的算例还带 `tecplot_precision = 17`，以便与精确解做逐位比对。

⚠️ 注意 TOML 的段作用域：键必须落在 `[output]` 段内。若把 `particle_format` 追加到文件末尾，
它会属于最后一个段（例如 `[dem.box]`），程序读不到、静默使用默认值。
另外同一键不能写两次（TOML 会直接报错退出）。

---

## 5. Tecplot 里的用法

### 5.1 载入
`File → Load Data`，选一个或多个 `.dat`。ASCII 文件无需 `preplot` 转换。

### 5.2 时间动画
所有 zone 都写了 `SOLUTIONTIME`，直接拖 `Solution Time` 滑块或用
`Plot → Animate → Time`。同一帧号的文件（`particles_5000.dat` 与 `bonds_5000.dat`）
时间一致，可以同步。

### 5.3 粒子（点云）
`particles_*.dat` 是 `FEPOINT` zone，用 **Scatter** 图（3D Cartesian + Symbols，关闭 Lines）。
建议着色变量：`bonds_broken_ratio`、`velocity`（用 `|v|` 表达式 `sqrt(vx**2+vy**2+vz**2)`）、
`temperature`。粒子尺寸可按 `radius` 缩放符号大小。

### 5.4 键合（线元）
`bonds_*.dat` 是 `FELINESEG` zone ⇒ 3D Cartesian 里显示 **Lines**（关掉 Symbols）。
- 着色 `strain` / `extension` 看法向受力；
- 着色 `max_compression` 看哪里在压碎；
- `Data → Extract → Threshold`，取 `active = 1` 只显示未断的键，取 `active = 0` 只看断裂面；
- 键太密时把 `Lines` 的线宽调细，或先 Threshold 出高应变段。

### 5.5 多数据集同帧
粒子/键合/流场的 `VARIABLES` 列表不同，Tecplot 会把它们当成不同**数据集**
（不是同一个数据集里的不同 zone）。Tecplot 360 支持在**同一帧**里显示多个数据集
（侧边栏的 Data Set 选择器 / `Frame → Data Set`）；如需严格单数据集，
就分别用 `Data Set` 切换，或只载入当前需要的那一个文件。

---

## 6. 实现上必须守住的四个格式约束

1. **FE zone 的 BLOCK 数据块顺序**：Tecplot 规定"带单元中心量的 zone 必须用 `BLOCK` 打包"，
   且数据顺序为 **① 全部节点变量（按 `VARIABLES` 顺序）→ ② 全部单元中心变量 → ③ 连接表**。
   `VARLOCATION=([4-23]=CELLCENTERED)` 声明第 4–23 列是单元中心量。
   顺序写错 Tecplot **不会报错**，只会把 x 坐标当 y 坐标画出来 —— 所以必须逐值比对验证。
2. **非有限值必须替换**：Tecplot 解析不了 `nan`/`inf`，写出一个就会导致**整个文件打不开**。
   实现里把非有限值写成 0 并计数告警。
3. **单行 ≤ 32000 字节**：Tecplot ASCII 的硬限制，超长会静默截断。实现里按固定列数换行。
4. **`Rep*Num` 重复记法**：`20000*1` 表示 20000 个 1。键合的 `active`/`source` 等常量块
   用它可以把文件缩小一个量级；但**粒子文件不用**（一行一粒子，保证脚本能直接解析）。

另外：`VARIABLES` / `ZONE` 参数 / `SOLUTIONTIME` 都写在**单独一行内**，
不做跨行续写（续行合法但容易踩解析歧义）。

---

## 7. 校验（不需要装 Tecplot）

```bash
python benchmarks/tecplot_ascii/check_tecplot_dat.py --case tests/dem_brazilian_disk_test
python benchmarks/tecplot_ascii/check_tecplot_dat.py --case tests/euler_uniform_sod --kind flow
```

做六项检查：

| 检查 | 内容 |
|---|---|
| C1 结构 | `VARIABLES` 唯一非空、zone 头参数齐备、N/E 与 zone 类型自洽 |
| C2 条数 | 节点块 = N、单元块 = E、连接表 = 2E（并核对 `VARLOCATION` 声明的列数）|
| C3 索引 | 连接表节点号全部落在 `[1, N]` |
| C4 行长 | 每行 ≤ 32000 字节 |
| C5 有限 | 全部数值有限 |
| C6 数值 | 与 `backup/iter=<N>/*.dat`(CSV) **逐值比对**（流场则核对结构化次序）|

**C6 是关键**：只有它能把"BLOCK 块顺序写错"这类 Tecplot 照单全收的错误查出来。
注意 BLOCK 是**变量优先**布局（同一变量的所有单元值连续），所以第 k 条键合的第 i 个单元量是
`cell[i*E + k]`，不是 `cell[k*20 + i]` —— **E = 1 时两者恰好等价**，
所以验证必须用多键合算例（这也是本脚本第一次跑单键合算例"假通过"的教训）。

---

## 8. 实测

**规模**（`tests/dem_brazilian_disk_10k`：10266 粒子 + 29635 键合，5 帧）

| 文件 | 大小 |
|---|---|
| `particles_<iter>.dat`（25 列）| ≈ 2.2 MB |
| `bonds_<iter>.dat`（23 列，含 `active`/`source` 常量块压缩）| ≈ 4.9 MB |

单帧写出 < 1 s，相对 10000 个 DEM 步（约 80 s）可忽略。

**校验结果**

| 算例 | 粒子 | 键合 | 判定 |
|---|---|---|---|
| `dem_bond_test`（2 粒子 / 1 键）| PASS | PASS | PASS |
| `dem_brazilian_disk_test`（139 / 222）| PASS | PASS | PASS |
| `dem_weibull_strength_test`（400 / 760）| PASS | PASS | PASS |
| `euler_uniform_sod`（200×3×3 网格）| — | — | PASS |

数值比对的最大归一化偏差在 1e-6 量级，与 `backup` CSV 只有 6 位有效数字相符。

---

## 9. 改动文件

| 文件 | 改动 |
|---|---|
| `inc/Basic/PostProcess/Visual.h` | 新增 `WriteParticleTecplotASCII` 两个重载 |
| `src/Basic/PostProcess/Visual.cpp` | 实现 DEM Tecplot ASCII 写出器；流场 ASCII 的 ZONE 参数收进单行 |
| `src/Main/DEMFieldSimulation.cpp` | 按 `output.particle_format` 分派，传入 `dem.current_time` |
| `src/Main/FieldSimulation.cpp` | 流场按 `output.tecplot_format` 分派（默认 ascii） |
| `benchmarks/euler_uniform/make_euler_case.py` | `--tecplot-format ascii\|binary`（默认 ascii），生成 `tecplot_precision = 17` |
| `benchmarks/tecplot_ascii/check_tecplot_dat.py` | 新增：六项结构/数值校验 |
| `tests/dem_spring_network_test/zaran.toml` | 显式 `particle_format = "both"`（其 verify.py 校验 VTP）|

**19/19 回归通过**；`dem_bond_test`、`dem_brazilian_disk_test`、`dem_weibull_strength_test`、
`dem_spring_network_test`、`euler_uniform_sod`、两相 sod（9 列）的 `.dat` 校验全部 PASS。
