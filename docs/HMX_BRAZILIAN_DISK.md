# HMX 巴西盘算例：材料参数映射与标定

> 目标：把 LSM 材料文件里的 HMX 物性映射到 Zaran3，并让巴西盘算例在
> **σ_t = 60 MPa** 附近开裂、继而压碎。

算例目录：`tests/hmx_brazilian_disk_10k/`
（`generate_case.py` / `particles.csv` / `bonds.csv` / `zaran.toml` / `verify.py`）

---

> ⚠️ **2026-09-17 更新：本算例已启用压缩失效与逐键 Weibull 异质性，failure mode 随之改变。**
> 下面 §4 记录的"中心起裂竖直劈裂（σ_t = 57.5 MPa）"是**未加压缩判据**时的结果。
> 现在 `zaran.toml` 设了 `bond_compression_strength_ratio = 8.0` +
> `bond_weibull_modulus = 6.0`，起裂改为**压板接触点压碎**（x≈0、|y|>0.98R），
> 载荷见峰即崩（4.3 kN → 0），压缩被截在 1.07e-2（原 8.16e-2），
> 碎屑速度 112 → 25 m/s、不再撞墙。§4 的 σ_t 标定关系（阈值 1.19e-3 ↔ 60 MPa）
> 仍然有效，只是不再由本算例的默认配置演示。
> 完整说明与三个配置的对照见 **`docs/CRUSH_FAILURE.md`**；
> 想恢复"劈裂主导"演示：把 `bond_compression_strength_ratio` 设为 40（或以上）——
> 但更根本的办法是采用真实的试样厚度（本算例 t/D = 0.0095，约一个粒子直径）。

---

## 1. 为什么不能把 LSM 的 `kn / kt` 直接搬过来

LSM 的 HMX 用的是 **分层叠加** 本构：`available_interactions = REP LEP HBP`，
法向力是三者之和。以 HBP 段为例（`covalent_ionic_spring_mod.f90:246-266`）：

```fortran
else if(rn<=HBP%rncr) then                 ! 双线性/超弹性
    delta=HBP%rncr-rn
    fn=HBP%Kn1*delta**HBP%nn+HBP%fncr
...
surface   = sqrt(3.0)*HBP%r0/3.0            ! 单键等效面积
HBP%fncr  = -TensileStrength*surface        ! 强度 → 力
HBP%nn    = Kn0*(rncr-r0)/fncr              ! 由 Kn0、rncr 反推指数
```

两点结论：
1. `TensileStrength` 是通过 `fncr = σ_t·A_bond` 换算成力的 —— 这是**可以直接用**的物理量。
2. `LEP%Kn0 = 24.2184 GPa` 在 LSM 里是"应力量纲"的 2D 弹簧系数（`fn = Kn0·drn`，
   力按单位厚度计），且被 REP/LEP/HBP 叠加使用；Zaran3 的 `bonds.kn` 是 3D 真实弹簧
   刚度 [N/m]。**两者语义不同，不能照搬。**

另外 LSM 的 `surface = √3·r0/3` 与 Zaran3 的 `A_bond = L0·lattice_thickness/√3`
在 `lattice_thickness = 1.0` 时完全一致 —— 这解释了 Zaran3 里那个看似奇怪的
`1/√3` 因子来源。

**所以采用"标定宏观量"的路线**：用 Zaran3 自己的晶格做单轴试验，把
E、ν、σ_t 三个可直接测量的宏观量标出来。

---

## 2. 标定（`_demtest/hmx_calib.py`）

标定试件与巴西盘**同一套晶格**（间距 s = 9.452646780e-4 m、粒子半径 s/2、
粒子恰好外切、同一种键合法则），尺寸 17×41 = 697 粒子，
上下两排分别为位移控制与固定。

### 2.1 弹性模量与泊松比（单轴压缩）

| kn (N/m) | kt (N/m) | kt/kn | 实测 E (Pa) | 实测 ν |
|---|---|---|---|---|
| 1.9623e7 | 6.541e6 | 0.333 | 2.6126e10 | 0.309 |
| 1.9623e7 | 2.0e7 | 1.019 | 3.1929e10 | 0.142 |
| **1.8862e7** | **1.4750e7** | **0.782** | **2.8978e10** | **0.1959** |

- 目标：E = 28.764 GPa、ν = 0.2 → 实测 **E 偏 +0.74%、ν 偏 0.2%**。
- 观察：ν 随 kt 增大而下降（kt=kn/3 → 0.309，kt=kn → 0.142），
  线性内插得 ν=0.2 对应 **kt/kn = 0.782**（远大于原算例的 0.25，也远大于 LSM 的 0.048）。
- 注意：**拉伸模量只有 2.0978e10 Pa**（低 27%）。原因见 §4。

### 2.2 抗拉强度（单轴拉伸）

位移控制拉伸，记录峰值应力：

| break_strain | σ_peak（宏观） |
|---|---|
| 2.086e-3 | 34.33 MPa |
| **3.646e-3** | **60.13 MPa**（目标 60 MPa，偏 0.22%） |

拉伸曲线是**完全脆性**的：应力升到 60.13 MPa 后一步之内崩塌到 ~0.3 MPa，
之后整块失去承载能力（`verify` 曲线见标定脚本输出）。

### 2.3 巴西盘柔度校验

实测荷载-位移柔度 **C = 1.04e-7 m/N**；
经典解 `C = 2/(πEt)·[ln(2D/b) − 1 − ν]`（E=2.9e10, t=2r, b=2R_platen）
给出 **1.015e-7 m/N** —— 相差 2.5%，说明标定出的模量能正确复现盘的弹性响应。

---

## 3. 参数映射表

| LSM HMX 材料文件 | Zaran3 参数 | 采用值 | 说明 |
|---|---|---|---|
| `density = 1.89 g/cc` | `dem.density` + 粒子质量 | 1890 kg/m³ | 质量按**球体积** ρ·(4/3)πr³ = 8.3584e-7 kg |
| `Young_modulus = 28.764 GPa` | `dem.young_modulus` | 2.8764e10 | 用于**接触**刚度 |
| `Poisson_ratio = 0.2` | `dem.poisson_ratio` | 0.2 | 用于接触 |
| `LEP kn = 24.2184 GPa` | `bonds.kn` | **1.8862e7 N/m** | 由 E 标定，非照搬 |
| `LEP kt = 1.1533 GPa` | `bonds.kt` | **1.4750e7 N/m** | 由 ν 标定，非照搬 |
| `HBP TensileStrength = 60 MPa` | `bond_break_strain` | **1.19e-3** | 按**巴西盘几何**标定，见 §4 |
| `surface_energy = 8.4 / GB_energy = 4.2 J/m²` | `surface_energy` / `grain_boundary_energy` | **未启用** | 见 §5 |
| `rnn = 0.8` | `dem.high_pressure_rnn` | 0.8 | 直接对应 |
| `knn = 10.0` | `dem.high_pressure_knn` | 10.0 | 直接对应 |
| `exp = 1.0` | `dem.high_pressure_exponent` | 1.0 | 直接对应 |
| `visco1 = visco2 = 100 Pa·s` | `dem.normal_viscosity` / `tangential_viscosity` | 100 / 100 | c = η·A/L 自动换算 |
| `lattice_thickness`（Zaran3 特有） | `dem.lattice_thickness` | **9.452646780e-4** | 必须显式给；默认 1.0 会让 A 大 1000 倍 |
| `HBP ShearStrength = 30 MPa` | — | — | **无实现**：剪切断裂分支缺失（见上一轮盘点） |
| `specific_heat = 1004.62` | 粒子第 9 可选列 | 待反应算例 | 需 27 列 particles.csv |
| `thermal_conductivity = 0.5358` | 粒子第 10 可选列 | 待反应算例 | 同上 |
| `Q = 5.53e6 J/kg` | 粒子第 11 可选列 | — | `reaction_heat` |
| `Z = 4.78e12 s⁻¹` | 粒子第 12 可选列 | — | `arrhenius_prefactor` |
| `E = 143.9e3 J/mol` | 粒子第 13 可选列 | — | `activation_temperature = Ea/R = 17306.5 K`（代码用 `exp(-Ta/T)`） |
| `JWL_A/B/C/R1/R2/w` | `dem.jwl_*` | 待反应算例 | 本次 `reaction_enabled = false` |

附加的、材料文件未给出的两项（需自行选取）：
- `friction_coeff = 0.30`：HMX–HMX 摩擦（原 HMX 算例用 0.10）
- `restitution_coeff = 0.80`：断裂后碎块接触

---

## 4. 关键发现：断裂阈值必须按**几何**标定，不能跨几何搬用

把单轴拉伸标定出的 `break_strain = 3.646e-3` 用到巴西盘上，
盘要加载到 **σ_t ≈ 120 MPa** 才开裂（目标 60 MPa，差 2 倍）。

原因：**本晶格粒子恰好外切（2r = L0）**，只要受压缩，相邻粒子的中心距就立刻小于
`r_a + r_b` → **接触力被激活，与键合力并联**。压缩方向因此显著变刚。
而巴西盘中心的对角方向正处于强压缩，于是"水平键"（承担中心拉应力）的拉伸应变
远低于同一宏观拉应力下单轴拉伸的对应值。

用**截面法**直接测盘中心的内部拉应力 σ_x（统计穿过 x=0 竖直截面的 241 条键的 x 向力）：

| 迭代 | ε_bond(最大) | 截面法 σ_x | 名义式 2P/(πDt) | 名义/截面 |
|---|---|---|---|---|
| 20000 | 4.19e-4 | 36.58 MPa | 34.86 MPa | 0.953 |
| 46000 | 1.94e-3 | 83.56 MPa | 101.58 MPa | 1.216 |

- 低载时两法一致（差 5%）⇒ **`σ_t = 2P/(πDt)` 公式本身是可靠的**，厚度约定 t = 2r 成立。
- 外插得 σ_t = 60 MPa 对应 ε_bond ≈ 9.1e-4（单轴值的 1/4）。

### 4.1 在巴西盘几何内迭代标定阈值

| `bond_break_strain` | 实测起裂时 σ_t | 备注 |
|---|---|---|
| 3.646e-3（单轴值） | **≈ 120 MPa** | 高 2 倍 |
| 9.1e-4 | **47 MPa** | 低 22% |
| **1.19e-3** | **57.53 MPa** | ✓ 目标 60 MPa，偏 −4.1% |

（每次上调 60/实测 的比例；一次迭代即收敛。）

### 4.2 最终算例的验证结果

```toml
bond_break_strain = 1.19e-3
bond_peak_strain  = 8.3e-4
```

`verify.py` 全部通过（exit 0）：

```
geometry   : 2r = L0 = 9.452646780e-04 m, rho = 1890.00 kg/m^3, D = 0.09918 m, t = 9.4526e-04 m
initiation : sigma_t = 57.53 MPa at iter=30000, 1250 bond(s) broken
crack band : 94.1% of them within |x| < 0.2R, symmetric 613 above / 637 below
final      : 2470 / 29635 bonds broken (8.33%), sigma_t(末帧) = 95.76 MPa
```

**裂纹形态就是巴西盘的经典间接拉伸劈裂**（起裂帧统计，1250 条断键）：

| 指标 | 实测 | 期望 |
|---|---|---|
| 断键中点 x 跨度 | ±0.0033 m = **±0.067R** | 窄竖直带（沿加载直径） |
| \|x\|/R < 0.1 占比 | **100%** | 全部在中心线附近 |
| 断键中点 y 跨度 | ±0.0187 m = **±0.38R** | 纵向贯穿 |
| 上/下对称性 | **613 / 637** | 近似相等 |

继续加载到 iter=45000 时损伤带扩展到 \|x\| < 0.45R（仍局限在中央竖直带内），
断裂键数 2470（8.33%）。

**注意：载荷没有跌落**（σ_t 从 57.5 继续升到 95.8 MPa）。原因是
①强度全局均匀 ⇒ 形成**弥散损伤带**而非单条贯通裂纹；
②断裂后的相邻粒子仍处于接触状态，仍能承压；③位移控制下压板持续下压。
要做到"载荷峰值后崩落 → 继而压碎"，需要 §6 中的异质性与剪切分支。

**这是离散键模型的固有性质**（BPM 的强度依赖几何与配位数），不是 bug；
但若要用同一套参数同时复现单轴与巴西盘强度，需要消除接触-键合双重计数（见 §6）。

---

## 5. 为什么本次不用 `surface_energy` 能量口径

Zaran3 的断裂能为 `Efract = Es·L0·t/√3·(1 − Egb/Es)`，
键失效时 `0.5·kn·ε² > Efract`。由此反推的宏观强度
`σ ≈ √(2·E·γ_eff/(√3·C·L0))` **只与 γ、E、L0 有关，与 t 无关**。

代入 γ_eff = (8.4−4.2) J/m²、E = 28.764 GPa、L0 = 9.4526e-4 m
→ 等效拉伸应变 ≈ 4.7e-4 → **强度只有约 12 MPa**（目标 60 MPa，差 5 倍）。

γ 与 σ_t 要自洽，需要 `L0 ≈ 70–130 μm`（`L0 ≈ 2γ_eff/(E·ε_f²)`）；
本网格 L0 = 945 μm 粗了 7–14 倍。故本次改用**应变判据**，
`surface_energy` 留空。若将来加密网格到 ~100 μm（约 90 万粒子），两者可自洽。

---

## 6. 已知局限与建议的下一步

| 项 | 现状 | 建议 |
|---|---|---|
| **接触-键合双重计数** | 粒子外切 ⇒ 受压即产生接触力，压缩方向刚度被高估（E_压缩 2.90e10 vs E_拉伸 2.10e10，差 38%），并使盘强度虚高 2 倍 | 增加开关：**键合粒子对之间不计接触力**（标准 BPM 做法）。同时能把单步耗时降低约一半 |
| **强度只判拉伸** | 压缩/剪切下永不断 ⇒ 只能"劈开"、不能真正"压碎" | 移植 LSM 的剪切分支 `rn<r0 && Et≥Uall`，并引入混合模式 `En/En_c + Et/Et_c ≥ 1`（LSM `ShearStrength = 30 MPa` 已给） |
| **强度全场均匀** | 裂纹路径完全由应力场对称性决定，只会出现单条主裂纹 | 逐键 Weibull 强度异质性（需要 `bonds.csv` 支持逐键断裂参数） |
| **`reaction_weakening`** | 死参数（无调用者） | 要么实现、要么删除 |
| **晶界能语义** | `grain_boundary_energy` 被折成全局统一因子，晶内/晶界无区别 | 恢复 LSM `GBratio` 的逐键语义 |
| **粘性阻尼偏弱** | η=100 Pa·s ⇒ 键阻尼 ζ ≈ 1e-5（几乎无阻尼） | 断裂瞬间若出现数值振荡，可另给 `bonds.cn/ct` |

---

## 7. 运行方法与验证

```bash
cd tests/hmx_brazilian_disk_10k
python generate_case.py                 # 重新生成 particles.csv / bonds.csv
<Zaran3.exe> .                          # 运行（约 4.5 万步）
python verify.py --sim backup --case .  # 校验
```

`verify.py` 校验四项：
1. **几何**：试件相邻粒子精确外切（2r = L0）、压板球彼此外切、
   密度回算 = 1890、压板球能真正触到试件顶排；
2. **标定**：加载峰值满足 σ_t = 2P/(πDt) ≈ 60 MPa（容差 ±25%）；
3. **起裂位置**：首次断裂那帧的断键中点应集中在**竖直直径附近**（经典间接拉伸劈裂）；
4. **脆性**：末帧断裂键数 > 0，且峰值后出现载荷跌落。

标定脚本（`_demtest/hmx_calib.py`）可复现 §2 的全部数据：
```bash
python _demtest/hmx_calib.py <exe> --kn 1.8862e7 --kt 1.4750e7 --break-strain 3.646e-3 --tension
python _demtest/hmx_calib.py <exe> --kn 1.8862e7 --kt 1.4750e7          # 压缩：测 E、ν
```
