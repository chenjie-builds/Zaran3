# 键合的压缩失效 + 逐键强度异质性（让巴西盘真正"被压碎"）

对应上一轮诊断（`docs/FRAGMENTATION_SPLASH.md`）的两条根因修复：

* **① 压缩/屈曲失效上限** —— 键原本只有拉伸断裂判据，压缩侧完全无上限，
  实测压缩应变累积到 **8.2%**，键储能的 **99.8% 是压缩能**（37.67 J vs 拉伸 0.06 J）。
  断裂瞬间这个虚高的水库一次性释放，碎屑以 112 m/s 飞溅。
* **② 逐键强度异质性（Weibull）** —— 强度全局均匀时损伤只能弥散、载荷不跌；
  引入离散度后破坏才能局部化。

## 1. 参数

| 参数 | 默认 | 含义 |
|---|---|---|
| `dem.bond_compression_strength_ratio` | **8.0** | 抗压/抗拉强度比 σ_c/σ_t；压缩断裂应变 = ratio × `bond_break_strain`。`0` = 关闭压缩判据 |
| `dem.bond_break_strain_compression` | 0 | 直接指定压缩断裂应变；>0 时优先于 ratio |
| `dem.bond_weibull_modulus` | 0（关闭） | Weibull 形状参数 m；`0` = 均质（所有键同一阈值） |
| `dem.bond_weibull_seed` | 20260917 | 抽样种子（保证可复现） |

**为什么默认 ratio = 8**：脆性固体的 σ_c/σ_t 典型为 8–15。因为阈值由抗拉阈值**派生**，
所以只有显式设了 `bond_break_strain` 的算例会受影响——本仓库里只有
`hmx_brazilian_disk_10k`（实测 8.2% 压缩）会改变行为，其余算例（`bond_break_strain`
为默认 1e30，或压缩仅 3e-4）逐位不变。

**为什么 Weibull 默认关闭**：它会把宏观强度按最弱链效应降下来（需要重新标定），
且会改动已按精确值验证过的破坏算例，所以做成可选材料属性。

## 2. 实现要点

### 2.1 压缩判据的位置与形式

放在 `CalcBondForce` 中 `bond.elastic_energy` 算出之后、能量准则之前：

```cpp
const double compression_limit = GetBondBreakStrainCompression() > 0.0
    ? GetBondBreakStrainCompression() * strength_scale
    : GetBondCompressionStrengthRatio() * limit;          // limit = bond_break_strain * s
if (compression_limit > 0.0 && compressive_strain > compression_limit)
{
    bond.dissipated_fracture_energy += bond.elastic_energy;
    bond.elastic_energy = 0.0;
    bond.damage = 1.0;  bond.active = false;  bond.force_a.setZero();
    continue;
}
```

三条设计决定：

1. **必须独立于 `fracture_energy` 分支。** 能量准则只在 `extension >= 0` 时生效
   （LSM LEP），若把压缩判据也塞进那个分支，凡是用 `surface_energy` 的能量口径算例
   就完全没有压缩强度。现在即使 `fracture_energy > 0`，压缩上限仍然生效。
2. **取脆性判据（达阈值即断），不给软化段。** 压缩侧不用双线性软化：粉碎区的
   渐进失效来自应变场的非均匀性，而不是材料软化。这样 `damage` 的语义仍专属于拉伸。
3. **断裂时把 `elastic_energy` 清零。** 这是顺带修掉的一个输出缺陷：原先断开的键
   会保留最后一步的 `elastic_energy`，导致对全部键求和时与

   `dissipated_fracture_energy` **重复计数**（实测两者会收敛到同一個数字：
   `Ebond = 0.8433` vs `Ediss = 0.8430`）。清零后 `Σelastic_energy` 只统计完好键，
   能量账才有意义。

### 2.2 逐键阈值缩放口径

`strength_scale = s` 时，应力与应变同乘 s，因此：

| 阈值 | 缩放 | 理由 |
|---|---|---|
| 拉伸断裂应变 | × s | σ = E·ε，线性 |
| 双线性峰值应变 | × s | 保持软化段在阈值中的相对位置 |
| 压缩断裂应变 | × s | 同上 |
| 断裂能（能量准则） | × s² | 能量密度 ∝ σ·ε ∝ s² |

### 2.3 确定性抽样

```cpp
key = SplitMix64(seed ^ (id_a * 0x9E3779B97F4A7C15) ^ ((id_b + 0x165667B19E3779F9) * 0xC2B2AE3D27D4EB4F));
u   = ((key >> 11) + 0.5) / 2^53            // (0,1)
s   = λ · (−ln(1−u))^{1/m},   λ = 1/Γ(1+1/m)
```

* 抽样键是**粒子 id 对**而不是数组下标 ⇒ 同一物理连接无论输入顺序如何都拿到同一强度。
* `λ = 1/Γ(1+1/m)` 把均值归一到 1 ⇒ **不改变**先前按宏观量标定出的整体强度，只加离散度。
* `s` 截断在 [0.01, 100] 作保护。

### 2.4 输出

* `bonds.dat` 末尾追加两列：`maximum_compressive_strain`、`strength_scale`
  （末尾追加，按列名解析的脚本不受影响）。
* VTP CellData 增加 `bond_max_compression`、`bond_strength_scale`（受 `output.bond_details` 控制）。

## 3. 单元级验证（3 个新算例，均纳入回归）

### 3.1 `dem_compression_fracture_test` —— 压缩断裂应变的解析校验

两个粒子**都是规定运动**、以 ±0.05 m/s 相向而行，所以压缩应变解析已知：

```
strain(t) = 2·v·t/L0 = 2×0.05×t/0.1 = t          [1/s · s]
```

两球半径 0.02 < 间距 0.1，全程不可能接触 ⇒ 只有键合判据可能触发。
`bond_break_strain = 1e-3`、ratio = 8 ⇒ 阈值 8e-3 ⇒ 预测在 t = 8 ms（iter 8000）断裂。

```
comp threshold   : ratio 8 x tensile 1.0e-03 = 8.0000e-03
strain rate      : 1.0000 1/s  -> predicted break at t = 8.0000 ms (iter 8000)
measured failure : maximum_compressive_strain = 8.000000000e-03
fracture energy  : dissipated 3.200000e-02 J  vs  0.5*kn*(strain*L0)^2 = 3.200000e-02 J
last intact frame: iter=8000 at strain 7.999000000e-03; first broken frame: iter=8500
```

**断裂应变与派生阈值精确吻合，耗散能等于解析储能。**

### 3.2 `dem_compression_fracture_off_test` —— 控件

同几何同载荷、`bond_compression_strength_ratio = 0`：键压缩到 **8.499e-03 = 8.5 倍抗拉阈值**
仍然存活、零耗散、零损伤 ⇒ 证明 3.1 的断裂确实由新判据触发，且关闭后旧行为逐位保留。

### 3.3 `dem_weibull_strength_test` —— 分布与可复现性

400 粒子、760 键、无载荷（纯材料属性）：`verify.py` 用 Python 重算
splitmix64 + Weibull 公式逐条比对。

```
bonds                : 760
mean strength_scale  : 1.00578   (target 1.0)
min / max            : 0.40047 / 1.53767
KS distance          : 0.0365   (Weibull m=6, 1% 临界值 ~0.059)
recipe agreement     : worst relative deviation 4.89e-06  (= 6 位有效数字精度)
```

⇒ 分布正确、均值归一、逐条可复现（种子确实进入哈希）。

## 4. 在 HMX 巴西盘上的效果

`_demtest/crush_compare.py` 的统一口径（`Ebond` 只统计完好键）：

| 配置 | 峰值载荷 | 峰后 min/peak | max 压缩应变 | KE 峰值 | Vmax 峰值 | 末 `rmax/R` | 压板功 | 断键 | 撞墙粒子 |
|---|---|---|---|---|---|---|---|---|---|
| 旧（无压缩判据） | 8.47→**27.0 kN 单调升** | **1.00（无跌落）** | **8.16e-02** | 8.570 J | **122.2 m/s** | 2.90 | 83.29 J | 21547 | **537** |
| A: ratio 8 | 4.337 kN | **0.000** | 1.073e-02 | 1.435 J | 65.0 m/s | 1.79 | 2.31 J | 5980 | **0** |
| B: ratio 8 + Weibull 6 | 3.828 kN | **0.000** | 1.410e-02 | 1.052 J | 49.4 m/s | 1.81 | 2.04 J | 5148 | **0** |
| C: ratio 40 | 19.38 kN（仍在升） | 1.00 | **4.760e-02** | 0.004 J | 12.1 m/s | 1.04 | 21.15 J | 3270 | **0** |

### 4.1 关键改善（A 配置）

* **压缩应变 8.16e-2 → 1.07e-2**（7.6 倍），不再有虚高的压缩水库。
* **载荷有了真正的峰与崩塌**：4.337 kN（iter 22500）→ 约 0，`min/peak = 0.000`。
  旧行为是单调升到 27 kN（3.2 倍）后爆炸、再反涨到 34.5 kN。
* **Vmax 112 → 25 m/s**，KE 8.57 → 1.27 J，压板做功 40.4 → 2.31 J（同 6 万步口径 17 倍）。
* **碎片不再撞墙**（537 → 0），`rmax/R` 2.90 → 1.79。
* 瞬时峰值 65 m/s 出现在崩塌那一帧，随后稳定在 25 m/s。

**能量账闭合**（A，iter=60000）：

```
W 2.3089 J = KE 1.2720 (55%) + Ebond(完好键) 0.00025 + 碎裂耗散 0.8430 (37%)
             + 阻尼/热 0.1936 (8.4%)
```

### 4.2 ② 的效果（B vs A）

Weibull 让破坏**更早发生、峰值更低**：崩塌从 iter 25000 提前到 22500，
峰值 4.337 → 3.828 kN，KE 1.435 → 1.052 J，起裂帧的断键 4 → 8 条。
定量上由 3.3 的分布校验保证。

### 4.3 ⚠️ 必须知道的一件副作用：failure mode 换了

| 配置 | 起裂位置 | 判定 |
|---|---|---|
| 旧 / C(ratio 40) | iter=30000，1250 条断键 **100% 落在 \|x\|<0.25R**、纵向跨 ±0.62R | ✅ 经典**中心竖直劈裂** |
| A / B（ratio 8） | iter=7500，4/8 条断键位于 **x≈0 且 \|y\|>0.98R** | ⚠️ 起裂改在**压板接触点**，随后损伤铺满全盘（弥散压碎） |

也就是说：**ratio 8 换来"压碎"和干净的载荷崩塌，但抢在经典巴西劈裂之前**。
（C 的 ratio 40 保留了劈裂：压缩被精确截在 4.760e-02 = 40 × 1.19e-3，
载荷、KE、`rmax/R` 都与旧行为几乎一致，只是水库上限降了 3 倍。）

### 4.4 两个选择的取舍，以及真正的根因

要让**劈裂先发生**，需要压缩阈值 > 起裂时的压缩应变 4.20e-2 ⇒
`ratio > 4.2e-2 / 1.19e-3 ≈ 35`。但脆性固体的 σ_c/σ_t 只有 8–15，
**35 已经超出物理范围**。

之所以需要这么高的比值，根因是**几何**：

* 试件厚度 `t = 2r = 0.945 mm`，直径 `D = 99 mm` ⇒ **t/D = 0.0095**（约一个粒子直径）。
  薄片的压板接触应力被极度放大（名义 `2P/πDt` 才 29 MPa 时，压板下方已压到 σ_c）。
  真实巴西盘标准试样是 **t/D ≈ 0.5**。
* 压板是 **133 个相切的离散球**（半径 3.78e-4），球下的局部应力集中进一步放大了压碎趋势。

所以正确的做法是：**用真实厚度（t/D ≈ 0.2–0.5）**，或把加载面做成连续面
（压板球相互重叠 / 加厚多排 / 平面墙）。届时 σ_c/σ_t = 8 就能让
"先中心劈裂、随后接触区压碎"按正确顺序出现。

## 5. 结论

* ① + ② 已实现、编入、并通过 **19/19 全量回归**（含 3 个新算例；`alpha=0.0956179`、
  `broken_bonds=240`、`volume=8.999992000`、`fracture_energy=0.577350` 等精确值逐位未变）。
* 巴西盘的**飞溅问题被消除**：压缩水库被截断，载荷出现峰值与崩塌，
  碎片速度进入物理范围，能量账闭合到 8%。
* 副作用是 failure mode 由"劈裂主导"变为"压碎主导"（ratio 8）；想保留劈裂需
  ratio ≥ 35，而这是几何过薄造成的，建议改试样厚度而不是调参数。

## 6. 复现命令

```bash
# 单元级
python _demtest/run_committed_tests.py bin/Release/Zaran3.10.2.exe \
    dem_compression_fracture_test dem_compression_fracture_off_test dem_weibull_strength_test

# 巴西盘对照（A/B/C）
python _demtest/crush_effect.py bin/Release/Zaran3.10.2.exe A_ratio8 --ratio 8 --weibull 0
python _demtest/crush_effect.py bin/Release/Zaran3.10.2.exe B_ratio8_weib6 --ratio 8 --weibull 6
python _demtest/crush_effect.py bin/Release/Zaran3.10.2.exe C_ratio40 --ratio 40 --weibull 0

# 汇总表
python _demtest/crush_compare.py "旧=tests/hmx_brazilian_disk_10k" "A=_demtest/crush_A_ratio8" ...
```
