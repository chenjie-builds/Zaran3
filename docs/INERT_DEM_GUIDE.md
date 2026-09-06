# Zaran3 惰性颗粒 DEM 入门

## 功能范围

当前模块求解三维球形惰性颗粒的平动、转动、颗粒接触、无限平面墙接触和永久
弹性键合。它没有热传导、化学反应或 CFD 双向耦合。永久键合移植了 LSM 中
“节点状态 + 连接拓扑 + 法切向弹性作用”的通用力学骨架，但不是全部 LSM 模型。

## 程序结构

一次 DEM 计算由以下对象组成：

| 层次 | 类型 | 职责 |
| --- | --- | --- |
| 应用入口 | `Application` | 读取 `zaran.toml`，按 `task.solver = "DEM"` 进入 DEM 分支 |
| 时间循环 | `DEMFieldSimulation` | 初始化、推进、记录时间和迭代数、定期输出 |
| 场对象 | `DEMField` | 按 Zaran 的 `Field` 接口组装数据、参数和求解器 |
| 状态容器 | `DEMFieldData` | 保存粒子、墙面和当前颗粒接触列表 |
| 求解器 | `DEMSolver` | 接触搜索、接触力、重力和运动积分 |
| 接触模型 | `LinearSpringDashpot`、`HertzMindlin` | 根据接触几何和相对速度计算法向、切向力 |
| 永久连接 | `DEMBond` | 保存端点拓扑、刚度、阻尼和当前连接状态 |
| 输入 | `ReadDEMParticle` | 从 CSV 建立初始粒子数组 |
| 连接输入 | `ReadDEMBond` | 用粒子 ID 读取永久连接并映射到数组索引 |

每个时间步的顺序为：清零力和力矩，计算永久键合力，搜索颗粒接触，计算颗粒接触力，计算墙面
接触力，清除已经结束的接触历史，叠加重力，最后用半隐式 Euler 法更新速度和
位置。

## 粒子数据

`DEMParticle` 直接保存一个球的状态：

| 字段 | 含义 | 单位 |
| --- | --- | --- |
| `id`, `group` | 用户编号与分组 | 无量纲 |
| `radius`, `mass` | 半径与质量 | m、kg |
| `inertia` | 球体转动惯量因子，实际转动惯量为 `inertia*m*r^2` | 无量纲 |
| `pos`, `vel`, `rotation`, `omega` | 位置、速度、累计转角、角速度 | m、m/s、rad、rad/s |
| `force`, `torque` | 当前步累积的合力和合力矩 | N、N·m |
| 四个材料参数 | 杨氏模量、泊松比、摩擦系数、恢复系数 | Pa、无量纲 |
| `active`, `kinematic` | 区分固定、动力学和运动学粒子 | 布尔值 |

输入 CSV 的必需列顺序为：

```text
id,group,radius,mass,px,py,pz,vx,vy,vz,ox,oy,oz
```

其后可以只增加 `motion_type`，也可以增加完整的五列：

```text
motion_type,young_modulus,poisson_ratio,friction_coeff,restitution_coeff
```

`motion_type=0` 为固定粒子，`1` 为动力学粒子，`2` 为按输入速度运动且轨迹不受
力改变的运动学粒子。因此每行只允许 13、14 或 18 个数值。`mass=0` 时，初始化阶段按

\[
m=\frac{4}{3}\pi r^3\rho
\]

计算质量，其中密度来自 `dem.density`。粒子行没有给出四个材料参数时，使用
控制文件中的全局默认值。

## 接触数据与历史

`DEMContact` 保存参与对象的数组索引、法向、重叠量、接触点、当前接触力和累计
切向位移。法向由粒子 A 指向粒子 B；墙面接触时由粒子指向墙。接触模型返回的
力是施加于粒子 A 的力。

切向弹簧需要跨时间步累计位移。邻域搜索会在每步重建接触列表，因此
`DEMSolver` 以 `(接触类型, A 索引, B/墙索引)` 为键保存切向历史，新接触从零
开始，已分离接触的历史在本步末删除。当前实现假定计算过程中不会改变粒子数组
的排列；如果以后支持运行中删除或重排粒子，接触键应改用稳定的粒子 ID。

## LSM 通用力学数据的 Zaran 映射

本次只映射与惰性力学直接相关的公共部分：

| LSM 数据 | Zaran3 数据 | 说明 |
| --- | --- | --- |
| `node_state_type%x,y` | `DEMParticle::pos.x/y` | 节点位置 |
| `node_state_type%theta` | `DEMParticle::rotation.z` | 二维累计转角；Zaran 内部使用三分量容器 |
| `node_state_type%vx,vy` | `DEMParticle::vel.x/y` | 平动速度 |
| `node_state_type%w` | `DEMParticle::omega.z` | 二维角速度 |
| `node_type%region,material` | `DEMParticle::group` | 当前只有一个分组字段 |
| `spring_item_type%ni,nj` | `DEMBond::idx_a,idx_b` | 运行时端点索引；输入使用粒子 ID |
| LEP `r0` | `DEMBond::rest_length` | 无应力长度，填 0 时取初始端点距离 |
| LEP `Kn*`, `Kt1` 的线弹性部分 | `normal_stiffness`, `tangential_stiffness` | 当前法向拉压共用一个刚度 |
| 节点相互作用力与力偶 | `DEMParticle::force`, `torque` | 每步由连接、接触和外力累加 |

当前连接模型按端点相对速度积分切向位移，并把切向力产生的力矩施加到两端。它
没有移植 LSM 中的非线性分支、拉压不同刚度、断裂能、材料切换、连接重建或热学
项，因此当前结果不能称为与 LSM 全模型数值等价。

## 永久弹性键合输入

在 `[dem]` 中设置 `bond_file = "bonds.csv"` 后，连接文件每行包含：

```text
id,particle_a_id,particle_b_id,rest_length,kn,kt,cn,ct,active
```

端点 ID 对应粒子 CSV 的 `id`。刚度单位为 N/m，阻尼单位为 N·s/m。程序会拒绝
不存在的端点、自连接、重复 ID、重复端点对和负参数。每次备份还会生成
`bonds.dat`，记录伸长量和当前连接力；VTP 中用线单元显示连接。

## 接触模型

线性模型的弹性法向力大小与重叠量成正比，阻尼系数由恢复系数、等效质量和法向
刚度得到。Hertz–Mindlin 模型的 Hertz 弹性法向力与重叠量的 3/2 次方成正比。
两种模型的切向力均受 Coulomb 上限约束：

\[
\lVert F_t\rVert \leq \mu\lVert F_n\rVert.
\]

固定粒子参与接触时，它被视为无限质量约束，等效质量取活动粒子的质量；两个
活动粒子接触时使用约化质量。墙面使用等效的大半径、大质量虚拟粒子向接触模型
提供材料参数。

## 控制文件

常用配置如下：

```toml
[task]
simulation = "SOLVE_FIELD"
solver = "DEM"

[dem]
time_step = 1e-5
end_time = 0.4
max_iter = 40000
output_iter = 1000
contact_model = "LinearSpringDashpot" # 或 HertzMindlin
particle_file = "particles.csv"
# 可选：bond_file = "bonds.csv"
density = 2500.0
young_modulus = 1e7
poisson_ratio = 0.3
friction_coeff = 0.4
restitution_coeff = 0.9
gravity_x = 0.0
gravity_y = -9.81
gravity_z = 0.0
```

求解器会拒绝非正时间步、非正密度或杨氏模量、非法泊松比、摩擦系数和恢复系数，
也会拒绝未知接触模型名称。

## 验证算例

`tests/dem_bounce_test` 验证单球自由下落和墙面反弹。当前 Debug 构建的首次反弹
峰值为 0.164044 m，解析参考为 0.163900 m。

`tests/dem_fixed_collision_test` 验证活动粒子可以与 `motion_type=0` 的固定粒子接触，
并确认固定粒子的状态不被积分修改。

`tests/dem_bond_test` 验证永久弹性键合、粒子 ID 到数组索引的映射、固定端约束、
连接状态备份及 VTP 连接线输出。

`tests/dem_brazilian_disk_test` 使用 85 个圆盘节点、222 条连接和上下运动学压头
验证弹性巴西圆盘响应。当前结果的竖向直径变化为 −2.458×10⁻⁴ m，圆心横向
平均应变为 1.526×10⁻³，上下压头反力分别为 ±5.20414 N，反力不平衡为 0。

## 当前限制

当前邻域搜索每步重建 KD 树；墙面接触仍对全部粒子和全部墙逐一检查。运动积分
为半隐式 Euler。累计转角尚未用于三维有限转动姿态更新。随机装填采用拒绝采样，适合小型
低体积分数初始场；高密度装填需要专门的生成或沉降过程。
