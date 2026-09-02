#!/usr/bin/env python3
"""
DEM 弹跳球验证脚本 – dem_bounce_test
======================================
算例描述
---------
一个球形粒子（半径 r=0.01 m，密度 ρ=2500 kg/m³）从高度 y=0.20 m 处
由静止开始自由下落，碰到地板（y=0.0）后弹起。

本脚本做两件事：
  1. 解析估算（不涉及接触刚度，仅用于量级参考）
  2. 逐步复现 C++ 代码中的 **完全相同算法**
     （LinearSpringDashpot + 前向欧拉 + 墙面法向惯例）
     并将结果写入 verify_result.csv，可与 Zaran 输出的
     backup/iter=xxx/particles.csv 逐行比对。

⚠️  重要：本脚本中的 `wall_normal` 约定与 C++ `CalcWallForce` 保持
    一致（法向 = wall.normal，指向计算域内侧，即朝上）。若物理上
    正确的法向应使粒子被推离壁面（朝上），但接触力是 -Fn*normal
    （向下），则说明代码存在符号错误，请对照输出结果判断。

用法
----
  python verify.py [--sim <backup_dir>]

  --sim <path>  可选：Zaran 仿真输出的 backup 根目录，
                脚本会读取其中最后一个 iter 文件夹的
                particles.csv 并与本地模拟结果比对。
示例
----
  python verify.py
  python verify.py --sim ./backup
"""

import math
import csv
import os
import sys
import argparse

# ──────────────────────────────────────────────
# 参数（与 zaran.toml / DEMWall 默认值保持一致）
# ──────────────────────────────────────────────
DT       = 1e-5          # 时间步长 (s)
T_END    = 0.30          # 终止时间 (s)
N_STEPS  = int(T_END / DT)

G        = -9.81         # 重力加速度 y 分量 (m/s²)
RHO      = 2500.0        # 粒子密度 (kg/m³)
R        = 0.01          # 粒子半径 (m)
Y0       = 0.20          # 粒子初始质心 y (m)

E_BALL   = 1e7           # 球杨氏模量 (Pa)
NU_BALL  = 0.30          # 球泊松比
E_WALL   = 1e9           # 墙杨氏模量 (Pa，DEMWall 默认)
NU_WALL  = 0.30          # 墙泊松比
R_WALL   = 1e10          # 虚拟墙半径 (m)
M_WALL   = 1e30          # 虚拟墙质量 (kg，视为无穷大)
EREST    = 0.90          # 恢复系数

Y_FLOOR  = 0.0           # 地板 y 坐标

PI = math.pi

# ──────────────────────────────────────────────
# 辅助函数：LinearSpringDashpot 参数计算
# （与 LinearSpringDashpot.cpp 完全对应）
# ──────────────────────────────────────────────
def calc_kn_cn(E_a, nu_a, m_a, r_a, E_b, nu_b, m_b, r_b, e_rest):
    """计算法向刚度 kn 和阻尼系数 cn"""
    Ea_eff = E_a * (1.0 - nu_a**2)
    Eb_eff = E_b * (1.0 - nu_b**2)
    E_star = 1.0 / (1.0/Ea_eff + 1.0/Eb_eff)

    R_star = (r_a * r_b) / (r_a + r_b)
    k_n    = 2.0 * E_star * R_star

    m_eff  = (m_a * m_b) / (m_a + m_b)
    e = max(e_rest, 1e-3)
    ln_e   = math.log(e)
    beta   = -ln_e / math.sqrt(PI**2 + ln_e**2)
    c_n    = 2.0 * beta * math.sqrt(m_eff * k_n)

    return k_n, c_n

# ──────────────────────────────────────────────
# 解析参考值（无接触刚度影响，仅自由落体+理想弹跳）
# ──────────────────────────────────────────────
def analytical_reference(r, y0, g, e_rest):
    h_drop = y0 - r          # 质心从 y0 落到刚接触地板（y=r）的距离
    v_impact = math.sqrt(2.0 * abs(g) * h_drop)   # 接触瞬间速度 (m/s，向下)
    v_bounce = e_rest * v_impact                   # 弹跳速度 (m/s，向上)
    h_bounce = v_bounce**2 / (2.0 * abs(g))        # 弹跳后质心高度 = r + h_bounce
    t_fall   = math.sqrt(2.0 * h_drop / abs(g))

    print("=" * 60)
    print("【解析参考值（理想自由落体 + 完全弹性恢复系数）】")
    print(f"  下落距离 h_drop      = {h_drop:.6f} m")
    print(f"  接触前速度 v_impact  = {v_impact:.6f} m/s  (向下)")
    print(f"  弹跳后速度 v_bounce  = {v_bounce:.6f} m/s  (向上)")
    print(f"  弹跳后最高质心 y     = {r + h_bounce:.6f} m")
    print(f"  理论落地时刻 t_fall  ≈ {t_fall:.6f} s")
    print("  注：实际数值解因时间步长、接触力积分而略有偏差")
    print("=" * 60)
    return t_fall, v_impact, v_bounce

# ──────────────────────────────────────────────
# 主模拟循环（镜像 C++ DEMSolver 算法）
# ──────────────────────────────────────────────
def simulate():
    # 粒子初始状态
    m = (4.0/3.0) * PI * R**3 * RHO
    inertia = 0.4   # I = inertia * m * r²

    pos_y  = Y0
    vel_y  = 0.0

    # 接触持久状态（切向位移 delta_t，这里 1D 无切向，忽略）
    k_n, c_n = calc_kn_cn(E_BALL, NU_BALL, m,     R,
                           E_WALL, NU_WALL, M_WALL, R_WALL,
                           EREST)

    print(f"\n【接触参数】")
    print(f"  等效法向刚度 k_n = {k_n:.4e} N/m")
    print(f"  法向阻尼系数 c_n = {c_n:.4e} N·s/m")
    print(f"  粒子质量      m  = {m:.6e} kg")
    print(f"  固有频率   ω_n  = {math.sqrt(k_n/m):.4f} rad/s")
    print(f"  稳定时间步   Δt < {2*math.sqrt(m/k_n):.2e} s   (当前 DT={DT:.0e})")
    print()

    # ──────────────────────────────────────────
    # 记录轨迹
    # ──────────────────────────────────────────
    records = []         # (time, pos_y, vel_y, force_y)
    t = 0.0
    bounce_count = 0
    max_y_after_bounce = Y0
    t_first_contact    = None
    vel_at_contact     = None

    for step in range(N_STEPS):
        # ── 1. ZeroForce ──
        force_y = 0.0

        # ── 2. CalcWallForce（地板 y=0，法向 n=(0,1,0)，即 +y 方向）──
        #   C++ 中 wall.normal 指向计算域内侧（上方）
        #   signed_dist = n · (pos - point) = pos_y - 0 = pos_y
        d       = pos_y        # 质心到地板的有符号距离
        overlap = R - d        # 压入量
        if overlap > 0.0:
            # 记录首次接触
            if t_first_contact is None:
                t_first_contact = t
                vel_at_contact  = vel_y

            # 相对速度（法向分量）：v_n_rel = (vel_a - vel_b) · n
            #   wall vel=0，n=(0,1,0) → v_n_rel = vel_y
            v_n_rel = vel_y

            Fn_mag = k_n * overlap - c_n * v_n_rel
            if Fn_mag < 0.0:
                Fn_mag = 0.0

            # ⚠️  C++ CalcWallForce 中 c.normal = wall.normal = (0,1,0) 上方
            #   force_n = -Fn_mag * normal = -Fn_mag * (+y)  → 向下(负 y)
            #   与物理直觉（地板应向上推球）相反 → 这是待验证的符号问题
            wall_contact_force_y = -Fn_mag   # 与 C++ 一致

            force_y += wall_contact_force_y

        # ── 3. CalcGravity ──
        force_y += m * G

        # ── 4. Integrate（前向欧拉）──
        acc_y  = force_y / m
        vel_y += acc_y  * DT
        pos_y += vel_y  * DT

        t += DT

        # 记录每 1000 步一次
        if step % 1000 == 0:
            records.append((t, pos_y, vel_y, force_y))

        # 检测弹跳（速度由负转正，且曾在接触中）
        if overlap > 0 and vel_y > 0 and bounce_count == 0:
            bounce_count += 1

        # 跟踪第一次弹跳后最大高度
        if bounce_count >= 1 and pos_y > max_y_after_bounce:
            max_y_after_bounce = pos_y

    return records, t_first_contact, vel_at_contact, max_y_after_bounce, k_n, c_n

# ──────────────────────────────────────────────
# 比对仿真输出
# ──────────────────────────────────────────────
def compare_with_sim(backup_dir, records):
    # 找最后一个 iter 文件夹
    iter_dirs = [d for d in os.listdir(backup_dir) if d.startswith("iter=")]
    if not iter_dirs:
        print(f"[比对] backup 目录 {backup_dir} 中未找到 iter= 文件夹")
        return
    iter_dirs.sort(key=lambda x: int(x.split("=")[1]))
    last_dir = os.path.join(backup_dir, iter_dirs[-1], "particles.csv")

    if not os.path.exists(last_dir):
        print(f"[比对] 文件不存在: {last_dir}")
        return

    with open(last_dir) as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    if not rows:
        print("[比对] particles.csv 为空")
        return

    p = rows[0]
    sim_py  = float(p["py"])
    sim_vy  = float(p["vy"])
    t_final = float(iter_dirs[-1].split("=")[1]) * DT * 1000  # 每 1000 步输出一次

    # 找本地模拟的最后一条记录
    local_t, local_py, local_vy, _ = records[-1]

    print("\n" + "=" * 60)
    print("【与 Zaran 仿真输出比对（最终时刻，粒子 0）】")
    print(f"  Zaran 输出   py={sim_py:.8f}  vy={sim_vy:.8f}")
    print(f"  Python 复现  py={local_py:.8f}  vy={local_vy:.8f}")
    dy = abs(sim_py - local_py)
    dv = abs(sim_vy - local_vy)
    print(f"  差值 |Δpy|={dy:.2e}  |Δvy|={dv:.2e}")
    if dy < 1e-6 and dv < 1e-6:
        print("  ✅ 两者高度吻合（数值算法完全一致）")
    else:
        print("  ❌ 差值较大，请检查算法或参数是否一致")
    print("=" * 60)

# ──────────────────────────────────────────────
# 主入口
# ──────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sim", default=None,
                        help="Zaran 仿真 backup 目录（可选，用于与本地结果比对）")
    args = parser.parse_args()

    # 1. 打印解析参考值
    t_fall_ref, v_impact_ref, v_bounce_ref = analytical_reference(R, Y0, G, EREST)

    # 2. 运行 Python 模拟
    records, t_contact, v_contact, max_y, k_n, c_n = simulate()

    print("【Python 复现结果（镜像 C++ 算法）】")
    if t_contact is not None:
        print(f"  首次接触时刻       t ≈ {t_contact:.6f} s  （解析参考 {t_fall_ref:.6f} s）")
        print(f"  接触时粒子速度  vy ≈ {v_contact:.6f} m/s （解析参考 {-v_impact_ref:.6f} m/s）")
    else:
        print("  ⚠️  在模拟时间内粒子未接触地板！请检查参数。")
    print(f"  第一次弹跳后最高点 y_max = {max_y:.6f} m")
    print(f"  解析参考弹跳高度 (质心)  = {R + v_bounce_ref**2/(2*abs(G)):.6f} m")

    print()
    print("⚠️  注意：本算法中 wall_contact_force_y = -Fn_mag（向下）。")
    print("   物理上正确的地板接触力应向上，因此请观察粒子是否真的能弹起。")
    print("   若 max_y 始终 ≤ Y0（未弹起），说明 CalcWallForce 中")
    print("   c.normal 符号存在问题（应使用 -wall.normal 作为接触法向）。")

    # 3. 写出 Python 模拟轨迹到 CSV，方便与 Zaran 输出对比
    out_csv = os.path.join(os.path.dirname(__file__), "verify_result.csv")
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["time", "pos_y", "vel_y", "force_y"])
        for row in records:
            writer.writerow([f"{v:.8e}" for v in row])
    print(f"\n  轨迹已写入: {out_csv}")

    # 4. 可选：与 Zaran 输出比对
    if args.sim:
        compare_with_sim(args.sim, records)

if __name__ == "__main__":
    main()
