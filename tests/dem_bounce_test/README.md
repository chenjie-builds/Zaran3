# DEM 测试算例：弹跳球 (dem_bounce_test)

## 问题描述

| 参数 | 值 |
|---|---|
| 粒子半径 r | 0.01 m |
| 粒子密度 ρ | 2500 kg/m³（`mass=0` 时自动计算质量 m ≈ 0.01047 kg）|
| 初始位置 | (0, 0.20, 0)  m |
| 初始速度 | 静止 |
| 重力 | (0, −9.81, 0)  m/s² |
| 地板位置 | y = 0.0 |
| 接触模型 | LinearSpringDashpot |
| 杨氏模量 | 1×10⁷ Pa |
| 泊松比 | 0.30 |
| 恢复系数 | 0.90 |
| 时间步长 dt | 1×10⁻⁵ s |
| 仿真终止时间 | 0.40 s |

## 目录结构

```
dem_bounce_test/
├── zaran.toml          # Zaran 输入配置
├── particles.csv       # 初始粒子数据
├── verify.py           # Python 验证脚本（镜像 C++ 算法）
└── README.md           # 本文件
```

## 如何运行

1. **跑 Zaran 仿真**（将程序工作目录指向本文件夹，输入文件为 `zaran.toml`）
2. **跑验证脚本**：

```bash
# 仅打印预期结果
python verify.py

# 同时与 Zaran 输出的 backup 目录比对
python verify.py --sim ./backup
```

## 预期物理结果（解析估算）

| 量 | 值 |
|---|---|
| 落地时刻 | ≈ 0.1968 s |
| 落地速度 | ≈ 1.93 m/s（向下）|
| 弹跳后速度 | ≈ 1.74 m/s（向上）|
| 弹跳后最高质心 y | ≈ **0.1639 m** |

如果代码正确，首次反弹最高质心位置应接近 **0.164 m**。运行程序后执行
`python verify.py --sim ./backup`；脚本读取每个输出时刻的 `particles.dat`，
确认发生反弹，并检查反弹高度与解析参考值的误差。
