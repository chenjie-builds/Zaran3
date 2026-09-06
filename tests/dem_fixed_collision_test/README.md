# DEM 测试算例：活动粒子碰撞固定粒子

粒子 0 以 1 m/s 沿 x 方向撞向粒子 1。粒子 1 的 CSV `motion_type` 字段为 0，
因此它参与接触计算但不进行运动积分。验证脚本检查粒子 0 反弹，且粒子 1
的位置和速度保持不变。

```bash
Zaran3.10.2.exe path/to/dem_fixed_collision_test
python verify.py --sim ./backup
```
