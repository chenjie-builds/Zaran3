# 10000 格点 HMX 圆盘中心点火算例

圆盘半径 20 mm、厚度 100 mm，使用严格中心反演对称的三角格点。材料参数、JWL、
气固共存核心燃烧和输出字段与小型 HMX 回归一致。圆心附近 1 mm 半径区域初始化为
1000 K，使给定 Arrhenius 参数在当前时间步下渐进点火，而不是首步直接燃尽。

```powershell
python generate_case.py
python verify.py
E:\ZARAN_DEM\ZARAN\Zaran3\bin\Release\Zaran3.10.2.exe
python verify.py
```

修正后输出写入 `result_lsm_fixed/backup_lsm_fixed`。原先 `result/backup` 中使用加速
燃速和反应弱化参数生成的结果予以保留，便于对比，但不应再作为当前算例结果读取。

`burn_a=0.639, burn_b=1.19` 来自 LSM `initializers.f90` 的 PBX-1 参数。算例启用
动态 Delaunay/Voronoi 气相控制体；每 100 步重构一次拓扑。连接断裂使用
`surface_energy=1680 J/m2` 与 `grain_boundary_energy=4.2 J/m2` 的 LSM LEP 储能判据，
`bond_peak_strain/bond_break_strain` 仅作为没有表面能时的后备准则。高压压缩采用
`rnn=0.8, knn=10, exp=1`，法向和切向黏度均为 `100 Pa s`。
