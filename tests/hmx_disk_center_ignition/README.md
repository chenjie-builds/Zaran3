# HMX 圆盘中心点火回归算例

该算例生成二维六角格点圆盘，在圆心半径 2.1 mm 内施加 1800 K 热点，检查
Arrhenius 热解、JWL 压力、气固共存固相核退缩燃烧、相间传热、沿连接表面燃烧、
反应波向外传播及反应弱化断键。

HMX 参数为：`rho=1890 kg/m3`、`E=28.764 GPa`、`nu=0.2`、
`c=1004.62 J/(kg K)`、`lambda=0.5358 W/(m K)`、`Q=5.53e6 J/kg`、
`Z=4.78e12 1/s`、`Ea=143.9 kJ/mol`（输入换算为 `Ea/R=17307.19 K`）。JWL 的
`A/B/R1/R2/omega` 取自 LSM 源码注释。用户没有提供压力燃速系数，`burn_a` 使用了
加速回归值，因此本算例验证算法链路，不用于预测真实点火延迟或爆速。

当前基准结果：91 个格点全部进入反应，圆心平均反应度为 1.0，外圈 60 个格点发生
反应，最高气压为 1.0e5 Pa（测试限幅），240 条连接中 234 条失效；输出中
`max_core_dalpha=6.002e-4`、`max_neighbor_dalpha=1.118e-3`，证明两种燃烧通道均被
执行；0–8000 步输出无 NaN/Inf。

```powershell
python generate_case.py
E:\ZARAN_DEM\ZARAN\Zaran3\bin\Release\Zaran3.10.2.exe
python verify.py
```
