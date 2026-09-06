# 10000 格点巴西圆盘算例

半径 50 mm 的惰性键合圆盘包含严格中心反演对称的 10000 个材料格点。上下运动学
粒子压头以 0.05 m/s 对称加载。压头粒子不计入“10000 格点”。

```powershell
python generate_case.py
python verify.py
E:\ZARAN_DEM\ZARAN\Zaran3\bin\Release\Zaran3.10.2.exe
python verify.py
```

当前参数用于验证大规模接触、键合、对称加载和 Release 计算能力；正式巴西劈裂
强度计算还需用试验标定键刚度、峰值应变、断裂应变和阻尼。
