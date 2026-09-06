# DEM 测试算例：永久弹性键合

两个粒子由永久线性弹性键连接。ID 10 的粒子固定，ID 20 的粒子初始距离
大于键的无应力长度。算例验证弹簧把活动端拉回，同时固定端保持不动。

`bonds.csv` 的列为：

```text
id,particle_a_id,particle_b_id,rest_length,kn,kt,cn,ct,active
```

`rest_length=0` 表示使用读取时两个粒子之间的距离。
