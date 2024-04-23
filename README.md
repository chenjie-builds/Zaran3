# Zaran -- A Total Automated CFD Solver
### Author： Chen Jie
### Email: chenjiework@live.com

## 第三方库
- Eigen3
- spdlog
- VTK


## Windows 下编译
- 编译环境： Visual Studio 2019/2022
- [1]安装[vcpkg](https://github.com/microsoft/vcpkg/blob/master/README_zh_CN.md)
- [2]利用vcpkg安装Eigen3，spdlog，VTK库
- [3]修改项目的[CMakeLists.txt](CMakeLists.txt)，将 ***CMAKE_TOOLCHAIN_FILE*** 修改为: ***[vcpkg 安装位置]]/scripts/buildsystems/vcpkg.cmake***
- [4]将[vcpkg安装位置]/installed/\x64-windows\bin 添加到环境变量，即可不用拷贝dll文件

## Reference
- [1]刘君, 韩芳, 夏冰. 有限差分法中几何守恒律的机理及算法[J]. 空气动力学学报, 2018, 036(006): 917-926.
- [2]刘君, 韩芳. 有限差分法中的贴体坐标变换[J]. 气体物理, 2018, 003(005): 18-29.
- [3]陈洁. 基于三相邻节点的非结构网格有限差分方法[D]. 大连理工大学, 2019.
- [4]刘君, 陈洁, 韩芳. 基于离散等价方程的非结构网格有限差分法[J]. 航空学报, 2020(1).
- [5]刘君, 魏雁昕, 陈洁. 基于非结构网格有限差分法的扎染算法[J]. 航空学报, 2021, 42(07): 258-270.
## [Change Log](change.md)