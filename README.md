# ForMat-Neque

这是一个包含MATLAB和Fortran代码的科学计算项目。

## 项目结构

- `k_solver.m` - 主要的求解器
- `k_solver_scan_B.m` - B参数扫描
- `k_solver_scan_omega.m` - omega参数扫描
- `solve_k_parallel.m` - 并行求解器
- `test_omega_ci_singularity.m` - 奇点测试
- `main.f90` - Fortran主程序
- `equations.tex` - 数学公式文档
- `faddeeva/` - Faddeeva函数库
- `Faddeeva_MATLAB/` - MATLAB版本的Faddeeva函数
- `zetaf/` - Zeta函数库

## 依赖项

- MATLAB
- Fortran编译器（如gfortran）

## 使用方法

1. 确保MATLAB已安装
2. 运行主求解器：`k_solver.m`
3. 根据需要运行其他脚本

## 许可证

请查看各个子目录中的license.txt文件了解具体许可信息。 