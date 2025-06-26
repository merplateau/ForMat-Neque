# Fortran调用MATLAB的方法

本项目提供了三种在Fortran中调用MATLAB `.m` 文件的方法。

## 方法1：MATLAB Engine API（推荐用于复杂交互）

### 优点：
- 直接内存交互
- 可以传递复杂数据结构
- 性能较好

### 缺点：
- 需要编译时链接MATLAB库
- 配置相对复杂

### 使用方法：

1. **编译程序**：
   ```bash
   make
   ```

2. **运行程序**：
   ```bash
   ./main_with_matlab
   ```

3. **如果需要修改MATLAB路径**，编辑 `Makefile` 中的 `MATLAB_ROOT` 变量。

## 方法2：系统调用（推荐用于简单调用）

### 优点：
- 简单易用
- 不需要特殊编译
- 适合一次性调用

### 缺点：
- 启动MATLAB较慢
- 无法直接传递复杂数据

### 使用方法：

1. **编译程序**：
   ```bash
   gfortran -o main_system_call main_system_call.f90
   ```

2. **运行程序**：
   ```bash
   ./main_system_call
   ```

## 方法3：函数化MATLAB代码

### 优点：
- 可以接受参数
- 返回计算结果
- 便于集成

### 使用方法：

1. **在MATLAB中测试函数**：
   ```matlab
   test_matlab_call
   ```

2. **在Fortran中调用**：
   ```fortran
   ! 创建MATLAB脚本
   open(unit=10, file='call_function.m', status='replace')
   write(10, *) "[result_real, result_imag] = k_solver_function('ne', 2e18);"
   write(10, *) "save('results.mat', 'result_real', 'result_imag');"
   write(10, *) "exit"
   close(10)
   
   ! 执行MATLAB
   call system("matlab -batch 'call_function'")
   ```

## 文件说明

- `main_with_matlab.f90` - 使用MATLAB Engine API的Fortran程序
- `main_system_call.f90` - 使用系统调用的Fortran程序
- `k_solver_function.m` - 函数化版本的k_solver
- `test_matlab_call.m` - 测试MATLAB函数的脚本
- `Makefile` - 编译配置文件

## 注意事项

1. **MATLAB路径**：确保MATLAB已正确安装并添加到PATH
2. **权限**：某些系统可能需要管理员权限来编译MATLAB Engine程序
3. **版本兼容性**：确保MATLAB版本与编译环境兼容

## 推荐使用场景

- **简单调用**：使用方法2（系统调用）
- **复杂交互**：使用方法1（Engine API）
- **参数化计算**：使用方法3（函数化）

## 故障排除

1. **编译错误**：检查MATLAB路径和库文件
2. **运行时错误**：确保MATLAB许可证有效
3. **路径问题**：确保所有文件在同一目录或正确路径下 