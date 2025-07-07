% setup_python.m - 配置MATLAB与Python环境的连接
% 用于调用scipy.special.wofz函数

clear; clc;


% 设置Python路径（根据您的Anaconda安装位置调整）
python_path = 'D:\ProgramData\anaconda3\envs\plasma\python.exe';

try
    % 检查Python路径是否存在
    if ~exist(python_path, 'file')
        fprintf('错误: Python路径不存在: %s\n', python_path);
        fprintf('请检查Anaconda安装路径和conda环境名称\n');
        return;
    end
    
    fprintf('Python路径: %s\n', python_path);
    
    % 配置Python环境
    pe = pyenv('Version', python_path);
    
    if isempty(pe.Executable)
        fprintf('错误: 无法配置Python环境\n');
        return;
    end
    
    % 导入scipy
    py.importlib.import_module('scipy');
    
    % 导入scipy.special
    py.importlib.import_module('scipy.special');
    
    % d导入wofz函数
    scipy_special = py.importlib.import_module('scipy.special');

    test_z = 1e03 + 1e03*i;
    result = 1i* sqrt(pi) * scipy_special.wofz(test_z);
    fprintf('result = %.6e + %.6ej\n', real(result), imag(result));
  
catch ME
    fprintf('错误: %s\n', ME.message);
end