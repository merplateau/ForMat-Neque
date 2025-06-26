% 测试k_solver_function
fprintf('=== 测试k_solver_function ===\n');

% 测试默认参数
fprintf('\n1. 使用默认参数:\n');
[result_real, result_imag] = k_solver_function();
fprintf('结果: Re = %.6f, Im = %.6f\n', result_real, result_imag);

% 测试自定义参数
fprintf('\n2. 使用自定义参数:\n');
[result_real, result_imag] = k_solver_function('ne', 2e18, 'B0', 0.3, 'Te_eV', 5);
fprintf('结果: Re = %.6f, Im = %.6f\n', result_real, result_imag);

% 测试参数扫描
fprintf('\n3. 参数扫描:\n');
ne_values = [1e18, 2e18, 5e18];
for i = 1:length(ne_values)
    [result_real, result_imag] = k_solver_function('ne', ne_values(i));
    fprintf('ne = %.1e: Re = %.6f, Im = %.6f\n', ne_values(i), result_real, result_imag);
end

fprintf('\n测试完成!\n'); 