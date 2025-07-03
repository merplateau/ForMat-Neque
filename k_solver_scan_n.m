% 基于k_solver.m的参数扫描程序
% 扫描等离子体密度，磁场固定为0.5T
clear; clc;

% 基本常数
e = 1.602e-19;           % 元电荷 (C)
me = 9.11e-31;           % 电子质量 (kg)
mi = 40*1.67e-27;        % 质子质量 (kg)
eps0 = 8.854e-12;        % 真空介电常数 (F/m)
c = 3e8;                 % 光速 (m/s)
kB = 1.38e-23;           % 玻尔兹曼常数 (J/K)

% VASIMR ICRH 典型参数
f = 190e3;               % 驱动频率 (Hz) - 固定为190kHz
Te_eV = 3;               % 电子温度 (eV)
Ti_eV = 0.3;             % 离子温度 (eV)
Te = Te_eV * e / kB;     % 电子温度 (K)
Ti = Ti_eV * e / kB;     % 离子温度 (K)
Tperp_e = 10*Te;         % 电子垂直温度 (K)
Tpara_e = 1*Te;          % 电子平行温度 (K)
Tperp_i = 10*Ti;         % 离子垂直温度 (K)
Tpara_i = 1*Ti;          % 离子平行温度 (K)

% 固定参数
omega = 2*pi*f;          % 角频率 (rad/s)
B0 = 0.5;                % 磁场强度 (T)

% 扫描参数
density_values = linspace(1e17, 1e19, 100); % 密度范围 (m^-3)
num_n = length(density_values);

% 存储结果
k_real_results = zeros(1, num_n);
k_imag_results = zeros(1, num_n);
omega_ci_results = zeros(1, num_n);
omega_ce_results = zeros(1, num_n);
success_flag = zeros(1, num_n);

fprintf('开始等离子体密度扫描...\n');
fprintf('磁场: %.2f T\n', B0);
fprintf('密度范围: %.1e - %.1e m^-3\n', min(density_values), max(density_values));
fprintf('扫描点数: %d\n\n', num_n);

% 色散函数Z(ξ)的实现 - 使用faddeeva.m
Z = @(xi)  1i*sqrt(pi)*faddeeva(xi,16);

% Z0(ξ)的实现
Z0 = @(xi, kpar) Z(xi) .* (real(kpar) > 0) - Z(-xi) .* (real(kpar) < 0);

% 谐波数
n_e = 1;                              % 电子谐波数
n_i = 1;                              % 离子谐波数
V_e = 0;                              % 电子漂移速度
V_i = 0;                              % 离子漂移速度

% 完整色散方程（包含电子和离子）- 归一化版本
function F = my_disp_eq(kpar, omega, omega_pe, omega_pi, c, A1e, A1i, xi_n_e, xi_n_i, Z0)
    Ae = A1e(kpar);
    Ai = A1i(kpar);
    xi_e = xi_n_e(kpar);
    xi_i = xi_n_i(kpar);
    Z_e = Z0(xi_e, kpar);
    Z_i = Z0(xi_i, kpar);
    F = (kpar.^2 * c^2 - omega^2 - omega_pe^2 * Ae - omega_pi^2 * Ai) / c^2;
end

% 主扫描循环
for i = 1:num_n
    ne = density_values(i);
    ni = ne; % 假设准中性
    
    % 计算当前密度下的参数
    omega_pe = sqrt(ne*e^2/(eps0*me));    % 电子等离子体频率
    omega_pi = sqrt(ni*e^2/(eps0*mi));    % 离子等离子体频率
    omega_ce = e*B0/me;                   % 电子回旋频率
    omega_ci = e*B0/mi;                   % 离子回旋频率
    
    % 热速 (K -> m/s)
    w_para_e = sqrt(2*kB*Te/me);          % 电子平行热速
    w_para_i = sqrt(2*kB*Ti/mi);          % 离子平行热速
    
    % xi_n^l(kpar)的实现
    xi_n_e = @(kpar) (omega - kpar*V_e - n_e*omega_ce) ./ (kpar*w_para_e);
    xi_n_i = @(kpar) (omega - kpar*V_i - n_i*omega_ci) ./ (kpar*w_para_i);
    
    % A_{+1}^l(kpar)的实现
    A1e = @(kpar) (Tperp_e-Tpara_e)/(omega*Tpara_e) + ...
        ( (xi_n_e(kpar)*Tperp_e)/(omega*Tpara_e) + n_e*omega_ce./(omega*kpar*w_para_e) ) .* Z0(xi_n_e(kpar), kpar);
    
    A1i = @(kpar) (Tperp_i-Tpara_i)/(omega*Tpara_i) + ...
        ( (xi_n_i(kpar)*Tperp_i)/(omega*Tpara_i) + n_i*omega_ci./(omega*kpar*w_para_i) ) .* Z0(xi_n_i(kpar), kpar);
    
    % 初始猜测
    k0 = (0.1 + 0.001i) * omega / c;
    
    % 尝试求解
    try
        options = optimset('Display','off','TolFun',1e-10,'TolX',1e-10);
        kpar_sol = fsolve(@(kpar) my_disp_eq(kpar, omega, omega_pe, omega_pi, c, A1e, A1i, xi_n_e, xi_n_i, Z0), k0, options);
        
        % 存储结果
        k_real_results(i) = real(kpar_sol);
        k_imag_results(i) = imag(kpar_sol);
        omega_ci_results(i) = omega_ci;
        omega_ce_results(i) = omega_ce;
        success_flag(i) = 1;
        
        fprintf('  成功: n = %.2e, k = %.6e + %.6ei 1/m\n', ne, real(kpar_sol), imag(kpar_sol));
        
    catch ME
        fprintf('  失败: n = %.2e, %s\n', ne, ME.message);
        k_real_results(i) = NaN;
        k_imag_results(i) = NaN;
        omega_ci_results(i) = omega_ci;
        omega_ce_results(i) = omega_ce;
        success_flag(i) = 0;
    end
end

% 输出结果
fprintf('\n=== 扫描结果 ===\n');
fprintf('n(m^-3)\t\tomega_ci(rad/s)\t\tomega/omega_ci\t\tk_real(1/m)\t\tk_imag(1/m)\t\t状态\n');
fprintf('--------\t\t-------------\t\t-------------\t\t----------\t\t----------\t\t----\n');

for i = 1:num_n
    if success_flag(i)
        fprintf('%.2e\t%.2e\t%.3f\t%.6e\t%.6e\t成功\n', ...
            density_values(i), omega_ci_results(i), omega/omega_ci_results(i), ...
            k_real_results(i), k_imag_results(i));
    else
        fprintf('%.2e\t%.2e\t%.3f\t%s\t%s\t失败\n', ...
            density_values(i), omega_ci_results(i), omega/omega_ci_results(i), ...
            'NaN', 'NaN');
    end
end

% 绘制结果
figure('Position', [100, 100, 1200, 800]);

% 子图1: k的实部
subplot(2, 3, 1);
plot(density_values, k_real_results, 'b-*', 'LineWidth', 1, 'MarkerSize', 2)
xlabel('等离子体密度 n (m^{-3})');
ylabel('k_{real} (1/m)');
title('k的实部 vs 密度');
grid on;

% 子图2: k的虚部
subplot(2, 3, 2);
plot(density_values, k_imag_results, 'b-*', 'LineWidth', 1, 'MarkerSize', 2);
xlabel('等离子体密度 n (m^{-3})');
ylabel('k_{imag} (1/m)');
title('k的虚部 vs 密度');
grid on;

% 子图3: omega/omega_ci比值
subplot(2, 3, 3);
omega_ci_ratio = omega ./ omega_ci_results;
plot(density_values, omega_ci_ratio, 'b-*', 'LineWidth', 1, 'MarkerSize', 2)
xlabel('等离子体密度 n (m^{-3})');
ylabel('\omega/\omega_{ci}');
title('\omega/\omega_{ci} vs 密度');
grid on;
yline(1, '--k', 'LineWidth', 1); % 标记共振线

% 子图4: 归一化的k实部
subplot(2, 3, 4);
k_real_norm = c * k_real_results ./ omega_ci_results;
plot(density_values, k_real_norm, 'b-*', 'LineWidth', 1, 'MarkerSize', 2)
xlabel('等离子体密度 n (m^{-3})');
ylabel('c*k_{real}/\omega_{ci}');
title('归一化k实部 vs 密度');
grid on;

% 子图5: 归一化的k虚部
subplot(2, 3, 5);
k_imag_norm = c * k_imag_results ./ omega_ci_results;
plot(density_values, k_imag_norm, 'b-*', 'LineWidth', 1, 'MarkerSize', 2)
xlabel('等离子体密度 n (m^{-3})');
ylabel('c*k_{imag}/\omega_{ci}');
title('归一化k虚部 vs 密度');
grid on;

% 子图6: 成功/失败状态
subplot(2, 3, 6);
bar(density_values, success_flag, 'FaceColor', 'b', 'EdgeColor', 'k');
xlabel('等离子体密度 n (m^{-3})');
ylabel('求解状态');
title('求解成功/失败状态');
ylim([0, 1.2]);
yticks([0, 1]);
yticklabels({'失败', '成功'});
grid on;

% 保存结果到文件
results_table = table(density_values', omega_ci_results', omega_ce_results', ...
    omega./omega_ci_results', k_real_results', k_imag_results', ...
    success_flag', ...
    'VariableNames', {'n_m3', 'omega_ci_rad_s', 'omega_ce_rad_s', ...
    'omega_omega_ci_ratio', 'k_real_1_m', 'k_imag_1_m', 'success'});

writetable(results_table, 'k_solver_scan_n_results.csv');
fprintf('\n结果已保存到 k_solver_scan_n_results.csv\n');

% 分析共振区域
fprintf('\n=== 共振分析 ===\n');
resonance_idx = find(abs(omega_ci_ratio - 1) < 0.1);
if ~isempty(resonance_idx)
    fprintf('在以下密度值附近接近共振 (|omega/omega_ci - 1| < 0.1):\n');
    for i = 1:length(resonance_idx)
        idx = resonance_idx(i);
        fprintf('  n = %.2e m^-3, omega/omega_ci = %.3f\n', ...
            density_values(idx), omega_ci_ratio(idx));
    end
else
    fprintf('在扫描范围内未发现接近共振的情况\n');
end

% faddeeva函数实现（移植自faddeeva.m）
function w = faddeeva(z,N)
if nargin<2, N = []; end
if isempty(N), N = 16; end
w = zeros(size(z)); % initialize output
% for purely imaginary-valued inputs, use erf as is if z is real
idx = real(z)==0; %
w(idx) = exp(-z(idx).^2).*erfc(imag(z(idx)));
if all(idx), return; end
idx = ~idx;
% for complex-valued inputs
% make sure all points are in the upper half-plane (positive imag. values)
idx1 = idx & imag(z)<0;
z(idx1) = conj(z(idx1));
M = 2*N;
M2 = 2*M;
k = (-M+1:1:M-1)'; % 注意这里的 ' 不能去掉
L = sqrt(N/sqrt(2)); % Optimal choice of L.
theta = k*pi/M;
t = L*tan(theta/2); % Variables theta and t.
f = exp(-t.^2).*(L^2+t.^2);
f = [0; f]; % Function to be transformed.
a = real(fft(fftshift(f)))/M2; % Coefficients of transform.
a = flipud(a(2:N+1)); % Reorder coefficients.
Z = (L+1i*z(idx))./(L-1i*z(idx));
p = polyval(a,Z); % Polynomial evaluation.
w(idx) = 2*p./(L-1i*z(idx)).^2 + (1/sqrt(pi))./(L-1i*z(idx)); % Evaluate w(z).
% convert the upper half-plane results to the lower half-plane if necesary
w(idx1) = conj(2*exp(-z(idx1).^2) - w(idx1));
end
