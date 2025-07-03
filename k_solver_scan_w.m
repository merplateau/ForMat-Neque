 % 基于k_solver.m的参数扫描程序
% 扫描频率omega，磁场强度固定
clear; clc;

% 基本常数
e = 1.602e-19;           % 元电荷 (C)
me = 9.11e-31;           % 电子质量 (kg)
mi = 40*1.67e-27;        % 质子质量 (kg)
eps0 = 8.854e-12;        % 真空介电常数 (F/m)
c = 3e8;                 % 光速 (m/s)
kB = 1.38e-23;           % 玻尔兹曼常数 (J/K)

% VASIMR ICRH 典型参数
ne = 1e18;               % 电子密度 (m^-3)
ni = ne;                 % 离子密度 (m^-3)，假设准中性
f0 = 190e3;              % 基准驱动频率 (Hz) - 190kHz
Te_eV = 3;               % 电子温度 (eV)
Ti_eV = 0.3;             % 离子温度 (eV)
Te = Te_eV * e / kB;     % 电子温度 (K)
Ti = Ti_eV * e / kB;     % 离子温度 (K)
Tperp_e = 10*Te;         % 电子垂直温度 (K)
Tpara_e = 1*Te;          % 电子平行温度 (K)
Tperp_i = 10*Ti;         % 离子垂直温度 (K)
Tpara_i = 1*Ti;          % 离子平行温度 (K)

% 固定磁场强度
B0 = 0.5;                % 固定磁场强度 (T)

% 扫描参数 - 在基准频率上下一个量级范围内扫描
f_min = f0 / 10;         % 最小频率 (19kHz)
f_max = f0 * 10;         % 最大频率 (1.9MHz)
f_values = logspace(log10(f_min), log10(f_max), 100); % 对数分布的频率点
num_f = length(f_values);

% 存储结果
k_real_results = zeros(1, num_f);
k_imag_results = zeros(1, num_f);
omega_ci_results = zeros(1, num_f);
omega_ce_results = zeros(1, num_f);
omega_results = zeros(1, num_f);
success_flag = zeros(1, num_f);

fprintf('开始频率扫描...\n');
fprintf('磁场强度: %.1f T\n', B0);
fprintf('频率范围: %.1f - %.1f kHz\n', f_values(1)/1e3, f_values(end)/1e3);
fprintf('扫描点数: %d\n\n', num_f);

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

% 计算固定磁场下的参数
omega_pe = sqrt(ne*e^2/(eps0*me));    % 电子等离子体频率
omega_pi = sqrt(ni*e^2/(eps0*mi));    % 离子等离子体频率
omega_ce = e*B0/me;                   % 电子回旋频率
omega_ci = e*B0/mi;                   % 离子回旋频率

fprintf('等离子体参数:\n');
fprintf('  电子等离子体频率: %.2e rad/s (%.1f MHz)\n', omega_pe, omega_pe/(2*pi*1e6));
fprintf('  离子等离子体频率: %.2e rad/s (%.1f MHz)\n', omega_pi, omega_pi/(2*pi*1e6));
fprintf('  电子回旋频率: %.2e rad/s (%.1f MHz)\n', omega_ce, omega_ce/(2*pi*1e6));
fprintf('  离子回旋频率: %.2e rad/s (%.1f kHz)\n', omega_ci, omega_ci/(2*pi*1e3));
fprintf('\n');

% 主扫描循环
for i = 1:num_f
    f = f_values(i);
    omega = 2*pi*f;      % 角频率 (rad/s)
    
    fprintf('处理 f = %.1f kHz (%.1f MHz) (%d/%d)...\n', f/1e3, f/1e6, i, num_f);
    
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
    
    % 初始猜测 - 根据当前频率调整
    k0 = (80000) * omega / c;
    
    % 尝试求解
    try
        options = optimset('Display','iter','TolFun',1e-10,'TolX',1e-15);
        kpar_sol = fsolve(@(kpar) my_disp_eq(kpar, omega, omega_pe, omega_pi, c, A1e, A1i, xi_n_e, xi_n_i, Z0), k0, options);
        
        % 存储结果
        k_real_results(i) = real(kpar_sol);
        k_imag_results(i) = imag(kpar_sol);
        omega_ci_results(i) = omega_ci;
        omega_ce_results(i) = omega_ce;
        omega_results(i) = omega;
        success_flag(i) = 1;
        
        fprintf('  成功: k = %.6e + %.6ei 1/m\n', real(kpar_sol), imag(kpar_sol));
        
    catch ME
        fprintf('  失败: %s\n', ME.message);
        k_real_results(i) = NaN;
        k_imag_results(i) = NaN;
        omega_ci_results(i) = omega_ci;
        omega_ce_results(i) = omega_ce;
        omega_results(i) = omega;
        success_flag(i) = 0;
    end
end

% 输出结果
fprintf('\n=== 扫描结果 ===\n');
fprintf('f(kHz)\t\tomega(rad/s)\t\tomega/omega_ci\t\tk_real(1/m)\t\tk_imag(1/m)\t\t状态\n');
fprintf('-----\t\t------------\t\t-------------\t\t----------\t\t----------\t\t----\n');

for i = 1:num_f
    if success_flag(i)
        fprintf('%.1f\t\t%.2e\t\t%.3f\t\t%.6e\t\t%.6e\t\t成功\n', ...
            f_values(i)/1e3, omega_results(i), omega_results(i)/omega_ci_results(i), ...
            k_real_results(i), k_imag_results(i));
    else
        fprintf('%.1f\t\t%.2e\t\t%.3f\t\t%s\t\t%s\t\t失败\n', ...
            f_values(i)/1e3, omega_results(i), omega_results(i)/omega_ci_results(i), ...
            'NaN', 'NaN');
    end
end

% 绘制结果
figure('Position', [100, 100, 1200, 800]);

% 子图1: k的实部
subplot(2, 3, 1);
semilogx(f_values/1e3, k_real_results, 'b-*', 'LineWidth', 1, 'MarkerSize', 2)
xlabel('频率 f (kHz)');
ylabel('k_{real} (1/m)');
title('k的实部 vs 频率');
grid on;

% 子图2: k的虚部
subplot(2, 3, 2);
semilogx(f_values/1e3, k_imag_results, 'b-*', 'LineWidth', 1, 'MarkerSize', 2);
xlabel('频率 f (kHz)');
ylabel('k_{imag} (1/m)');
title('k的虚部 vs 频率');
grid on;

% 子图3: omega/omega_ci比值
subplot(2, 3, 3);
omega_ci_ratio = omega_results ./ omega_ci_results;
semilogx(f_values/1e3, omega_ci_ratio, 'b-*', 'LineWidth', 1, 'MarkerSize', 2)
xlabel('频率 f (kHz)');
ylabel('\omega/\omega_{ci}');
title('\omega/\omega_{ci} vs 频率');
grid on;
yline(1, '--k', 'LineWidth', 1); % 标记共振线

% 子图4: 归一化的k实部
subplot(2, 3, 4);
k_real_norm = c * k_real_results ./ omega_results;
semilogx(f_values/1e3, k_real_norm, 'b-*', 'LineWidth', 1, 'MarkerSize', 2)
xlabel('频率 f (kHz)');
ylabel('c*k_{real}/\omega');
title('归一化k实部 vs 频率');
grid on;

% 子图5: 归一化的k虚部
subplot(2, 3, 5);
k_imag_norm = c * k_imag_results ./ omega_results;
semilogx(f_values/1e3, k_imag_norm, 'b-*', 'LineWidth', 1, 'MarkerSize', 2)
xlabel('频率 f (kHz)');
ylabel('c*k_{imag}/\omega');
title('归一化k虚部 vs 频率');
grid on;

% 子图6: 成功/失败状态
subplot(2, 3, 6);
semilogx(f_values/1e3, success_flag, 'b-*', 'LineWidth', 1, 'MarkerSize', 2);
xlabel('频率 f (kHz)');
ylabel('求解状态');
title('求解成功/失败状态');
ylim([0, 1.2]);
yticks([0, 1]);
yticklabels({'失败', '成功'});
grid on;

% 保存结果到文件
results_table = table(f_values', omega_results', omega_ci_results', omega_ce_results', ...
    omega_results'./omega_ci_results', k_real_results', k_imag_results', ...
    success_flag', ...
    'VariableNames', {'f_Hz', 'omega_rad_s', 'omega_ci_rad_s', 'omega_ce_rad_s', ...
    'omega_omega_ci_ratio', 'k_real_1_m', 'k_imag_1_m', 'success'});

writetable(results_table, 'k_solver_scan_w_results.csv');
fprintf('\n结果已保存到 k_solver_scan_w_results.csv\n');

% 分析共振区域
fprintf('\n=== 共振分析 ===\n');
resonance_idx = find(abs(omega_ci_ratio - 1) < 0.1);
if ~isempty(resonance_idx)
    fprintf('在以下频率附近接近共振 (|omega/omega_ci - 1| < 0.1):\n');
    for i = 1:length(resonance_idx)
        idx = resonance_idx(i);
        fprintf('  f = %.1f kHz, omega/omega_ci = %.3f\n', ...
            f_values(idx)/1e3, omega_ci_ratio(idx));
    end
else
    fprintf('在扫描范围内未发现接近共振的情况\n');
end

% 分析等离子体频率关系
fprintf('\n=== 等离子体频率分析 ===\n');
fprintf('电子等离子体频率: %.1f MHz\n', omega_pe/(2*pi*1e6));
fprintf('离子等离子体频率: %.1f MHz\n', omega_pi/(2*pi*1e6));
fprintf('电子回旋频率: %.1f MHz\n', omega_ce/(2*pi*1e6));
fprintf('离子回旋频率: %.1f kHz\n', omega_ci/(2*pi*1e3));

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