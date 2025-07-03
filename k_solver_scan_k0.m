% 基于k_solver.m的参数扫描程序（扫描初值版本）
% B固定为0.5T，初值k0从(0.1)*omega/c到(10)*omega/c
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
f = 190e3;               % 驱动频率 (Hz) - 固定为190kHz
Te_eV = 3;               % 电子温度 (eV)
Ti_eV = 0.3;             % 离子温度 (eV)
Te = Te_eV * e / kB;     % 电子温度 (K)
Ti = Ti_eV * e / kB;     % 离子温度 (K)
Tperp_e = 10*Te;         % 电子垂直温度 (K)
Tpara_e = 1*Te;          % 电子平行温度 (K)
Tperp_i = 10*Ti;         % 离子垂直温度 (K)
Tpara_i = 1*Ti;          % 离子平行温度 (K)

% 固定频率
omega = 2*pi*f;          % 角频率 (rad/s)

% 固定磁场
B0 = 0.5;                % 磁场强度 (T)

% 扫描初值
k0_values = linspace(1, 1000, 100); % 50个初值点
num_k0 = length(k0_values);

% 存储结果
k_real_results = zeros(1, num_k0);
k_imag_results = zeros(1, num_k0);
success_flag = zeros(1, num_k0);

% 计算当前磁场下的参数
omega_pe = sqrt(ne*e^2/(eps0*me));    % 电子等离子体频率
omega_pi = sqrt(ni*e^2/(eps0*mi));    % 离子等离子体频率
omega_ce = e*B0/me;                   % 电子回旋频率
omega_ci = e*B0/mi;                   % 离子回旋频率

% 热速 (K -> m/s)
w_para_e = sqrt(2*kB*Te/me);          % 电子平行热速
w_para_i = sqrt(2*kB*Ti/mi);          % 离子平行热速

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

% xi_n^l(kpar)的实现
xi_n_e = @(kpar) (omega - kpar*V_e - n_e*omega_ce) ./ (kpar*w_para_e);
xi_n_i = @(kpar) (omega - kpar*V_i - n_i*omega_ci) ./ (kpar*w_para_i);

% A_{+1}^l(kpar)的实现
A1e = @(kpar) (Tperp_e-Tpara_e)/(omega*Tpara_e) + ...
    ( (xi_n_e(kpar)*Tperp_e)/(omega*Tpara_e) + n_e*omega_ce./(omega*kpar*w_para_e) ) .* Z0(xi_n_e(kpar), kpar);
A1i = @(kpar) (Tperp_i-Tpara_i)/(omega*Tpara_i) + ...
    ( (xi_n_i(kpar)*Tperp_i)/(omega*Tpara_i) + n_i*omega_ci./(omega*kpar*w_para_i) ) .* Z0(xi_n_i(kpar), kpar);

fprintf('B固定为%.2f T，初值k0从%.2f到%.2f (omega/c)\n', B0, k0_values(1)*c/omega, k0_values(end)*c/omega);

for i = 1:num_k0
    k0 = k0_values(i);
    fprintf('处理 k0 = %.2f * omega/c (%d/%d)...\n', k0*c/omega, i, num_k0);
    try
        options = optimset('Display','off','TolFun',1e-30,'TolX',1e-30);
        kpar_sol = fsolve(@(kpar) my_disp_eq(kpar, omega, omega_pe, omega_pi, c, A1e, A1i, xi_n_e, xi_n_i, Z0), k0, options);
        k_real_results(i) = real(kpar_sol);
        k_imag_results(i) = imag(kpar_sol);
        success_flag(i) = 1;
        fprintf('  成功: k = %.6e + %.6ei 1/m\n', real(kpar_sol), imag(kpar_sol));
    catch ME
        k_real_results(i) = NaN;
        k_imag_results(i) = NaN;
        success_flag(i) = 0;
        fprintf('  失败: %s\n', ME.message);
    end
end

% 输出结果
fprintf('\n=== 扫描结果 ===\n');
fprintf('k0(omega/c)\tk_real(1/m)\tk_imag(1/m)\t状态\n');
for i = 1:num_k0
    if success_flag(i)
        fprintf('%.2f\t%.6e\t%.6e\t成功\n', k0_values(i)*c/omega, k_real_results(i), k_imag_results(i));
    else
        fprintf('%.2f\tNaN\tNaN\t失败\n', k0_values(i)*c/omega);
    end
end

% 绘图
figure;
subplot(1,2,1);
plot(k0_values*c/omega, k_real_results, 'b-o', 'LineWidth', 1, 'MarkerSize', 3);
xlabel('初值 k_0 (\omega/c)'); ylabel('k_{real} (1/m)');
title('k_{real} vs 初值k_0'); grid on;

subplot(1,2,2);
plot(k0_values*c/omega, k_imag_results, 'r-o', 'LineWidth', 1, 'MarkerSize', 3);
xlabel('初值 k_0 (\omega/c)'); ylabel('k_{imag} (1/m)');
title('k_{imag} vs 初值k_0'); grid on;

% 保存结果到文件
results_table = table((k0_values'*c/omega), k_real_results', k_imag_results', success_flag', ...
    'VariableNames', {'k0_omega_c', 'k_real_1_m', 'k_imag_1_m', 'success'});
writetable(results_table, 'k_solver_scan_k0_results.csv');
fprintf('\n结果已保存到 k_solver_scan_k0_results.csv\n');

% faddeeva函数实现（移植自faddeeva.m）
function w = faddeeva(z,N)
if nargin<2, N = []; end
if isempty(N), N = 16; end
w = zeros(size(z)); % initialize output
idx = real(z)==0; %
w(idx) = exp(-z(idx).^2).*erfc(imag(z(idx)));
if all(idx), return; end
idx = ~idx;
idx1 = idx & imag(z)<0;
z(idx1) = conj(z(idx1));
M = 2*N;
M2 = 2*M;
k = (-M+1:1:M-1)';
L = sqrt(N/sqrt(2));
theta = k*pi/M;
t = L*tan(theta/2);
f = exp(-t.^2).*(L^2+t.^2);
f = [0; f];
a = real(fft(fftshift(f)))/M2;
a = flipud(a(2:N+1));
Z = (L+1i*z(idx))./(L-1i*z(idx));
p = polyval(a,Z);
w(idx) = 2*p./(L-1i*z(idx)).^2 + (1/sqrt(pi))./(L-1i*z(idx));
w(idx1) = conj(2*exp(-z(idx1).^2) - w(idx1));
end 