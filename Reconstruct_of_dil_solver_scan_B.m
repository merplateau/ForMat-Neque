% 基于dil_solver_m.m的磁场扫描程序
% 扫描磁场强度，比较冷等离子体和热等离子体介电张量分量
clear; clc;


% 设置Python路径（根据您的Anaconda安装位置调整）
python_path = 'D:\ProgramData\anaconda3\envs\plasma\python.exe';


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
f = 190e3;               % 驱动频率 (Hz) - 固定s为190kHz
Te_eV = 0.0003;               % 电子温度 (eV)
Ti_eV = 0.00003;             % 离子温度 (eV)
Te = Te_eV * e / kB;     % 电子温度 (K)
Ti = Ti_eV * e / kB;     % 离子温度 (K)
Tperp_e = 1*Te;        % 电子垂直温度 (K)
Tpara_e = 1*Te;          % 电子平行温度 (K)
Tperp_i = 1*Ti;        % 离子垂直温度 (K)
Tpara_i = 1*Ti;          % 离子平行温度 (K)

% 固定频率
omega = 2*pi*f;          % 角频率 (rad/s)

% 扫描参数
B_values = 0.3:0.001:0.6; % 磁场强度范围 (T)
num_B = length(B_values);

% 存储结果
k_real_results = zeros(1, num_B);
k_imag_results = zeros(1, num_B);
K_perp_hot = zeros(1, num_B);
K_phi_hot = zeros(1, num_B);
K_par_hot = zeros(1, num_B);
K_perp_cold = zeros(1, num_B);
K_phi_cold = zeros(1, num_B);
K_par_cold = zeros(1, num_B);
omega_ci_results = zeros(1, num_B);
omega_ce_results = zeros(1, num_B);
success_flag = zeros(1, num_B);
k_real_cold_results = zeros(1, num_B);
k_imag_cold_results = zeros(1, num_B);

fprintf('开始磁场强度扫描...\n');
fprintf('频率: %.1f kHz\n', f/1e3);
fprintf('磁场范围: %.1f - %.1f T\n', min(B_values), max(B_values));
fprintf('扫描点数: %d\n\n', num_B);

% 色散函数Z(ξ)的实现 - 使用faddeeva.m
%Z = @(xi)  1i*sqrt(pi)*faddeeva(xi,16);
%Z = @(xi) 1i*sqrt(pi)*wTrap(xi,16);
%Z = @(xi) zetaf(xi);
%Z = @(xi) funZl(xi);


% 修改Z函数，添加NaN检查
%Z = @(xi) safe_Z_wrapper(xi);

%function result = safe_Z_wrapper(xi)
%    result = funZl(xi);
%    % 如果结果是NaN或Inf，返回0
%    if any(isnan(result)) || any(isinf(result))
%        result = 0;
%    end
%end

%Z = @(xi) 1i * sqrt(pi) *scipy_special.wofz(xi);


% Z0(ξ)的实现
%Z0 = @(xi, kpar) Z(xi) .* (real(kpar) > 0) - Z(-xi) .* (real(kpar) < 0);
Z0 = @(xi, kpar) Z(real(xi) + abs(imag(xi))*1i);

% 谐波数
n_e = 1;                              % 电子谐波数
n_i = 1;                              % 离子谐波数
V_e = 10000;                              % 电子漂移速度
V_i = 10000;                              % 离子漂移速度

% xi_n^l(kpar)的实现
xi_n = @(omega, kpar, V_l, n_l, omega_cl, w_par_l) (omega - kpar*V_l - n_l*omega_cl) ./ (kpar*w_par_l);

% A_l^n(kpar)的实现
A_l_n = @(Tperp_l, Tpara_l, omega, kpar, V_l, n_l, omega_cl, w_par_l) ...
    (Tperp_l-Tpara_l)/(omega*Tpara_l) + ...
    ( (xi_n(omega, kpar, V_l, n_l, omega_cl, w_par_l)*Tperp_l)/(omega*Tpara_l) + ...
      n_l*omega_cl./(kpar*omega*w_par_l) ) .* Z0(xi_n(omega, kpar, V_l, n_l, omega_cl, w_par_l), kpar);

% B_l^0(kpar)的实现
B_l_0 = @(Tperp_l, Tpara_l, omega, omega_cl, kpar, V_l, w_par_l) ...
    (omega*Tperp_l - kpar*V_l*Tpara_l)/(omega*kpar*Tpara_l) + ...
    (xi_n(omega, kpar, V_l, 0, omega_cl, w_par_l)*Tperp_l)./(kpar*Tpara_l) .* ...
    Z0(xi_n(omega, kpar, V_l, 0, omega_cl, w_par_l), kpar);

% 色散方程
function F = my_disp_eq(kpar, omega, omega_pe, omega_pi, c, A1e, A1i, Z0, V_e, V_i, n_e, n_i, omega_ce, omega_ci, w_par_e, w_par_i)
    xi_i = (omega - kpar*V_i - n_i*omega_ci) ./ (kpar*w_par_i);
    xi_e = (omega - kpar*V_e - n_e*omega_ce) ./ (kpar*w_par_e);
    Z0_i = Z0(xi_i,kpar);
    Z0_e = Z0(xi_e,kpar);
    F = (kpar^2 * c^2 - omega^2 - omega_pe^2 * omega * A1e(kpar) - omega_pi^2 * omega * A1i(kpar));
    fprintf('xi_r_i = %.6e, xi_i_i = %.6e, xi_r_e = %.6e, xi_i_e = %.6e \n', real(xi_i), imag(xi_i), real(xi_e), imag(xi_e));
    %fprintf('A1e = %.6e, A1i = %.6e \n', A1e(kpar), A1i(kpar));
    fprintf('Zr(xi_i) = %.6e, Zi(xi_i) = %.6e, Zr(xi_e) = %.6e, Zi(xi_e) = %.6e \n', real(Z0_i), imag(Z0_i), real(Z0_e), imag(Z0_e));
end

% 主扫描循环
for i = 1:num_B
    B0 = B_values(i);
    
    fprintf('处理 B = %.2f T (%d/%d)...\n', B0, i, num_B);
    
    % 计算当前磁场下的参数
    omega_pe = sqrt(ne*e^2/(eps0*me));    % 电子等离子体频率
    omega_pi = sqrt(ni*e^2/(eps0*mi));    % 离子等离子体频率
    omega_ce = e*B0/me;                   % 电子回旋频率
    omega_ci = e*B0/mi;                   % 离子回旋频率
    
    % 热速 (K -> m/s)
    w_para_e = sqrt(2*kB*Te/me);          % 电子平行热速
    w_para_i = sqrt(2*kB*Ti/mi);          % 离子平行热速
    w_perp_e = sqrt(2*kB*Tperp_e/me);     % 电子垂直热速
    w_perp_i = sqrt(2*kB*Tperp_i/mi);     % 离子垂直热速
    
    % 存储频率
    omega_ci_results(i) = omega_ci;
    omega_ce_results(i) = omega_ce;
    
    % 计算冷等离子体介电常数分量
    K_perp_cold(i) = 1 - (omega_pe^2/(omega^2 - omega_ce^2)) - (omega_pi^2/(omega^2 - omega_ci^2));
    K_phi_cold(i) = (omega_ce*omega_pe^2/(omega^3 - omega*omega_ce^2)) + (omega_ci*omega_pi^2/(omega^3 - omega*omega_ci^2));
    K_par_cold(i) = 1 - (omega_pe^2/omega^2) - (omega_pi^2/omega^2);
    
    % 先解冷等离子体色散方程，得到k_par_cold
    cold_disp_eq = @(kpar) kpar.^2 * c^2 - (omega^2 - ...
        omega_pe^2 * (omega - kpar*V_e) / (omega - kpar*V_e - omega_ce) ...
        - omega_pi^2 * (omega - kpar*V_i) / (omega - kpar*V_i - omega_ci));
    k0_cold = 10+0.1i;
    try
        options_cold = optimset('Display','off','TolFun',1e-10,'TolX',1e-10);
        kpar_cold_sol = fsolve(cold_disp_eq, k0_cold, options_cold);
        k_real_cold_results(i) = real(kpar_cold_sol);
        k_imag_cold_results(i) = imag(kpar_cold_sol);
    catch
        k_real_cold_results(i) = NaN;
        k_imag_cold_results(i) = NaN;
    end
    
    % 尝试求解热等离子体色散方程
    try
        % A1e和A1i函数
        A1e = @(kpar) A_l_n(Tperp_e, Tpara_e, omega, kpar, V_e, n_e, omega_ce, w_para_e);
        A1i = @(kpar) A_l_n(Tperp_i, Tpara_i, omega, kpar, V_i, n_i, omega_ci, w_para_i)

   
        
        % 初始猜测n_{\|}^2=1-\sum_s \frac{\omega_{p s}^2}{\omega^2} \frac{\omega-k_{\|} V_s}{\omega-k_{\|} V_s-\Omega_s}
        k0 = kpar_cold_sol;
        %k0 = -1 + 1i;
        
        
        % 求解kpar
        options = optimset('Display','iter','TolFun',1e-30,'TolX',1e-20);
        kpar_sol = fsolve(@(kpar) my_disp_eq(kpar, omega, omega_pe, omega_pi, c, A1e, A1i, Z0, V_e, V_i, n_e, n_i, omega_ce, omega_ci, w_para_e, w_para_i), k0, options);
        
        % 存储k结果
        k_real_results(i) = real(kpar_sol);
        k_imag_results(i) = imag(kpar_sol);
        success_flag(i) = 1;
        
        % 计算热等离子体介电张量分量
        % 电子参数
        params_e = struct('omega_p', omega_pe, 'omega_c', omega_ce, 'Tperp', Tperp_e, 'Tpara', Tpara_e, ...
            'V', V_e, 'n', n_e, 'w_par', w_para_e, 'w_perp', w_perp_e);
        % 离子参数
        params_i = struct('omega_p', omega_pi, 'omega_c', omega_ci, 'Tperp', Tperp_i, 'Tpara', Tpara_i, ...
            'V', V_i, 'n', n_i, 'w_par', w_para_i, 'w_perp', w_perp_i);
        
        % 初始化
        K_perp = 1;
        K_phi = 0;
        K_par = 1;
        
        for l = [params_e, params_i]
            % A_l^+1, A_l^-1
            A_p1 = A_l_n(l.Tperp, l.Tpara, omega, kpar_sol, l.V, +1, l.omega_c, l.w_par);
            A_m1 = A_l_n(l.Tperp, l.Tpara, omega, kpar_sol, l.V, -1, l.omega_c, l.w_par);
            % B_l^0
            B_0 = B_l_0(l.Tperp, l.Tpara, omega, l.omega_c, kpar_sol, l.V, l.w_par);
            % K_perp
            K_perp = K_perp + (l.omega_p^2/(2*omega)) * (A_m1 + A_p1);
            % K_phi
            K_phi = K_phi + (l.omega_p^2/(2*omega)) * (A_m1 - A_p1);
            % K_par
            K_par = K_par + (2*l.omega_p^2/(kpar_sol*l.w_perp^2)) * (l.V/omega + B_0);
        end
        
        % 存储热等离子体结果
        K_perp_hot(i) = K_perp;
        K_phi_hot(i) = K_phi;
        K_par_hot(i) = K_par;
        
        fprintf('  成功: k = %.6e + %.6ei 1/m\n', real(kpar_sol), imag(kpar_sol));
        
    catch ME
        fprintf('  失败: %s\n', ME.message);
        k_real_results(i) = NaN;
        k_imag_results(i) = NaN;
        K_perp_hot(i) = NaN;
        K_phi_hot(i) = NaN;
        K_par_hot(i) = NaN;
        success_flag(i) = 0;
    end
end

% 计算共振磁场（omega_ci = omega）
B_resonance = omega * mi / e;
fprintf('\n共振磁场: B = %.3f T (omega_ci = omega)\n', B_resonance);

% 输出结果
fprintf('\n=== 扫描结果 ===\n');
fprintf('B(T)\t\tomega_ci(rad/s)\t\tomega/omega_ci\t\tk_real(1/m)\t\tk_imag(1/m)\t\t状态\n');
fprintf('----\t\t-------------\t\t-------------\t\t----------\t\t----------\t\t----\n');

for i = 1:num_B
    if success_flag(i)
        fprintf('%.2f\t\t%.2e\t\t%.3f\t\t%.6e\t\t%.6e\t\t成功\n', ...
            B_values(i), omega_ci_results(i), omega/omega_ci_results(i), ...
            k_real_results(i), k_imag_results(i));
    else
        fprintf('%.2f\t\t%.2e\t\t%.3f\t\t%s\t\t%s\t\t失败\n', ...
            B_values(i), omega_ci_results(i), omega/omega_ci_results(i), ...
            'NaN', 'NaN');
    end
end

% 绘制结果
figure('Position', [100, 100, 1400, 1000]);

% 子图1: K_par 比较
subplot(3, 3, 1);
hold on;
%plot(B_values, imag(K_par_hot), 'r--*', 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', '热等离子体(虚部)');
plot(B_values, K_par_cold, 'b-*', 'LineWidth', 1.0, 'MarkerSize', 0.1, 'DisplayName', '冷等离子体');
plot(B_values, real(K_par_hot), 'r-.', 'LineWidth', 0.1, 'MarkerSize', 10, 'DisplayName', '热等离子体(实部)');
xline(B_resonance, '--k', 'LineWidth', 2, 'DisplayName', '共振磁场');
xlabel('磁场强度 B (T)');
ylabel('K_{par}');
title('K_{par} vs 磁场强度');
legend('Location', 'best');
grid on;

% 子图2: K_phi 比较
subplot(3, 3, 2);
hold on;
%plot(B_values, imag(K_phi_hot), 'r--*', 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', '热等离子体(虚部)');
plot(B_values, K_phi_cold, 'b-*', 'LineWidth', 1.0, 'MarkerSize', 0.1, 'DisplayName', '冷等离子体');
plot(B_values, real(K_phi_hot), 'r-.', 'LineWidth', 0.1, 'MarkerSize', 10, 'DisplayName', '热等离子体(实部)');
xline(B_resonance, '--k', 'LineWidth', 2, 'DisplayName', '共振磁场');
xlabel('磁场强度 B (T)');
ylabel('K_{phi}');
title('K_{phi} vs 磁场强度');
legend('Location', 'best');
grid on;

% 子图3: K_perp 比较
subplot(3, 3, 3);
hold on;
%plot(B_values, imag(K_perp_hot), 'r--*', 'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', '热等离子体(虚部)');
plot(B_values, K_perp_cold, 'b-*', 'LineWidth', 1.0, 'MarkerSize', 0.1, 'DisplayName', '冷等离子体');
plot(B_values, real(K_perp_hot), 'r-.', 'LineWidth', 0.1, 'MarkerSize', 10, 'DisplayName', '热等离子体(实部)');
xline(B_resonance, '--k', 'LineWidth', 2, 'DisplayName', '共振磁场');
xlabel('磁场强度 B (T)');
ylabel('K_{perp}');
title('K_{perp} vs 磁场强度');
legend('Location', 'best');
grid on;

% 子图1: K_par 比较
subplot(3, 3, 4);
hold on;
plot(B_values, imag(K_par_hot), 'r-.', 'LineWidth', 0.1, 'MarkerSize', 10, 'DisplayName', '热等离子体(虚部)');
%plot(B_values, K_par_cold, 'b-*', 'LineWidth', 1.0, 'MarkerSize', 0.1, 'DisplayName', '冷等离子体');
%plot(B_values, real(K_par_hot), 'r-.', 'LineWidth', 0.1, 'MarkerSize', 10, 'DisplayName', '热等离子体(实部)');
xline(B_resonance, '--k', 'LineWidth', 2, 'DisplayName', '共振磁场');
xlabel('磁场强度 B (T)');
ylabel('K_{par}');
title('K_{par} vs 磁场强度');
legend('Location', 'best');
grid on;

% 子图4: k的实部
subplot(3, 3, 5);
plot(B_values, k_real_results, 'r:.', 'LineWidth', 0.1, 'MarkerSize', 10, 'DisplayName', '热等离子体');
hold on;
plot(B_values, k_real_cold_results, 'b:.', 'LineWidth', 1.0, 'MarkerSize', 0.1, 'DisplayName', '冷等离子体');
xline(B_resonance, '--k', 'LineWidth', 2);
xlabel('磁场强度 B (T)');
ylabel('k_{real} (1/m)');
title('k的实部 vs 磁场强度');
legend('Location', 'best');
grid on;

% 子图5: k的虚部
subplot(3, 3, 6);
plot(B_values, k_imag_results, 'r:.', 'LineWidth', 0.1, 'MarkerSize', 10, 'DisplayName', '热等离子体');
hold on;
plot(B_values, k_imag_cold_results, 'b:.', 'LineWidth', 1.0, 'MarkerSize', 0.1, 'DisplayName', '冷等离子体');
xline(B_resonance, '--k', 'LineWidth', 2);
xlabel('磁场强度 B (T)');
ylabel('k_{imag} (1/m)');
title('k的虚部 vs 磁场强度');
legend('Location', 'best');
grid on;

% 子图6: omega/omega_ci比值
subplot(3, 3, 7);
omega_ci_ratio = omega ./ omega_ci_results;
plot(B_values, omega_ci_ratio, 'b-*', 'LineWidth', 1.5, 'MarkerSize', 3);
xline(B_resonance, '--k', 'LineWidth', 2);
yline(1, '--r', 'LineWidth', 1); % 标记共振线
xlabel('磁场强度 B (T)');
ylabel('\omega/\omega_{ci}');
title('\omega/\omega_{ci} vs 磁场强度');
grid on;

% 子图7: 归一化的k实部
subplot(3, 3, 8);
k_real_norm = c * k_real_results ./ omega_ci_results;
plot(B_values, k_real_norm, 'b-*', 'LineWidth', 1.5, 'MarkerSize', 3);
xline(B_resonance, '--k', 'LineWidth', 2);
xlabel('磁场强度 B (T)');
ylabel('c*k_{real}/\omega_{ci}');
title('归一化k实部 vs 磁场强度');
grid on;

% 子图8: 归一化的k虚部
subplot(3, 3, 9);
k_imag_norm = c * k_imag_results ./ omega_ci_results;
plot(B_values, k_imag_norm, 'b-*', 'LineWidth', 1.5, 'MarkerSize', 3);
xline(B_resonance, '--k', 'LineWidth', 2);
xlabel('磁场强度 B (T)');
ylabel('c*k_{imag}/\omega_{ci}');
title('归一化k虚部 vs 磁场强度');
grid on;

% 子图9: 成功/失败状态
%subplot(3, 3, 9);
%bar(B_values, success_flag, 'FaceColor', 'b', 'EdgeColor', 'k');
%xline(B_resonance, '--k', 'LineWidth', 2);
%xlabel('磁场强度 B (T)');
%ylabel('求解状态');
%title('求解成功/失败状态');
%ylim([0, 1.2]);
%yticks([0, 1]);
%yticklabels({'失败', '成功'});
%grid on;

% 保存结果到文件
results_table = table(B_values', omega_ci_results', omega_ce_results', ...
    omega./omega_ci_results', k_real_results', k_imag_results', ...
    real(K_perp_hot)', imag(K_perp_hot)', real(K_phi_hot)', imag(K_phi_hot)', ...
    real(K_par_hot)', imag(K_par_hot)', K_perp_cold', K_phi_cold', K_par_cold', ...
    success_flag', ...
    'VariableNames', {'B_T', 'omega_ci_rad_s', 'omega_ce_rad_s', ...
    'omega_omega_ci_ratio', 'k_real_1_m', 'k_imag_1_m', ...
    'K_perp_hot_real', 'K_perp_hot_imag', 'K_phi_hot_real', 'K_phi_hot_imag', ...
    'K_par_hot_real', 'K_par_hot_imag', 'K_perp_cold', 'K_phi_cold', 'K_par_cold', ...
    'success'});

writetable(results_table, 'dil_solver_scan_B_results.csv');
fprintf('\n结果已保存到 dil_solver_scan_B_results.csv\n');

% 分析共振区域
fprintf('\n=== 共振分析 ===\n');
resonance_idx = find(abs(omega_ci_ratio - 1) < 0.1);
if ~isempty(resonance_idx)
    fprintf('在以下磁场值附近接近共振 (|omega/omega_ci - 1| < 0.1):\n');
    for i = 1:length(resonance_idx)
        idx = resonance_idx(i);
        fprintf('  B = %.2f T, omega/omega_ci = %.3f\n', ...
            B_values(idx), omega_ci_ratio(idx));
    end
else
    fprintf('在扫描范围内未发现接近共振的情况\n');
end

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

function w=zetaf(z)
%--------------------------------------------------
%               2015-04-20 Tsing Wong 
%               2015-05-12 update
% compute the faddeeva function:
%       w(z)=exp(-z^2)erfc(-iz);
% Reference: 
%        [1] M.Zaghloul&A. Ali. "Algorithm 916:
%            computing the Faddeyeva and Voigt func-
%            tions". ACM Trans. Math. Soft. 38 (2), 
%            15 (2011). Preprint available at arXiv:
%            1106.0151
% --------------------------------------------------
x=real(z);y=imag(z);
a=0.5;% a-convergence parameter
x1=x(:);
nz=20+max(max(ceil(abs(x1/a))));% n- truncation parameter
%n=20+ceil(abs(x/a));
s1=0;s2=0;s3=0;s4=0;s5=0;
for k=1:nz
    %sk=1./(a^2*k^2+y.*y).*exp(-(a^2*k^2+x.*x));
    s1=s1+1./(a^2*k^2+y.*y).*exp(-(a^2*k^2+x.*x));
    s2=s2+1./(a^2*k^2+y.*y).*exp(-(a*k+x).^2);
    s3=s3+1./(a^2*k^2+y.*y).*exp(-(a*k-x).^2);
    s4=s4+a*k./(a^2*k^2+y.*y).*exp(-(a*k+x).^2);
    s5=s5+a*k./(a^2*k^2+y.*y).*exp(-(a*k-x).^2);
end
sin1=sin(x.*y)./(x.*y);
sin1(x.*y==0)=1;% there will be a bug (1 or -1)
sin2=sin(2*x.*y)./(2*x.*y);
sin2(x.*y==0)=1;% there will be a bug (1 or -1)
% caculate the real part of w(z)
Rw1=exp(-x.*x).*erfcx(y).*cos(2*x.*y);
Rw2=2*a*x.*sin(x.*y).*exp(-x.*x).*sin1/pi;
Rw3=(2*a/pi)*(-y.*cos(2*x.*y).*s1+y.*s2/2+y.*s3/2);
Rw=Rw1+Rw2+Rw3;
% caculate the imaginary part of w(z)
Iw1=-exp(-x.*x).*erfcx(y).*sin(2*x.*y);
Iw2=(2*a*x/pi).*exp(-x.*x).*sin2;
Iw3=(2*a/pi)*(y.*sin(2*x.*y).*s1-s4/2+s5/2);
Iw=Iw1+Iw2+Iw3;
% caculate  w(z)
w=Rw+1i*Iw;
% the relationship between the plasma dispersion function 
% and the faddeeva function : zeta=i*sqrt(pi)w(z);
w=1i*sqrt(pi)*w;
end

function w = wTrapWCP(z,N)     
% Computes Faddeeva function w(z) for arbitrary complex z.
% See 'Computation of the complex error function using modified trapezoidal
% rules', M Al Azah and S N Chandler-Wilde, 2020, for details of method.
% Inputs: z complex (scalar, vector, or array).
%         N positive integer (N+1 is the number of quadrature points used);
%           the choice N = 11 is recommended for absolute errors (and 
%           relative errors in the upper half-plane) < 1.7frace-15  
% Output: w complex array of same dimensions as z containing 
%           estimates for values of w(z)
selectXneg = real(z) < 0; selectYneg = imag(z) < 0;
selectNotBoth = xor(selectXneg,selectYneg);
zYneg = z(selectYneg);
z(selectXneg) = -z(selectXneg); z(selectNotBoth) = conj(z(selectNotBoth));
w = wTrap(z,N);
w(selectNotBoth) = conj(w(selectNotBoth));
w(selectYneg) = 2*exp(-zYneg.*zYneg) - w(selectYneg);
end

function w = wTrap(z,N)     
% Computes Faddeeva function w(z) for z = x+iy with x,y >=0.
% See 'Computation of the complex error function using modified trapezoidal
% rules', M Al Azah and S N Chandler-Wilde, 2021, for details of method.
% Inputs: z complex (scalar, vector, or array).
%         N positive integer (N+1 is the number of quadrature points used);
%           the choice N = 11 is recommended for absolute and relative 
%           errors < 2e-15  
% Output: w complex array of same dimensions as z containing 
%           estimates for values of w(z)
h = sqrt(pi/(N+1)); 
H = pi/h;  
rz = real(z); rzh = rz/h; iz = imag(z); 
buff = abs(rzh-floor(rzh)-0.5);
select1 = iz >= max(rz,H); 
select2 = (iz < rz) & (buff <= 0.25);
select3 = ~(select1|select2);
w = zeros(size(z)); 
w(select1) = wM(z(select1),N);
w(select2) = wMT(z(select2),N);
w(select3) = wMM(z(select3),N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ih = wM(z,N) 
% Midpoint rule approximation
z2 = z.*z; az = (2i/H)*z;
t = h*((N:-1:1)+0.5); t2 = t.^2; et2 = exp(-t2); h0 = 0.5*h;
Sum = exp(-h0^2)./(z2-h0^2); 
for n = 1 : N
    Sum = Sum + et2(n)./(z2-t2(n));
end
Ih = az.*Sum;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ihstar = wMM(z,N)
% Modified midpoint rule approximation
z2 = z.*z; az = (2i/H)*z;
t = h*((N:-1:1)+0.5); t2 = t.^2; et2 = exp(-t2); h0 = 0.5*h;
Sum = exp(-h0.^2)./(z2-h0.^2); 
for n = 1 : N
    Sum = Sum + et2(n)./(z2-t2(n));
end
Ch = 2./(exp(z2).*(1+exp((-2i*H)*z))); % This is the modification
Ihstar = az.*Sum + Ch;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ihstar = wMT(z,N)   
% Modified trapezium rule approximation
z2 = z.*z; az = (2i/H)*z;
tau = h*(N:-1:1); tau2 = tau.^2; et2 = exp(-tau2);  
Sum = et2(1)./(z2-tau2(1));
for n = 2 : N
    Sum = Sum + et2(n)./(z2-tau2(n));
end
Ch2 = 2./(exp(z2).*(1-exp((-2i*H)*z))); % This is the modification
Ihstar = (1i/H)./z + az.*Sum + Ch2; 
end
end


function Z=funZl(z,l,N)
% Hua-sheng XIE, huashengxie@gmail.com, IFTS-ZJU, 2013-05-26 23:37
% Calculate GPDF, Generalized Plasma Dispersion Function, see [Xie2013]
%
% Modify the parameters in the code for your own usage
%
% Ref: [1] H. S. Xie, Generalized Plasma Dispersion Function:
% One-Solve-All Treatment, Visualizations and Application to
% Landau Damping, PoP, 2013. (Or http://arxiv.org/abs/1305.6476)
%
% Examples:
%  1. z=-1:0.01:1; [Zp,Z]=zetaph(z); plot(z,real(Zp),z,imag(Zp));
%  2. z=-1:0.01:1; F = '1./(1+v.^2)/pi'; [Zp,Z]=zetaph(z,0,F);
%     plot(z,real(Z),z,imag(Z));
%  3. [x,y]=meshgrid(-1:0.1:1,-1:0.1:1); z=x+1i*y; [Zp,Z]=zetaph(z,1);
%     subplot(121);surf(x,y,real(Zp)); subplot(122);surf(x,y,imag(Zp));
% 24-12-28 07:58 modify for v^l*exp(-v^2)/sqrt(pi)
if nargin<3, N = 128; end

% Default for usual Plasma Dispersion Function
if (nargin<2)
    l=0;
end
F = ['v.^',num2str(l),'.*exp(-v.^2)/sqrt(pi)'];
Z=calZ(z,F,N);
end

function Z=calZ(z,F,N,del)
if nargin<4, del = 1; end
Z=hilb(z,F,N,del);
ind1=find(isnan(Z));
z1=z(ind1)+1e-10;         % avoid NaN of higher order singular point
Z(ind1)=hilb(z1,F,N,del);% e.g., z=ia for Lorentzian F = '1./(a^2+v.^2)'
end

function Z=hilb(z,F,N,del,L)
% note: 1. in fact, a_n need calcualte only once for a fixed F, so you can
%   rewrite the code to speed up.
%       2. f(z) in analytic continuation can also be replaced by
%   sum{a_n*rho_n}
%       3. usually, del=1, but for flat-top and triangular distributions,
%   it seems we should set del=0 for correct analytic continuation

if nargin<5, L = sqrt(N/sqrt(2)); end  % optimal choice of L
%     L=10;
if nargin<4, del = 1; end

% 1. Define initial parameters
Z = zeros(size(z)); % initialize output

% 2. Calculate
idx=find(imag(z) == 0);
idx1=find(imag(z) ~= 0);
% 2.1 real line
for id=idx
    n  = [-N:N-1]';                   % Compute the collocation points
    v  = L*tan(pi*(n+1/2)/(2*N));
    FF = eval(F);             % sample the function
    FF(isnan(FF))=0;
    Fz = eval(['@(v)',F]);    % define F(z), for analytic continuation

    % Weideman95 method to calculate the real line Hilbert transform
    a  = fft(fftshift(FF.*(L-1i*v)));     % These three lines compute
    a  = exp(-1i*n*pi/(2*N)).*fftshift(a);   % expansion coefficients
    a  = flipud(1i*(sign(n+1/2).*a))/(2*N);
    t  = (L+1i*z(id))./(L-1i*z(id)); % The evaluation of the transform
    h  = polyval(a,t)./(t.^N.*(L-1i*z(id)));   % reduces to polynomial
    % evaluation in the variable t
    Z(id) = h + 1i.*Fz(z(id));
end
% 2.2. upper and lower half plane
for id=idx1
    M = 2*N;  M2 = 2*M;
    k = [-M+1:1:M-1]';% M2 = no. of sampling points
    theta = k*pi/M; v = L*tan(theta/2);     % define theta & v
    FF = eval(F);             % sample the function
    FF(isnan(FF))=0;
    Fz = eval(['@(v)',F]);    % define F(z), for analytic continuation

    % Weideman94 method to calculate upper half plane Hilbert transform
    W = (L^2+v.^2);                         % default weight function
    FF = FF.*W; FF = [0; FF];             % function to be transformed
    a = (fft(fftshift(FF)))/M2;           % coefficients of transform
    a0 = a(1); a = flipud(a(2:N+1));            % reorder coefficients
    z1 = (imag(z(id))>0).*z(id)+(imag(z(id))<0).*conj(z(id));
    t = (L+1i*z1)./(L-1i*z1); p = polyval(a,t); % polynomial evaluation
    h = 1i*(2*p./(L-1i*z1).^2+(a0/L)./(L-1i*z1));  % Evaluate h(f(v))
    % convert upper half-plane results to lower half-plane if necesary
    %         Z(id) = h.*(imag(z(id))>0)+conj(h-2i.*Fz(z1)).*(imag(z(id))<0);
    Z(id) = h.*(imag(z(id))>0)+(conj(h)+...
        del*2i.*Fz(z(id))).*(imag(z(id))<0);
end
Z=Z.*pi;
end
