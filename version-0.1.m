clear; clc;


% 设置Python路径
python_path = 'D:\ProgramData\anaconda3\envs\plasma\python.exe';

    % 检查Python路径是否存在
    if ~exist(python_path, 'file')
        fprintf('Python路径不存在');
        return;
    end
    
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
n = 1e18;               % 电子密度 (m^-3)
ne = n;
ni = n;                 % 离子密度 (m^-3)，假设准中性
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
omega_c_i_results = zeros(1, num_B);
omega_c_e_results = zeros(1, num_B);
success_flag = zeros(1, num_B);
k_real_cold_results = zeros(1, num_B);
k_imag_cold_results = zeros(1, num_B);


% 谐波数
n_e = 1;                              % 电子谐波数
n_i = 1;                              % 离子谐波数
V_e = 10000;                              % 电子漂移速度
V_i = 10000;                              % 离子漂移速度

ele = struct('omega_p', omega_p_e, 'omega_c', omega_c_e, 'Tperp', Tperp_e, 'Tpara', Tpara_e, 'V', V_e, 'w_par', w_para_e, 'w_perp', w_perp_e, 'm', me);
ion = struct('omega_p', omega_p_i, 'omega_c', omega_c_i, 'Tperp', Tperp_i, 'Tpara', Tpara_i, 'V', V_i, 'w_par', w_para_i, 'w_perp', w_perp_i, 'm', mi);

function result = hotPDE(kpar, omega, ele, ion)
    result = kpar^2 * c^2 - omega^2 - ele.omega_p^2 * omega * funA(omega, kpar, 1, ele) - ion.omega_p^2 * omega * funA(omega, kpar, 1, ion);
end

function result = coldPDE(kpar, omega, ele, ion)
    result = kpar.^2 * c^2 - (omega^2 - ...
        ele.omega_p^2 * (omega - kpar*ele.V) / (omega - kpar*ele.V - ele.omega_c) ...
        - ion.omega_p^2 * (omega - kpar*ion.V) / (omega - kpar*ion.V - ion.omega_c));
end

function result = Z0(xi, kpar)
    z = real(xi) + abs(imag(xi))*1i;
    result = 1i*sqrt(pi)*scipy_special.wofz(z);
end

function result = funA(omega, kpar, n, l)
    result = (l.Tperp-l.Tpara)/(omega*l.Tpara) + ...
        ( (funXi(omega, kpar, n, l)*l.Tperp)/(omega*l.Tpara) + ...
          n*l.omega_c./(kpar*omega*l.w_par) ) .* Z0(funXi(omega, kpar, n, l), kpar);
end

function result = funB0(omega, kpar, l)
    result = (omega*l.Tperp - kpar*l.V*l.Tpara)/(omega*kpar*l.Tpara) + ...
        (funXi(omega, kpar, 0, l)*l.Tperp)./(kpar*l.Tpara) .* ...
        Z0(funXi(omega, kpar, 0, l), kpar);
end

function result = funXi(omega, kpar, n, l)
    result = (omega - kpar*l.V - n*l.omega_c) ./ (kpar*l.w_par);
end

% 主扫描循环
for i = 1:num_B
    B0 = B_values(i);
    
    fprintf('处理 B = %.2f T (%d/%d)...\n', B0, i, num_B);
    
    for l = [ele, ion]
        l.omega_p = sqrt(n*e^2/(eps0*l.m));    % 等离子体频率
        l.omega_c = e*B0/l.m;                  % 回旋频率
        l.w_para = sqrt(2*kB*l.T/l.m);          % 平行热速
        l.w_perp = sqrt(2*kB*l.Tperp/l.m);     % 垂直热速
    end
    
    % 存储频率
    omega_c_i_results(i) = ion.omega_c;
    omega_c_e_results(i) = ele.omega_c;

    % 计算冷等离子体介电常数分量
    K_perp_cold(i) = 1 - (ele.omega_p^2/(omega^2 - ele.omega_c^2)) - (ion.omega_p^2/(omega^2 - ion.omega_c^2));
    K_phi_cold(i) = (ele.omega_c*ele.omega_p^2/(omega^3 - omega*ele.omega_c^2)) + (ion.omega_c*ion.omega_p^2/(omega^3 - omega*ion.omega_c^2));
    K_par_cold(i) = 1 - (ele.omega_p^2/omega^2) - (ion.omega_p^2/omega^2);
    
    % 先解冷等离子体色散方程，得到k_par_cold
    k0_cold = 10+0.1i;
    try
        options_cold = optimset('Display','off','TolFun',1e-10,'TolX',1e-10);
        kpar_cold_sol = fsolve(@(kpar) coldPDE(kpar, omega, ele, ion), k0_cold, options_cold);
        k_real_cold_results(i) = real(kpar_cold_sol);
        k_imag_cold_results(i) = imag(kpar_cold_sol);
    catch
        k_real_cold_results(i) = NaN;
        k_imag_cold_results(i) = NaN;
    end
    
    % 尝试求解热等离子体色散方程
    try
        step = 10;
        k0_hot = k0_cold;
        options = optimset('Display','iter','TolFun',1e-30,'TolX',1e-20);
        for istep = 1:step
            for l = [ele, ion]
                l.T_perp = Tperp_e * (istep/step);
                l.T_para = Tpara_e * (istep/step);
            end
            kpar_sol = fsolve(@(kpar) hotPDE(kpar, omega, ele, ion), k0_hot, options);
            k0_hot = kpar_sol;
        end
        
        % 存储k结果
        k_real_results(i) = real(kpar_sol);
        k_imag_results(i) = imag(kpar_sol);
        success_flag(i) = 1;
        
        % 计算热等离子体介电张量分量
        % 初始化
        K_perp = 1;
        K_phi = 0;
        K_par = 1;
        
        for l = [ele, ion]
            K_perp = K_perp + (l.omega_p^2/(2*omega)) * (funA(omega, kpar_sol, -1, l) + funA(omega, kpar_sol, 1, l));
            K_phi = K_phi + (l.omega_p^2/(2*omega)) * (funA(omega, kpar_sol, -1, l) - funA(omega, kpar_sol, 1, l));
            K_par = K_par + (2*l.omega_p^2/(kpar_sol*l.w_perp^2)) * (l.V/omega + funB0(omega, kpar_sol, l));
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

% 计算共振磁场（omega_c_i = omega）
B_resonance = omega * ion.m / e;
fprintf('\n共振磁场: B = %.3f T (omega_c_i = omega)\n', B_resonance);

% 输出结果
fprintf('\n=== 扫描结果 ===\n');
fprintf('B(T)\t\tomega_c_i(rad/s)\t\tomega/omega_c_i\t\tk_real(1/m)\t\tk_imag(1/m)\t\t状态\n');
fprintf('----\t\t-------------\t\t-------------\t\t----------\t\t----------\t\t----\n');

for i = 1:num_B
    if success_flag(i)
        fprintf('%.2f\t\t%.2e\t\t%.3f\t\t%.6e\t\t%.6e\t\t成功\n', ...
            B_values(i), omega_c_i_results(i), omega/omega_c_i_results(i), ...
            k_real_results(i), k_imag_results(i));
    else
        fprintf('%.2f\t\t%.2e\t\t%.3f\t\t%s\t\t%s\t\t失败\n', ...
            B_values(i), omega_c_i_results(i), omega/omega_c_i_results(i), ...
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

% 子图6: omega/omega_c_i比值
subplot(3, 3, 7);
omega_c_i_ratio = omega ./ omega_c_i_results;
plot(B_values, omega_c_i_ratio, 'b-*', 'LineWidth', 1.5, 'MarkerSize', 3);
xline(B_resonance, '--k', 'LineWidth', 2);
yline(1, '--r', 'LineWidth', 1); % 标记共振线
xlabel('磁场强度 B (T)');
ylabel('\omega/\omega_{ci}');
title('\omega/\omega_{ci} vs 磁场强度');
grid on;

% 子图7: 归一化的k实部
subplot(3, 3, 8);
k_real_norm = c * k_real_results ./ omega_c_i_results;
plot(B_values, k_real_norm, 'b-*', 'LineWidth', 1.5, 'MarkerSize', 3);
xline(B_resonance, '--k', 'LineWidth', 2);
xlabel('磁场强度 B (T)');
ylabel('c*k_{real}/\omega_{ci}');
title('归一化k实部 vs 磁场强度');
grid on;

% 子图8: 归一化的k虚部
subplot(3, 3, 9);
k_imag_norm = c * k_imag_results ./ omega_c_i_results;
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
results_table = table(B_values', omega_c_i_results', omega_c_e_results', ...
    omega./omega_c_i_results', k_real_results', k_imag_results', ...
    real(K_perp_hot)', imag(K_perp_hot)', real(K_phi_hot)', imag(K_phi_hot)', ...
    real(K_par_hot)', imag(K_par_hot)', K_perp_cold', K_phi_cold', K_par_cold', ...
    success_flag', ...
    'VariableNames', {'B_T', 'omega_c_i_rad_s', 'omega_c_e_rad_s', ...
    'omega_omega_c_i_ratio', 'k_real_1_m', 'k_imag_1_m', ...
    'K_perp_hot_real', 'K_perp_hot_imag', 'K_phi_hot_real', 'K_phi_hot_imag', ...
    'K_par_hot_real', 'K_par_hot_imag', 'K_perp_cold', 'K_phi_cold', 'K_par_cold', ...
    'success'});

writetable(results_table, 'dil_solver_scan_B_results.csv');
fprintf('\n结果已保存到 dil_solver_scan_B_results.csv\n');
