% 基本常数
e = 1.602e-19;           % 元电荷 (C)
me = 9.11e-31;           % 电子质量 (kg)
mi = 1.67e-27;           % 质子质量 (kg)
eps0 = 8.854e-12;        % 真空介电常数 (F/m)
c = 3e8;                 % 光速 (m/s)
kB = 1.38e-23;           % 玻尔兹曼常数 (J/K)

% VASIMR ICRH 典型参数
ne = 1e18;               % 电子密度 (m^-3)
ni = ne;                 % 离子密度 (m^-3)，假设准中性
B0 = 0.5;                % 磁场 (T)
f = 30e6;                % 驱动频率 (Hz)
Te_eV = 3;               % 电子温度 (eV)
Ti_eV = 0.3;             % 离子温度 (eV)
Te = Te_eV;              % 电子温度 (eV)
Ti = Ti_eV;              % 离子温度 (eV)
Tperp_e = Te;            % 电子垂直温度 (eV)
Tpara_e = Te;            % 电子平行温度 (eV)
Tperp_i = Ti;            % 离子垂直温度 (eV)
Tpara_i = Ti;            % 离子平行温度 (eV)

% 频率
omega = 2*pi*f;
omega_pe = sqrt(ne*e^2/(eps0*me));    % 电子等离子体频率
omega_pi = sqrt(ni*e^2/(eps0*mi));    % 离子等离子体频率
omega_ce = e*B0/me;                   % 电子回旋频率
omega_ci = e*B0/mi;                   % 离子回旋频率

% 热速 (eV -> K -> m/s)
Te_K = Te_eV * 11604;                 % 1eV = 11604K
Ti_K = Ti_eV * 11604;
w_para_e = sqrt(2*kB*Te_K/me);        % 电子平行热速
w_para_i = sqrt(2*kB*Ti_K/mi);        % 离子平行热速

% 谐波数
n_e = 1;                              % 电子谐波数
n_i = 1;                              % 离子谐波数
V_e = 0;                              % 电子漂移速度
V_i = 0;                              % 离子漂移速度

% 色散函数Z(ξ)的实现
Z = @(xi) plasma_dispersion_func(xi);

% Z0(ξ)的实现
Z0 = @(xi, kpar) (kpar>0).*Z(xi) - (kpar<0).*Z(-xi);

% xi_n^l(kpar)的实现
xi_n_e = @(kpar) (omega - kpar*V_e - n_e*omega_ce) ./ (kpar*w_para_e);
xi_n_i = @(kpar) (omega - kpar*V_i - n_i*omega_ci) ./ (kpar*w_para_i);

% A_{+1}^l(kpar)的实现
A1e = @(kpar) (Tperp_e-Tpara_e)/(omega*Tpara_e) + ...
    ( (xi_n_e(kpar)*Tperp_e)/(omega*Tpara_e) + n_e*omega_ce./(omega*kpar*w_para_e*Tpara_e*e) ) .* Z0(xi_n_e(kpar), kpar);

A1i = @(kpar) (Tperp_i-Tpara_i)/(omega*Tpara_i) + ...
    ( (xi_n_i(kpar)*Tperp_i)/(omega*Tpara_i) + n_i*omega_ci./(omega*kpar*w_para_i*Tpara_i*e) ) .* Z0(xi_n_i(kpar), kpar);

% 完整色散方程（包含电子和离子）
disp_eq = @(kpar) kpar.^2 * c^2 - omega^2 - omega_pe^2 * A1e(kpar) - omega_pi^2 * A1i(kpar);

% 初始猜测 - 尝试多个值
k0_1 = omega/c;                    % 光速近似
k0_2 = omega/sqrt(c^2 - omega_pe^2/omega^2);  % 冷等离子体近似
k0_3 = omega_ce/c;                 % 回旋频率近似

% 求解选项
options = optimset('Display','iter','TolFun',1e-12,'TolX',1e-12,'MaxIter',100);

% 尝试多个初始值
kpar_sol = [];
result = [];

try
    kpar_sol = fsolve(disp_eq, k0_1, options);
    result = c * kpar_sol / omega_ce;
    fprintf('初始猜测1成功: c*k_parallel/omega_ce = %.6f\n', result);
catch
    try
        kpar_sol = fsolve(disp_eq, k0_2, options);
        result = c * kpar_sol / omega_ce;
        fprintf('初始猜测2成功: c*k_parallel/omega_ce = %.6f\n', result);
    catch
        try
            kpar_sol = fsolve(disp_eq, k0_3, options);
            result = c * kpar_sol / omega_ce;
            fprintf('初始猜测3成功: c*k_parallel/omega_ce = %.6f\n', result);
        catch
            fprintf('所有初始猜测都失败\n');
            return;
        end
    end
end

% 输出各组分贡献
fprintf('电子贡献: %.6e\n', omega_pe^2 * A1e(kpar_sol));
fprintf('离子贡献: %.6e\n', omega_pi^2 * A1i(kpar_sol));

% 输出参数检查
fprintf('=== 参数检查 ===\n');
fprintf('omega = %.2e rad/s\n', omega);
fprintf('omega_pe = %.2e rad/s\n', omega_pe);
fprintf('omega_pi = %.2e rad/s\n', omega_pi);
fprintf('omega_ce = %.2e rad/s\n', omega_ce);
fprintf('omega_ci = %.2e rad/s\n', omega_ci);
fprintf('w_para_e = %.2e m/s\n', w_para_e);
fprintf('w_para_i = %.2e m/s\n', w_para_i);
fprintf('omega/omega_ce = %.2f\n', omega/omega_ce);
fprintf('omega/omega_ci = %.2f\n', omega/omega_ci);

% 测试冷等离子体极限（简化版本）
fprintf('\n=== 测试冷等离子体极限 ===\n');
k_cold = sqrt((omega^2 - omega_pe^2 - omega_pi^2)/c^2);
result_cold = c * k_cold / omega_ce;
fprintf('冷等离子体近似: c*k_parallel/omega_ce = %.6f\n', result_cold);

% 检查A_{+1}^l函数在初始猜测处的值
fprintf('\n=== A_{+1}^l函数检查 ===\n');
k_test = omega/c;
A1e_test = A1e(k_test);
A1i_test = A1i(kpar_sol);
fprintf('A1e(k_test) = %.6e\n', A1e_test);
fprintf('A1i(k_test) = %.6e\n', A1i_test);
fprintf('xi_n_e(k_test) = %.6e\n', xi_n_e(k_test));
fprintf('xi_n_i(k_test) = %.6e\n', xi_n_i(k_test));

% 运行基本测试
test_basic_calculations();

% =====================
% 自定义等离子体色散函数（数值积分近似）
function z = plasma_dispersion_func(xi)
    % 数值积分近似Z函数
    z = zeros(size(xi));
    for k = 1:numel(xi)
        if abs(xi(k)) > 10
            % 大参数渐近展开
            z(k) = -1/xi(k) * (1 + 1/(2*xi(k)^2) + 3/(4*xi(k)^4));
        else
            % 数值积分
            fun = @(t) exp(-t.^2) ./ (t - xi(k));
            try
                z(k) = 1i*sqrt(pi) * integral(fun, -Inf, Inf, 'ArrayValued', true, 'RelTol',1e-8,'AbsTol',1e-10);
            catch
                % 如果积分失败，用渐近展开
                z(k) = -1/xi(k) * (1 + 1/(2*xi(k)^2) + 3/(4*xi(k)^4));
            end
        end
    end
end

% =====================
% 测试函数
function test_basic_calculations()
    fprintf('\n=== 基本计算测试 ===\n');
    
    % 测试色散函数
    xi_test = 1.0;
    z_test = plasma_dispersion_func(xi_test);
    fprintf('Z(1.0) = %.6f + %.6fi\n', real(z_test), imag(z_test));
    
    % 测试Z0函数
    z0_test = (1>0)*z_test - (1<0)*plasma_dispersion_func(-xi_test);
    fprintf('Z0(1.0, 1) = %.6f + %.6fi\n', real(z0_test), imag(z0_test));
    
    % 测试渐近展开
    xi_large = 15.0;
    z_large = plasma_dispersion_func(xi_large);
    fprintf('Z(15.0) = %.6f + %.6fi\n', real(z_large), imag(z_large));
end 