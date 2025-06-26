% 分析求解过程中遇到奇异时各项的情况
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
B0 = 0.5;                % 磁场 (T)
Te_eV = 3;               % 电子温度 (eV)
Ti_eV = 0.3;             % 离子温度 (eV)
Te = Te_eV * e / kB;     % 电子温度 (K)
Ti = Ti_eV * e / kB;     % 离子温度 (K)
Tperp_e = 10*Te;         % 电子垂直温度 (K)
Tpara_e = 1*Te;          % 电子平行温度 (K)
Tperp_i = 10*Ti;         % 离子垂直温度 (K)
Tpara_i = 1*Ti;          % 离子平行温度 (K)

% 频率
omega_pe = sqrt(ne*e^2/(eps0*me));    % 电子等离子体频率
omega_pi = sqrt(ni*e^2/(eps0*mi));    % 离子等离子体频率
omega_ce = e*B0/me;                   % 电子回旋频率
omega_ci = e*B0/mi;                   % 离子回旋频率

% 热速 (K -> m/s)
w_para_e = sqrt(2*kB*Te/me);          % 电子平行热速
w_para_i = sqrt(2*kB*Ti/mi);          % 离子平行热速

% 谐波数
n_e = 1;                              % 电子谐波数
n_i = 1;                              % 离子谐波数
V_e = 0;                              % 电子漂移速度
V_i = 0;                              % 离子漂移速度

% 色散函数Z(ξ)的实现 - 使用faddeeva.m
Z = @(xi)  1i*sqrt(pi)*faddeeva(xi,16);

% Z0(ξ)的实现
Z0 = @(xi, kpar) Z(xi) .* (real(kpar) > 0) - Z(-xi) .* (real(kpar) < 0);

fprintf('=== 基本参数 ===\n');
fprintf('omega_ci = %.2e rad/s\n', omega_ci);
fprintf('omega_pe = %.2e rad/s\n', omega_pe);
fprintf('omega_pi = %.2e rad/s\n', omega_pi);
fprintf('w_para_i = %.2e m/s\n', w_para_i);

% 测试不同的omega值，包括接近omega_ci的情况
omega_ratios = [0.99, 0.999, 1.0, 1.001, 1.01];

 function F = my_disp_eq(kpar, omega, omega_pe, omega_pi, c, A1e, A1i, xi_n_e, xi_n_i, Z0)
        Ae = A1e(kpar);
        Ai = A1i(kpar);
        xi_e = xi_n_e(kpar);
        xi_i = xi_n_i(kpar);
        Z_e = Z0(xi_e, kpar);
        Z_i = Z0(xi_i, kpar);
        
        % 详细分析各项
        term1 = kpar.^2 * c^2;
        term2 = omega^2;
        term3 = omega_pe^2 * Ae;
        term4 = omega_pi^2 * Ai;
        
        fprintf('kpar = %.6e + %.6ei\n', real(kpar), imag(kpar));
        fprintf('  xi_e = %.6e + %.6ei\n', real(xi_e), imag(xi_e));
        fprintf('  xi_i = %.6e + %.6ei\n', real(xi_i), imag(xi_i));
        fprintf('  Z_e = %.6e + %.6ei\n', real(Z_e), imag(Z_e));
        fprintf('  Z_i = %.6e + %.6ei\n', real(Z_i), imag(Z_i));
        fprintf('  Ae = %.6e + %.6ei\n', real(Ae), imag(Ae));
        fprintf('  Ai = %.6e + %.6ei\n', real(Ai), imag(Ai));
        fprintf('  各项贡献:\n');
        fprintf('    k^2*c^2 = %.6e + %.6ei\n', real(term1), imag(term1));
        fprintf('    omega^2 = %.6e + %.6ei\n', real(term2), imag(term2));
        fprintf('    omega_pe^2*Ae = %.6e + %.6ei\n', real(term3), imag(term3));
        fprintf('    omega_pi^2*Ai = %.6e + %.6ei\n', real(term4), imag(term4));
        
        F = (term1 - term2 - term3 - term4) / c^2;
        fprintf('  F = %.6e + %.6ei\n', real(F), imag(F));
        fprintf('  |F| = %.6e\n', abs(F));
        fprintf('  ---\n');
    end

for i = 1:length(omega_ratios)
    omega = omega_ratios(i) * omega_ci;
    
    fprintf('\n========================================\n');
    fprintf('测试 omega/omega_ci = %.3f\n', omega_ratios(i));
    fprintf('========================================\n');
    
    % xi_n^l(kpar)的实现
    xi_n_e = @(kpar) (omega - kpar*V_e - n_e*omega_ce) ./ (kpar*w_para_e);
    xi_n_i = @(kpar) (omega - kpar*V_i - n_i*omega_ci) ./ (kpar*w_para_i);
    
    % A_{+1}^l(kpar)的实现
    A1e = @(kpar) (Tperp_e-Tpara_e)/(omega*Tpara_e) + ...
        ( (xi_n_e(kpar)*Tperp_e)/(omega*Tpara_e) + n_e*omega_ce./(omega*kpar*w_para_e) ) .* Z0(xi_n_e(kpar), kpar);
    
    A1i = @(kpar) (Tperp_i-Tpara_i)/(omega*Tpara_i) + ...
        ( (xi_n_i(kpar)*Tperp_i)/(omega*Tpara_i) + n_i*omega_ci./(omega*kpar*w_para_i) ) .* Z0(xi_n_i(kpar), kpar);
    
    % 完整色散方程（包含电子和离子）- 归一化版本
    
    % 初始猜测
    k0 = 1 * omega / c;
    
    fprintf('初始猜测: k0 = %.6e + %.6ei\n', real(k0), imag(k0));
    
    % 尝试求解
    try
        options = optimset('Display','iter','TolFun',1e-10,'TolX',1e-10);
        kpar_sol = fsolve(@(kpar) my_disp_eq(kpar, omega, omega_pe, omega_pi, c, A1e, A1i, xi_n_e, xi_n_i, Z0), k0, options);
        
        % 输出归一化结果
        result_real = c * real(kpar_sol) / omega_ci;
        result_imag = c * imag(kpar_sol) / omega_ci;
        fprintf('求解成功!\n');
        fprintf('Re[c*k_parallel/omega_ci] = %.6f\n', result_real);
        fprintf('Im[c*k_parallel/omega_ci] = %.6f\n', result_imag);
        
    catch ME
        fprintf('求解失败: %s\n', ME.message);
        
        % 分析失败原因
        fprintf('\n分析失败原因:\n');
        
        % 检查xi_i是否接近0
        xi_i_test = xi_n_i(k0);
        fprintf('在初始猜测处 xi_i = %.6e + %.6ei\n', real(xi_i_test), imag(xi_i_test));
        
        if abs(xi_i_test) < 0.1
            fprintf('*** xi_i接近0，接近离子回旋共振 ***\n');
        end
        
        % 检查Z函数在xi_i接近0时的行为
        if abs(xi_i_test) < 1
            fprintf('Z函数在xi_i接近0时的行为:\n');
            xi_test_values = linspace(-0.1, 0.1, 21);
            for j = 1:length(xi_test_values)
                xi_test = xi_test_values(j);
                Z_test = Z(xi_test);
                fprintf('  xi = %.3f: Z = %.6e + %.6ei, |Z| = %.6e\n', ...
                    xi_test, real(Z_test), imag(Z_test), abs(Z_test));
            end
        end
        
        % 检查色散方程各项的相对大小
        Ae_test = A1e(k0);
        Ai_test = A1i(kpar);
        term1_test = k0^2 * c^2;
        term2_test = omega^2;
        term3_test = omega_pe^2 * Ae_test;
        term4_test = omega_pi^2 * Ai_test;
        
        fprintf('\n色散方程各项相对大小:\n');
        fprintf('|k^2*c^2| = %.2e\n', abs(term1_test));
        fprintf('|omega^2| = %.2e\n', abs(term2_test));
        fprintf('|omega_pe^2*Ae| = %.2e\n', abs(term3_test));
        fprintf('|omega_pi^2*Ai| = %.2e\n', abs(term4_test));
        
        % 检查是否有某项异常大
        max_term = max([abs(term1_test), abs(term2_test), abs(term3_test), abs(term4_test)]);
        if abs(term4_test) > max_term * 0.1
            fprintf('*** 警告: omega_pi^2*Ai项异常大，可能导致数值不稳定 ***\n');
        end
    end
end

% 分析物理机制
fprintf('\n========================================\n');
fprintf('物理机制分析\n');
fprintf('========================================\n');
fprintf('当omega接近omega_ci时:\n');
fprintf('1. xi_i = (omega - omega_ci)/(k_parallel*w_para_i) 接近0\n');
fprintf('2. 等离子体色散函数Z(xi_i)在xi_i≈0时表现出奇异行为\n');
fprintf('3. 这导致Ai项变得很大且数值不稳定\n');
fprintf('4. 色散方程变得病态，数值求解器无法收敛\n');
fprintf('5. 这是离子回旋共振的典型特征\n');


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