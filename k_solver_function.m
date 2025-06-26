function [result_real, result_imag] = k_solver_function(varargin)
% K_SOLVER_FUNCTION - 函数版本的k_solver，可以接受参数并返回结果
% 用法：
%   [result_real, result_imag] = k_solver_function()
%   [result_real, result_imag] = k_solver_function('param', value, ...)
%
% 可选参数：
%   'ne' - 电子密度 (m^-3)，默认 1e18
%   'B0' - 磁场 (T)，默认 0.5
%   'f' - 驱动频率 (Hz)，默认 13.56e6
%   'Te_eV' - 电子温度 (eV)，默认 3
%   'Ti_eV' - 离子温度 (eV)，默认 0.3

% 解析输入参数
p = inputParser;
addParameter(p, 'ne', 1e18);
addParameter(p, 'B0', 0.5);
addParameter(p, 'f', 13.56e6);
addParameter(p, 'Te_eV', 3);
addParameter(p, 'Ti_eV', 0.3);
parse(p, varargin{:});

% 获取参数
ne = p.Results.ne;
B0 = p.Results.B0;
f = p.Results.f;
Te_eV = p.Results.Te_eV;
Ti_eV = p.Results.Ti_eV;

% 基本常数
e = 1.602e-19;           % 元电荷 (C)
me = 9.11e-31;           % 电子质量 (kg)
mi = 40*1.67e-27;        % 质子质量 (kg)
eps0 = 8.854e-12;        % 真空介电常数 (F/m)
c = 3e8;                 % 光速 (m/s)
kB = 1.38e-23;           % 玻尔兹曼常数 (J/K)

% 设置参数
ni = ne;                 % 离子密度 (m^-3)，假设准中性
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
omega = 0.999*omega_ci;

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

% xi_n^l(kpar)的实现
xi_n_e = @(kpar) (omega - kpar*V_e - n_e*omega_ce) ./ (kpar*w_para_e);
xi_n_i = @(kpar) (omega - kpar*V_i - n_i*omega_ci) ./ (kpar*w_para_i);

% A_{+1}^l(kpar)的实现
A1e = @(kpar) (Tperp_e-Tpara_e)/(omega*Tpara_e) + ...
    ( (xi_n_e(kpar)*Tperp_e)/(omega*Tpara_e) + n_e*omega_ce./(omega*kpar*w_para_e) ) .* Z0(xi_n_e(kpar), kpar);

A1i = @(kpar) (Tperp_i-Tpara_i)/(omega*Tpara_i) + ...
    ( (xi_n_i(kpar)*Tperp_i)/(omega*Tpara_i) + n_i*omega_ci./(omega*kpar*w_para_i) ) .* Z0(xi_n_i(kpar), kpar);

% 初始猜测
k0 = 1 * omega / c;

% 求解
options = optimset('Display','off','TolFun',1e-10,'TolX',1e-10);
kpar_sol = fsolve(@(kpar) my_disp_eq(kpar, omega, omega_pe, omega_pi, c, A1e, A1i, xi_n_e, xi_n_i, Z0), k0, options);

% 输出归一化结果
result_real = c * real(kpar_sol) / omega_ci;
result_imag = c * imag(kpar_sol) / omega_ci;

% 输出结果（可选）
if nargout == 0
    fprintf('Re[c*k_parallel/omega_ci] = %.6f\n', result_real);
    fprintf('Im[c*k_parallel/omega_ci] = %.6f\n', result_imag);
end

end

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