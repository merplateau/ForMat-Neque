function kpar_sol=solvek(M,ni,ne,B0,f,Te,Ti,Tperp_e,Tpara_e,Tperp_i,Tpara_i,V_e,V_i)

    % 基本常数
    e = 1.602e-19;           % 元电荷 (C)
    me = 9.11e-31;           % 电子质量 (kg)
    mi = M*1.67e-27;         % 离子质量 (kg)，使用参数M
    eps0 = 8.854e-12;        % 真空介电常数 (F/m)
    c = 3e8;                 % 光速 (m/s)
    kB = 1.38e-23;           % 玻尔兹曼常数 (J/K)

    % VASIMR ICRH 典型参数
    %ne = 1e18;               % 电子密度 (m^-3)
    %ni = ne;                 % 离子密度 (m^-3)，假设准中性
    %B0 = 0.5;                % 磁场 (T)
    %f = 13.56e6;                % 驱动频率 (Hz)
    %Te_eV = 3;              % 电子温度 (eV)
    %Ti_eV = 0.3;               % 离子温度 (eV)
    %Te = Te_eV * e / kB;     % 电子温度 (K)
    %Ti = Ti_eV * e / kB;     % 离子温度 (K)
    %Tperp_e = 10*Te;         % 电子垂直温度 (K)
    %Tpara_e = 1*Te;         % 电子平行温度 (K)
    %Tperp_i = 10*Ti;         % 离子垂直温度 (K)
    %Tpara_i = 1*Ti;         % 离子平行温度 (K)

    % 频率
    omega = 2*pi*f;
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
    %V_e = 0;                              % 电子漂移速度
    %V_i = 0;                              % 离子漂移速度

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

end

% 输出归一化结果
%result_real = c * real(kpar_sol) / omega_ci;
%result_imag = c * imag(kpar_sol) / omega_ci;

% =====================
% zetaf函数 - 等离子体色散函数
function w=zetaf(z)
% 已弃用，使用本地faddeeva实现
w = 1i*sqrt(pi)*faddeeva(z);
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

% 完整色散方程（包含电子和离子）- 归一化版本
function F = my_disp_eq(kpar, omega, omega_pe, omega_pi, c, A1e, A1i, xi_n_e, xi_n_i, Z0)
    Ae = A1e(kpar);
    Ai = A1i(kpar);
    xi_e = xi_n_e(kpar);
    xi_i = xi_n_i(kpar);
    Z_e = Z0(xi_e, kpar);
    Z_i = Z0(xi_i, kpar);
    fprintf('当前kpar = %.6e + %.6ei, xi_e = %.6e + %.6ei, xi_i = %.6e + %.6ei, Z_e = %.6e + %.6ei, Z_i = %.6e + %.6ei, Ae = %.6e + %.6ei, Ai = %.6e + %.6ei\n', ...
        real(kpar), imag(kpar), real(xi_e), imag(xi_e), real(xi_i), imag(xi_i), real(Z_e), imag(Z_e), real(Z_i), imag(Z_i), real(Ae), imag(Ae), real(Ai), imag(Ai));
    F = (kpar.^2 * c^2 - omega^2 - omega_pe^2 * Ae - omega_pi^2 * Ai) / c^2;
end