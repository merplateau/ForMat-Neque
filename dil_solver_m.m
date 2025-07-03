 % 介电张量求解程序，基于k_solver.m
% 1. 先求解色散方程，获得kpar
% 2. 用kpar代入介电张量公式，输出K_perp, K_phi, K_par

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
f = 190e3;             % 驱动频率 (Hz)
Te_eV = 0.3;               % 电子温度 (eV)
Ti_eV = 0.3;             % 离子温度 (eV)
Te = Te_eV * e / kB;     % 电子温度 (K)
Ti = Ti_eV * e / kB;     % 离子温度 (K)
Tperp_e = 1*Te;         % 电子垂直温度 (K)
Tpara_e = 1*Te;          % 电子平行温度 (K)
Tperp_i = 1*Ti;         % 离子垂直温度 (K)
Tpara_i = 1*Ti;          % 离子平行温度 (K)

% 频率
omega = 2*pi*f;
omega_pe = sqrt(ne*e^2/(eps0*me));    % 电子等离子体频率
omega_pi = sqrt(ni*e^2/(eps0*mi));    % 离子等离子体频率
omega_ce = e*B0/me;                   % 电子回旋频率
omega_ci = e*B0/mi;                   % 离子回旋频率

% 热速 (K -> m/s)
w_para_e = sqrt(2*kB*Te/me);          % 电子平行热速
w_para_i = sqrt(2*kB*Ti/mi);          % 离子平行热速
w_perp_e = sqrt(2*kB*Tperp_e/me);     % 电子垂直热速
w_perp_i = sqrt(2*kB*Tperp_i/mi);     % 离子垂直热速

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
xi_n = @(omega, kpar, V_l, n_l, omega_cl, w_par_l) (omega - kpar*V_l - n_l*omega_cl) ./ (kpar*w_par_l);

% A_l^n(kpar)的实现
A_l_n = @(Tperp_l, Tpara_l, omega, kpar, V_l, n_l, omega_cl, w_par_l) ...
    (Tperp_l-Tpara_l)/(omega*Tpara_l) + ...
    ( (xi_n(omega, kpar, V_l, n_l, omega_cl, w_par_l)*Tperp_l)/(omega*Tpara_l) + ...
      n_l*omega_cl./(omega*kpar*w_par_l) ) .* Z0(xi_n(omega, kpar, V_l, n_l, omega_cl, w_par_l), kpar);

% B_l^0(kpar)的实现
B_l_0 = @(Tperp_l, Tpara_l, omega, kpar, V_l, w_par_l) ...
    (omega*Tperp_l - kpar*V_l*Tpara_l)/(omega*kpar*Tpara_l) + ...
    (xi_n(omega, kpar, V_l, 0, 0, w_par_l)*Tperp_l)./(kpar*Tpara_l) .* ...
    Z0(xi_n(omega, kpar, V_l, 0, 0, w_par_l), kpar);

% 色散方程
A1e = @(kpar) A_l_n(Tperp_e, Tpara_e, omega, kpar, V_e, n_e, omega_ce, w_para_e);
A1i = @(kpar) A_l_n(Tperp_i, Tpara_i, omega, kpar, V_i, n_i, omega_ci, w_para_i);

function F = my_disp_eq(kpar, omega, omega_pe, omega_pi, c, A1e, A1i)
    F = (kpar.^2 * c^2 - omega^2 - omega_pe^2 * A1e(kpar) - omega_pi^2 * A1i(kpar)) / c^2;
end

% 初始猜测
k0 = 1000 * omega / c;

% 求解kpar
options = optimset('Display','iter','TolFun',1e-10,'TolX',1e-10);
kpar_sol = fsolve(@(kpar) my_disp_eq(kpar, omega, omega_pe, omega_pi, c, A1e, A1i), k0, options);

fprintf('求解得到 k_parallel = %.6e + %.6ei 1/m\n', real(kpar_sol), imag(kpar_sol));

% ========== 计算介电张量分量 ==========
% 电子参数
params_e = struct('omega_p', omega_pe, 'omega_c', omega_ce, 'Tperp', Tperp_e, 'Tpara', Tpara_e, ...
    'V', V_e, 'n', n_e, 'w_par', w_para_e, 'w_perp', w_perp_e);
% 离子参数
params_i = struct('omega_p', omega_pi, 'omega_c', omega_ci, 'Tperp', Tperp_i, 'Tpara', Tpara_i, ...
    'V', V_i, 'n', n_i, 'w_par', w_para_i, 'w_perp', w_perp_i);

% K_perp, K_phi, K_par
K_perp = 1;
K_phi = 0;
K_par = 1;
for l = [params_e, params_i]
    % A_l^+1, A_l^-1
    A_p1 = A_l_n(l.Tperp, l.Tpara, omega, kpar_sol, l.V, +1, l.omega_c, l.w_par);
    A_m1 = A_l_n(l.Tperp, l.Tpara, omega, kpar_sol, l.V, -1, l.omega_c, l.w_par);
    % B_l^0
    B_0 = B_l_0(l.Tperp, l.Tpara, omega, kpar_sol, l.V, l.w_par);
    % K_perp
    K_perp = K_perp + (l.omega_p^2/(2*omega)) * (A_m1 + A_p1);
    % K_phi
    K_phi = K_phi + (l.omega_p^2/(2*omega)) * (A_m1 - A_p1);
    % K_par
    K_par = K_par + (2*l.omega_p^2/(kpar_sol*l.w_perp^2)) * (l.V/omega + B_0);
end

fprintf('\n介电张量分量：\n');
fprintf('K_perp = %.6e + %.6ei\n', real(K_perp), imag(K_perp));
fprintf('K_phi  = %.6e + %.6ei\n', real(K_phi), imag(K_phi));
fprintf('K_par  = %.6e + %.6ei\n', real(K_par), imag(K_par));

% ========== 冷等离子体介电常数分量 ==========
K_perp_cold = 1 - (omega_pe^2/(omega^2 - omega_ce^2)) - (omega_pi^2/(omega^2 - omega_ci^2));
K_phi_cold = (omega_ce*omega_pe^2/(omega^3 - omega*omega_ce^2)) + (omega_ci*omega_pi^2/(omega^3 - omega*omega_ci^2));
K_par_cold = 1 - (omega_pe^2/omega^2) - (omega_pi^2/omega^2);

fprintf('\n冷等离子体介电常数分量：\n');
fprintf('K_perp_cold = %.6e\n', K_perp_cold);
fprintf('K_phi_cold  = %.6e\n', K_phi_cold);
fprintf('K_par_cold  = %.6e\n', K_par_cold);

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
