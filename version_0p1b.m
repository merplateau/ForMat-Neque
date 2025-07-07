clear; clc;

% 设置Python路径
python_path = 'D:/ProgramData/anaconda3/envs/plasma/python.exe';
if ~exist(python_path, 'file')
    error('Python路径不存在');
end
pe = pyenv('Version', python_path);
if isempty(pe.Executable)
    error('无法配置Python环境');
end
scipy_special = py.importlib.import_module('scipy.special');

% 基本常数
e = 1.602e-19;
me = 9.11e-31;
mi = 40*1.67e-27;
eps0 = 8.854e-12;
c = 3e8;
kB = 1.38e-23;

% VASIMR ICRH 典型参数
n = 1e18;
ne = n;
ni = n;
f = 190e3;
Te_eV = 10;
Ti_eV = 1;
Te = Te_eV * e / kB;
Ti = Ti_eV * e / kB;
Tperp_e = Te;
Tpara_e = Te;
Tperp_i = Ti;
Tpara_i = Ti;

omega = 2*pi*f;
B0 = 0.48; % 固定一个B

% 电子/离子参数
V_e = 0;
V_i = 0;
ele = struct('omega_p', sqrt(n*e^2/(eps0*me)), 'omega_c', e*B0/me, 'T_perp', Tperp_e, 'T_para', Tpara_e, 'V', V_e, 'w_para', sqrt(2*kB*Tpara_e/me), 'w_perp', sqrt(2*kB*Tperp_e/me), 'm', me);
ion = struct('omega_p', sqrt(n*e^2/(eps0*mi)), 'omega_c', e*B0/mi, 'T_perp', Tperp_i, 'T_para', Tpara_i, 'V', V_i, 'w_para', sqrt(2*kB*Tpara_i/mi), 'w_perp', sqrt(2*kB*Tperp_i/mi), 'm', mi);

% kpar扫描范围
k_real = linspace(2.6324698, 2.63246985, 500);
k_imag = linspace(2.5612955, 2.5612956, 500);
[K_real, K_imag] = meshgrid(k_real, k_imag);
Kpar = K_real + 1i*K_imag;

hotPDE_val = zeros(size(Kpar));
for idx = 1:numel(Kpar)
    kpar = Kpar(idx);
    hotPDE_val(idx) = abs(hotPDE(kpar, omega, ele, ion, c, scipy_special));
end
hotPDE_val = reshape(hotPDE_val, size(Kpar));

figure;
surf(K_real, K_imag, log10(hotPDE_val+1e-20));
shading interp;
colorbar;
xlabel('k_{par} 实部 (1/m)');
ylabel('k_{par} 虚部 (1/m)');
zlabel('log_{10}|hotPDE|');
title('hotPDE在k_{par}复平面上的分布');

% hotPDE函数及相关子函数
function result = hotPDE(kpar, omega, ele, ion, c, scipy_special)
    result = kpar^2 * c^2 - omega^2 - ele.omega_p^2 * omega * funA(omega, kpar, 1, ele, c, scipy_special) - ion.omega_p^2 * omega * funA(omega, kpar, 1, ion, c, scipy_special);
end

function result = funA(omega, kpar, n, l, c, scipy_special)
    result = (l.T_perp-l.T_para)/(omega*l.T_para) + ...
        ( (funXi(omega, kpar, n, l)*l.T_perp)/(omega*l.T_para) + ...
          (n * l.omega_c)/(kpar*omega*l.w_para) ) .* Z0(funXi(omega, kpar, n, l), kpar, scipy_special);
end

function result = funXi(omega, kpar, n, l)
    result = (omega - kpar*l.V - n*l.omega_c) ./ (kpar*l.w_para);
end

function result = Z0(xi, kpar, scipy_special)
    z1 = xi * (real(kpar)>0) - xi *(real(kpar)<0);
    z2 = real(z1) + abs(imag(z1))*1i;
    result = 1i*sqrt(pi)*scipy_special.wofz(z2);
end 