% k_solver_scan_omega.m
% 等离子体色散关系扫描程序
% 对 omega/omega_ce 从 0 到 1 进行扫描，步长 0.01

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
B0 = 5e-2;                % 磁场 (T)
f = 10e9;                % 驱动频率 (Hz)
Te_eV = 3;              % 电子温度 (eV)
Ti_eV = 0.3;               % 离子温度 (eV)
Te = Te_eV * e / kB;     % 电子温度 (K)
Ti = Ti_eV * e / kB;     % 离子温度 (K)
Tperp_e = Te;         % 电子垂直温度 (K)
Tpara_e = Te;         % 电子平行温度 (K)
Tperp_i = Ti;         % 离子垂直温度 (K)
Tpara_i = Ti;         % 离子平行温度 (K)

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

% 色散函数Z(ξ)的实现 - 使用zetaf.m
Z = @(xi)  zetaf(xi);

% Z0(ξ)的实现
Z0 = @(xi, kpar) Z(xi) .* (kpar > 0) - Z(-xi) .* (kpar < 0);

% 扫描参数
omega_ratios = 0.1:0.1:1;  % omega/omega_ce 从0到1，步长0.01
results_real = zeros(size(omega_ratios));  % 实部
results_imag = zeros(size(omega_ratios));  % 虚部
electron_contributions = zeros(size(omega_ratios));
ion_contributions = zeros(size(omega_ratios));

fprintf('开始扫描 omega/omega_ce 从 0 到 1...\n');
fprintf('总扫描点数: %d\n', length(omega_ratios));

for i = 1:length(omega_ratios)
    omega = omega_ratios(i) * omega_ce;
    
    % xi_n^l(kpar)的实现
    xi_n_e = @(kpar) (omega - kpar*V_e - n_e*omega_ce) ./ (kpar*w_para_e);
    xi_n_i = @(kpar) (omega - kpar*V_i - n_i*omega_ci) ./ (kpar*w_para_i);
    
    % A_{+1}^l(kpar)的实现
    A1e = @(kpar) (Tperp_e-Tpara_e)/(omega*Tpara_e) + ...
        ( (xi_n_e(kpar)*Tperp_e)/(omega*Tpara_e) + n_e*omega_ce./(omega*kpar*w_para_e*Tpara_e*e) ) .* Z0(xi_n_e(kpar), kpar);
    
    A1i = @(kpar) (Tperp_i-Tpara_i)/(omega*Tpara_i) + ...
        ( (xi_n_i(kpar)*Tperp_i)/(omega*Tpara_i) + n_i*omega_ci./(omega*kpar*w_para_i*Tpara_i*e) ) .* Z0(xi_n_i(kpar), kpar);
    
    % 完整色散方程（包含电子和离子）
    disp_eq = @(kpar)  (kpar.^2 * c^2 - omega^2 - omega_pe^2 * A1e(kpar) - omega_pi^2 * A1i(kpar)) / omega^2;
    
    % 初始猜测
    k0 = 5 * omega / c;
    
    % 求解
    try
        fprintf('iteration+1');
        options = optimset('Display','off','TolFun',1e-10,'TolX',1e-10);
        kpar_sol = fsolve(disp_eq, k0, options);
        
        % 存储结果
        results_real(i) = c * real(kpar_sol) / omega_ce;
        results_imag(i) = c * imag(kpar_sol) / omega_ce;
        electron_contributions(i) = omega_pe^2 * A1e(kpar_sol);
        ion_contributions(i) = omega_pi^2 * A1i(kpar_sol);
        
        % 检查解的实部和虚部
        if mod(i, 1) == 0  % 每5个点检查一次
            fprintf('omega/omega_ce = %.2f: k = %.6e + %.6ei\n', ...
                omega_ratios(i), real(kpar_sol), imag(kpar_sol));
            
            % 检查色散方程的值
            disp_val = disp_eq(kpar_sol);
            fprintf('  色散方程值 = %.6e + %.6ei\n', real(disp_val), imag(disp_val));
        end
        
        if mod(i, 20) == 0
            fprintf('进度: %.1f%% (omega/omega_ce = %.2f, c*k/omega_ce = %.4f)\n', ...
                100*i/length(omega_ratios), omega_ratios(i), results_real(i));
        end
    catch
        fprintf('omega/omega_ce = %.2f 求解失败\n', omega_ratios(i));
        results_real(i) = NaN;
        results_imag(i) = NaN;
        electron_contributions(i) = NaN;
        ion_contributions(i) = NaN;
    end
end

% 绘制结果
figure('Position', [100, 100, 1200, 800]);

% 子图1: c*k_parallel/omega_ce vs omega/omega_ce
subplot(2, 2, 1);
plot(omega_ratios, results_real, 'b-', 'LineWidth', 2, 'DisplayName', '实部');
hold on;
plot(omega_ratios, results_imag, 'r-', 'LineWidth', 2, 'DisplayName', '虚部');
xlabel('\omega/\omega_{ce}');
ylabel('c k_{\parallel}/\omega_{ce}');
title('等离子体色散关系');
legend('Location', 'best');
grid on;

% 子图2: 电子贡献
subplot(2, 2, 2);
plot(omega_ratios, electron_contributions, 'r-', 'LineWidth', 2);
xlabel('\omega/\omega_{ce}');
ylabel('\omega_{pe}^2 A_{1e}');
title('电子贡献');
grid on;

% 子图3: 离子贡献
subplot(2, 2, 3);
plot(omega_ratios, ion_contributions, 'g-', 'LineWidth', 2);
xlabel('\omega/\omega_{ce}');
ylabel('\omega_{pi}^2 A_{1i}');
title('离子贡献');
grid on;

% 子图4: 归一化贡献比较
subplot(2, 2, 4);
plot(omega_ratios, electron_contributions/omega_ce^2, 'r-', 'LineWidth', 2, 'DisplayName', '电子');
hold on;
plot(omega_ratios, ion_contributions/omega_ce^2, 'g-', 'LineWidth', 2, 'DisplayName', '离子');
xlabel('\omega/\omega_{ce}');
ylabel('归一化贡献');
title('归一化贡献比较');
legend('Location', 'best');
grid on;

% 输出一些关键结果
fprintf('\n=== 扫描结果摘要 ===\n');
fprintf('omega/omega_ce = 0.5 时: c*k_parallel/omega_ce = %.6f\n', results_real(omega_ratios == 0.5));
fprintf('omega/omega_ce = 0.8 时: c*k_parallel/omega_ce = %.6f\n', results_real(omega_ratios == 0.8));
fprintf('omega/omega_ce = 1.0 时: c*k_parallel/omega_ce = %.6f\n', results_real(omega_ratios == 1.0));

% 保存结果到文件
save('k_solver_scan_results.mat', 'omega_ratios', 'results_real', 'results_imag', 'electron_contributions', 'ion_contributions');
fprintf('结果已保存到 k_solver_scan_results.mat\n');

% 显示物理参数
fprintf('\n=== 物理参数 ===\n');
fprintf('omega_ce = %.2e rad/s\n', omega_ce);
fprintf('omega_pe = %.2e rad/s\n', omega_pe);
fprintf('omega_pi = %.2e rad/s\n', omega_pi);
fprintf('omega_ci = %.2e rad/s\n', omega_ci);
fprintf('w_para_e = %.2e m/s\n', w_para_e);
fprintf('w_para_i = %.2e m/s\n', w_para_i); 


% =====================
% zetaf函数 - 等离子体色散函数
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
n=20+max(max(ceil(abs(x1/a))));% n- truncation parameter
%n=20+ceil(abs(x/a));
s1=0;s2=0;s3=0;s4=0;s5=0;
for k=1:n
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