%% 主脚本：水滴生长基本参数定义

T_air = 28;         % 空气温度 (℃)
RH = 0.85;        % 相对湿度（85%）

rho_air = 1.18;    % 空气密度 (kg/m^3)
mu_air = 1.85e-5;   % 空气动力黏度 (Pa·s)
sigma = 72.00e-3;   %空气与水滴表面的张力系数，近似为常数，(N/m)
C_d = 0.0063 ; %气流拖拽系数，在u<4.0m/s时为常数，是无因次数
k_a = 2.53e-2 ;%空气的导热系数，W/m*K
k_w = 0.606 ; %水的导热系数
Cp_air = 1005 ;%空气比热
Cp_w = 4186 ;%水的比热
D_aw = 2.6e-5;%空气中水蒸气扩散系数,m^2/s
M_w = 18.01528e-3 ;%水的摩尔质量，kg/mol
delta_eff = 1e-5;     % 等效液膜厚度 (m)

P = 101325;%大气压，Pa


D0 = 1.0e-4;          %水滴初始直径、半径(m)
R0 = D0/2;

R_gas = 8.314; % 通用气体常数，J/(mol·K)

rho_w = 997;    % kg/m^3
g = 9.81;        % m/s^2

Pr = 0.707 ; %空气的普朗特数，近似为常数





T_w0 = 14;   % 水滴初始温度 (℃)，接近冷凝表面温度


%% ========= 工况参数 =========
params.T_air = T_air;
params.RH = RH;

params.rho_air = rho_air;
params.mu_air = mu_air;
params.sigma = sigma;
params.k_a = k_a;
params.Pr = Pr;
params.D_aw = D_aw;
params.M_w = M_w;
params.R_gas = R_gas;

params.rho_w = rho_w;
params.g = g;
params.C_d = C_d;
params.Cp_w = Cp_w;

params.T_sur = 12;


%% ========= 来流速度参数组 =========
u_list = [0.5, 0.8, 1.2, 1.6, 2.0];   % m/s


params.theta_e = pi/5;   % 平衡接触角（弧度）
params.We_crit = 10;   % 临界 Weber 数（课程作业级完全够）















%% ========= ODE 初始条件 =========
Y0 = [R0; T_w0];   % [初始半径; 初始水滴温度]

%% ========= 时间区间 =========
tspan = [0 20000];    % 秒（会被事件函数提前终止）


%% ========= 作图窗口初始化 =========
figure_R = figure; hold on; grid on;
xlabel('t (s)'); ylabel('R (mm)');
title('不同来流速度下水滴半径演化');

figure_V = figure; hold on; grid on;
xlabel('t (s)'); ylabel('V (mm^3)');
title('不同来流速度下水滴体积演化');

figure_T = figure; hold on; grid on;
xlabel('t (s)'); ylabel('T_w (℃)');
title('不同来流速度下水滴温度演化');




%% ========= 不同来流速度下的 ODE 求解 =========

for iu = 1:length(u_list)

    params.u = u_list(iu);

    options = odeset( ...
        'Events', @(t,Y) droplet_events(t, Y, params), ...
        'RelTol',1e-6, ...
        'AbsTol',1e-9, ...
        'MaxStep',1.0 );

    [t, Y] = ode15s(@(t,Y) droplet_ode(t, Y, params), ...
                    tspan, Y0, options);


 % ======== ★ 就插在这里 ★ ========
    if ~isempty(ie)
        switch ie(1)
            case 1
                fprintf('u = %.1f m/s：水滴发生滑移，t = %.2f s\n', params.u, te(1));
            case 2
                fprintf('u = %.1f m/s：水滴发生脱附，t = %.2f s\n', params.u, te(1));
            case 3
                fprintf('u = %.1f m/s：水滴发生破碎，t = %.2f s\n', params.u, te(1));
        end
    else
        fprintf('u = %.1f m/s：时间到达上限，水滴仍处于生长阶段\n', params.u);
    end










    R  = Y(:,1);
    Tw = Y(:,2);
    Nt = length(t);
    V  = zeros(Nt,1);

    u        = params.u;
    mu_air  = params.mu_air;
    sigma   = params.sigma;
    theta_e = params.theta_e;
    Lm = 1e-9;

    for k = 1:Nt
        R_k = max(R(k),1e-8);
        A_cv = 9 * mu_air * u / sigma;

        fun_theta = @(th) ...
            th.^3 ...
          - theta_e^3 ...
          - A_cv .* log( (2*pi*R_k.*max(sin(th),1e-6)) / Lm );

        theta_k = fzero(fun_theta, theta_e);

        V(k) = (pi * R_k^3 / 3) * ...
               (2 - 3*cos(theta_k) + cos(theta_k)^3);
    end
    figure(figure_R);
    plot(t, R*1e3, 'LineWidth',1.5, ...
        'DisplayName', sprintf('u = %.1f m/s', params.u));

    figure(figure_T);
    plot(t, Tw, 'LineWidth',1.5, ...
        'DisplayName', sprintf('u = %.1f m/s', params.u));

    figure(figure_V);
    plot(t, V*1e9, 'LineWidth',1.5, ...
        'DisplayName', sprintf('u = %.1f m/s', params.u));

end
figure(figure_R); legend show;
figure(figure_T); legend show;
figure(figure_V); legend show;
