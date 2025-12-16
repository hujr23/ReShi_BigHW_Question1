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
params.u = 0.8;

params.theta_e = pi/5;   % 平衡接触角（弧度）
params.We_crit = 10;   % 临界 Weber 数（课程作业级完全够）





%% ========= 失稳模式相图 =========




%% ========= 相图参数范围（密集） =========
u_list    = linspace(0.5, 2.0, 25);    % 25 个速度点
Tsur_list = linspace(8, 15, 25);       % 25 个表面温度点

Nu = length(u_list);
Nt = length(Tsur_list);

mode_map  = zeros(Nt, Nu);
tcrit_map = NaN(Nt, Nu);



%% ========= 构建失稳模式相图 =========
for it = 1:Nt
    for iu = 1:Nu

        params.T_sur = Tsur_list(it);
        params.u     = u_list(iu);

        Y0 = [R0; T_w0];
        tspan = [0 600];

        options = odeset( ...
            'Events', @(t,Y) droplet_events(t, Y, params), ...
            'RelTol',1e-6, ...
            'AbsTol',1e-9, ...
            'MaxStep',1.0 );

        [~, ~, te, ~, ie] = ode15s( ...
            @(t,Y) droplet_ode(t, Y, params), ...
            tspan, Y0, options);

        if ~isempty(ie)
            mode_map(it, iu)  = ie(1);   % 1滑移 2脱附 3破碎
            tcrit_map(it, iu) = te(1);
        else
            mode_map(it, iu)  = 0;
        end

    end
end

figure;
imagesc(u_list, Tsur_list, mode_map);
set(gca,'YDir','normal');





colormap([ ...
    0 0 1;    % 1 滑移 - 蓝
    1 0 0;    % 2 脱附 - 红
    0.6 0.6 0.6]);  % 3 破碎 - 灰

clim([1 3]);
colorbar('Ticks',[1 2 3], ...
         'TickLabels',{'滑移','脱附','破碎'});

xlabel('来流速度 u (m/s)');
ylabel('翅片表面温度 T_{sur} (℃)');
title('水滴失稳模式相图（密集扫描）');
grid on;


hold on;
contour(u_list, Tsur_list, mode_map, ...
        [1.5 1.5], 'w--', 'LineWidth',2);



t_levels = [100 150 200 250 300 350 400];  % 秒，可按你结果改


[C,h] = contour(u_list, Tsur_list, tcrit_map, ...
                t_levels, ...
                'w', 'LineWidth',1.2);

clabel(C, h, 'FontSize',9, 'Color','w');


