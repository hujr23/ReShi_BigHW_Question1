function dYdt = droplet_ode(t, Y , params)
% droplet_ode
% 水滴生长的微分方程组
%
% 输入：
%   t      时间 (s)
%   Y(1)   水滴曲率半径 R (m)
%   Y(2)   水滴温度 Tw (℃)
%   params 结构体，存放所有物性和工况参数
%
% 输出：
%   dYdt(1) = dR/dt
%   dYdt(2) = dTw/dt

% ========= 取出状态量 =========
R  = Y(1);

R = max(R, 1e-7);

T_w = Y(2);

% ========= 取出参数 =========
Ta   = params.T_air;
Tsur = params.T_sur;
u    = params.u;
rho_air = params.rho_air;
mu_air  = params.mu_air;
sigma   = params.sigma;
k_a     = params.k_a;
Pr      = params.Pr;
D_aw    = params.D_aw;
Pr = params.Pr;
RH = params.RH;
Rgas = params.R_gas;

rho_w = params.rho_w;   % 水密度
g     = params.g;       % 重力加速度
C_d   = params.C_d;









% ======== 变量关联物理公式 ========

%% --- 动态接触角（Cox–Voinov 模型） ---

theta_e = pi/5;    % 平衡接触角（弧度）
Lm = 1e-9;         % 分子尺度长度 (m)

A_cv = 9 * mu_air * u / sigma;

% 固定 sin(theta) = 1/sqrt(2)
sin_theta_eff = 1/sqrt(2);

theta_cubed = theta_e^3 ...
    + A_cv * log( (2*pi*R*sin_theta_eff) / Lm );

% 数值保护
theta_cubed = max(theta_cubed, 1e-6);

theta = theta_cubed^(1/3);

% 再加一道物理限幅（非常重要）
theta = max(min(theta, 0.99*pi), 1e-2);


%% --- 球形帽状水滴的几何关系 ---

% 体积
V = (pi * R^3 / 3) * ( ...
      2 - 3*cos(theta) + cos(theta)^3 );

% 表面积
A = 2 * pi * R^2 * (1 - cos(theta));

% 投影面积
Ap = pi * R^2 * sin(theta)^2;

% 接触半径
Rc = R * sin(theta);

% 湿润长度（按给定公式）
L = 2 * R * sin(theta);





%% --- 雷诺数 Re ---
D = 2 * Rc;   % 特征长度：水滴直径
Re = rho_air * u * D / mu_air;



%% --- 努塞尔数 Nu ---
Nu = 2 + 0.6 * Re^0.5 * Pr^(1/3);


%% --- 对流换热系数 h_c ---
h_c = Nu * k_a / D;



%% --- 施密特数 Sc ---
Sc = mu_air / (rho_air * D_aw);




%% --- 舍伍德数 Sh ---
Sh = 2 + 0.6 * Re^0.5 * Sc^(1/3);



%% --- 对流传质系数 h_m ---
h_m = Sh * D_aw / D;



%% --- 韦伯数 We ---
We = rho_air * u^2 * D / sigma;




%% --- 水滴表面饱和水蒸气压 P_sat(T_w) ---

% 水滴温度（K）
Tw_K = T_w + 273.15;
% 空气温度（K）
Ta_K = Ta + 273.15;



% ASHRAE 饱和水蒸气压公式常数
c1 = -5.8002206e3;
c2 =  1.3914993;
c3 = -4.8640239e-2;
c4 =  4.1764768e-5;
c5 = -1.4452093e-8;
c6 =  6.5459673;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% ln(P_sat) 计算
ln_Psat_air = c1 ./ Ta_K ...
            + c2 ...
            + c3 .* Ta_K ...
            + c4 .* Ta_K.^2 ...
            + c5 .* Ta_K.^3 ...
            + c6 .* log(Ta_K);
P_sat = exp(ln_Psat_air);

% 空气侧水蒸气分压
P_v = RH * P_sat;

P_s = 101325;
X_w = P_v / P_s;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% --- 水蒸气摩尔通量 N_w ---



Tw = Y(2);                % ℃
P_sat_w = 610.78 * exp(17.27 * Tw / (Tw + 237.3));  % Pa



N_w = h_m * ( ...
      P_v / (Rgas * Ta_K) ...
    - P_sat_w / (Rgas * Tw_K) );




%%%%%%%%%%%下面是原来的，可能改错了，这里没有删除！！！

%N_w = h_m * ( ...
 %     P_sat / (Rgas * Tw_K) ...
 %   - X_w  * P_s / (Rgas * Ta_K) );
%
%% --- 水的相变潜热 h_fg(T_w) ---

% h_fg: kJ/kg → J/kg
h_fg = (2500 - 2.35 * T_w) *1e3;








%% --- 水滴所受力（力学模块） ---

% 重力（方向：竖直向下）
F_g = rho_w * g * V;

% 表面张力（沿接触线的附着力）
F_sigma = sigma * L * sin(theta);

% 气流拖拽力（沿气流方向）
F_d = 0.5 * C_d * rho_air * Ap * u^2;



% ======== 微分方程（ODE 主体） ========


%% --- 水滴质量 ---
m_w = rho_w * V;

%% --- 微分方程 1：水滴半径变化 dR/dt ---

dRdt =  N_w * params.M_w / ...
       ( rho_w ...
       * (2 - 3*cos(theta) + cos(theta)^3) ...
       * sin(theta)^2 );
% 防止数值爆炸
if ~isfinite(dRdt)
    dRdt = 0;
end



%% --- 微分方程 2：水滴温度变化 dTw/dt ---



k_w = 0.606;          % 水导热系数 (W/m·K)
delta_eff = 1e-5;     % 等效液膜厚度 (m)，1~10 μm 都合理

A_c = pi * Rc^2;      % 接触面积

% 防止 Rc → 0
A_c = max(A_c, 1e-12);

R_ws = delta_eff / (k_w * A_c);

Q_ws = (T_w - Tsur) / R_ws;   % W，正值表示向翅片放热







% 空气-水滴对流换热
Q_aw = h_c * A * (Ta - T_w);

% 冷凝潜热（正号：冷凝放热）
Q_lat = N_w * Ap * params.M_w * h_fg;

% 总能量平衡
dTwdt = ( Q_aw ...
        - Q_ws ...      % 向翅片导热（冷却项）
        + Q_lat ...     % 冷凝放热
        ) / ( m_w * params.Cp_w );


% ========= 输出 =========
dYdt = [dRdt; dTwdt];
end
