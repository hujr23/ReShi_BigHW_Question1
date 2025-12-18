function [value, isterminal, direction] = droplet_events(t, Y, params)
% droplet_events
% 统一事件函数：
%   1 - 滑移
%   2 - 脱附
%   3 - 破碎
%
% 任何一个事件触发 → 中止积分

R   = Y(1);
T_w = Y(2); %#ok<NASGU>

%% ========= 统一计算接触角 θ（与你 ODE 完全一致） =========
u        = params.u;
mu_air  = params.mu_air;
sigma   = params.sigma;
theta_e = params.theta_e;

Lm = 1e-9;
A_cv = 9 * params.mu_air * params.u / params.sigma;



fun_theta = @(th) ...
    th.^3 ...
  - theta_e^3 ...
  - A_cv .* log( (2*pi*R.*max(sin(th),1e-6)) / Lm );

% 用平衡接触角作为初值，非常稳
theta = fzero(fun_theta, theta_e);
theta = max(min(theta, 0.99*pi), 1e-2);



%% ========= 几何 =========
Ap = pi * R^2 * sin(theta)^2;
Rc = R * sin(theta);
L  = 2 * Rc;

V = (pi * R^3 / 3) * ...
    (2 - 3*cos(theta) + cos(theta)^3);

%% ========= 事件 1：滑移 =========
rho_air = params.rho_air;
C_d = params.C_d;

F_d = 0.5 * C_d * rho_air * Ap * u^2;
F_tau = mu_air * Ap * (u / Rc);

value_slide = F_tau - F_d;

%% ========= 事件 2：脱附 =========
rho_w = params.rho_w;
g     = params.g;

F_g = rho_w * g * V;
F_sigma = sigma * L * sin(theta);

value_detach = F_sigma - F_g;

%% ========= 事件 3：破碎（Weber 数） =========
We = rho_air * u^2 * (2*R) / sigma;
We_crit = params.We_crit;

value_break = We_crit - We;

%% ========= 汇总 =========
value = [value_slide;
         value_detach;
         value_break];

isterminal = [1; 1; 1];   % 全部终止
direction  = [-1; -1; -1];

end

