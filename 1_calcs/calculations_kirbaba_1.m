%% initial data
f = 50; % Hz
u_1 = 230; % V
u_2 = 10; % V
b_m = 1.5; % T
j = 2.5 * 10^6; % A/m^2
k_fe = 0.9;
k_cu = 0.25;
delta = 0.05 * 10^(-3); % m

% core sizes in meters
H = 38 * 10^(-3);
L = 44 * 10^(-3); 
a = 12 * 10^(-3);
b = 8 * 10^(-3);
c = 12 * 10^(-3);
h = 22 * 10^(-3);

%% task 1

% necessary geometric distances in meters
e = (H - h - delta) / 2;
d = (L - a - 2 * b) / 2;
l_e = (a + d) / 2 + b + e; % 

% hole area in m^2
s_b = b * h;

% section areas in m^2
s_a = a * c;
s_d = d * c;
s_e = e * c;

% nominal current
i_n = s_b * k_cu * j / 2;

% find the maximum induction in the rod 'b_m_cur' using the output loop 
% under the condition of equality of the no-load current to 40% of the
% nominal current value
b_m_cur = b_m;

while 1
    % calculate inductions in different parts
    b_a = b_m_cur;
    b_d = b_m_cur * s_a * k_fe / (2 * s_d * k_fe);
    b_e = b_m_cur * s_a * k_fe / (2 * s_e * k_fe);

    % obtain the field strength from induction
    h_a = h_b_lin_approx(b_a); 
    h_d = h_b_lin_approx(b_d);
    h_e = h_b_lin_approx(b_e);

    h_oa = b_a / (sqrt(2) * 4 * pi * 10 ^ (-7)); % A/m
    h_od = b_d / (sqrt(2) * 4 * pi * 10 ^ (-7)); % A/m
    
    u_o = (h_oa + h_od) * delta;
    u_a = h_a * h;
    u_d = h_d * h;
    u_e = 2 * h_e * l_e;
    
    i_mu = u_o + u_a + u_d + u_e; % magnetized current
       
    k = i_mu / i_n;

    if  k > 0.4
        b_m_cur = b_m_cur - 0.001; % step
    else
        break
    end
end

%% task 2
s = 2.22 * f * b_m_cur * s_a * k_cu * s_b * k_fe * j;
k_tr = u_1/u_2;
w_1 = u_1 / (4.44 * f * b_m * s_a * k_fe);
w_2 = w_1 / k_tr;

i_1n = i_n/w_1;
i_2n = i_n/w_2;

s_10 = i_1n / j;
s_20 = i_2n / j;

rho_1 = rho_s_lin_approx(s_10);
rho_2 = rho_s_lin_approx(s_20);

l_1 = 2 * (a + c + 3*b);
l_2 = 2 * (a + c + b);

r_1 = l_1 * w_1 * rho_1;
r_2 = l_2 * w_2 * rho_2;

%% task 3

gamma_fe = 7800;
v_a = a * h * c;
v_d = d * h * c;
v_e = e * c * (a + 2*b + 2*d);

g_a = v_a * gamma_fe * k_fe;
g_d = v_d * gamma_fe * k_fe;
g_e = v_e * gamma_fe * k_fe;

p_a = p_b_lin_approx(b_a);
q_a = q_b_lin_approx(b_a);
p_d = p_b_lin_approx(b_d);
q_d = q_b_lin_approx(b_d);
p_e = p_b_lin_approx(b_e);
q_e = q_b_lin_approx(b_e);

p_fe = g_a * p_a + 2 * g_d * p_d + 2 * g_e * p_e;
q_fe = g_a * q_a + 2 * g_d * q_d + 2 * g_e * q_e;

p_cu = r_1 * i_1n^2 + r_2 * i_2n^2;

i_o = i_mu / w_1;
p = r_1 * i_o^2 + p_fe;

cosphi = 1 / sqrt(1 + (q_fe/p)^2);
beta_max = sqrt(p_fe / q_fe);
eta_n = u_1 * i_1n / (u_1 * i_1n + p_fe + p_cu);
eta_max = (beta_max * u_1 * i_1n) / (beta_max * u_1 * i_1n + p_fe + p_cu * beta_max^2);

%% functions

% find the dependence H(B) from the tabulated data
function h = h_b_lin_approx(b) 
    b_data = [0 0.5 1 1.25 1.5 1.6 1.7 1.8 1.9 2];
    h_data = [0 0.3 0.6 1.6 4.8 7.2 8.6 14.2 24 40] .* 100;
    
    if b < 0
        h = 0;
        return
    end
    
    l_index = 0;
    for index = 1:1:length(b_data)
        if b > b_data(index)
            l_index = index;
        end
    end

    k = polyfit([b_data(l_index), b_data(l_index + 1)], [h_data(l_index), h_data(l_index + 1)], 1);
    h = b .* k(1) + k(2);  
end

% get the resistivity dependence rho(S) from the cross section of the conductor
function rho = rho_s_lin_approx(s)
    s_data = [0.00368 0.00502 0.00636 0.00785 0.00850 0.01131 0.01327 0.01539 0.01767 0.02011 0.0227 0.02545 0.02835 0.03142 0.03464 0.04155 0.04909 0.05726 0.06605 0.07548 0.08853 0.09621 0.11341 0.13202 0.15205 0.17349 0.18848 0.20428 0.22051 0.23578 0.25565 0.27340 0.30191 0.32170 0.35256 0.37393 0.40715 0.43008 0.46556 0.50265 0.54060 0.58088 0.63617 0.67920 0.72382 0.78540 0.84950 0.91610 0.98520 1.0568 1.1310 1.2272 1.3273 1.4314 1.5394 1.6513 1.7670 1.9113 2.0612 2.2167 2.3780] .* 10.^(-6);

    rho_data = [4.4 3.63 2.86 2.24 1.85 1.55 1.32 1.14 0.994 0.873 0.773 0.688 0.618 0.558 0.507 0.423 0.357 0.306 0.266 0.233 0.205 0.182 0.155 0.133 0.115 0.101 0.0931 0.0859 0.0793 0.0739 0.0687 0.0643 0.0579 0.0546 0.0497 0.0469 0.0430 0.0408 0.0376 0.0349 0.0324 0.0302 0.0275 0.0258 0.0242 0.0224 0.0206 0.0192 0.0177 0.0166 0.0155 0.0143 0.0132 0.0122 0.0114 0.0106 0.00989 0.00918 0.00850 0.00792 0.00736];

    if s < 0
        rho = 0;
        return
    end
    
    l_index = 0;
    for index = 1:1:length(s_data)
        if s > s_data(index)
            l_index = index;
        end
    end

    k = polyfit([s_data(l_index), s_data(l_index + 1)], [rho_data(l_index), rho_data(l_index + 1)], 1);
    rho = s .* k(1) + k(2);  
end

% get the specific reactive power q(b) by induction
function q = q_b_lin_approx(b)
    b_data = [0 0.5 1 1.25 1.5 1.6 1.7 1.8 1.9 2];
    q_data = [0 0.43 1.7 5.7 20.5 32.5 41.5 73 110 180];

    if b < 0
        q = 0;
        return
    end
    
    l_index = 0;
    for index = 1:1:length(b_data)
        if b > b_data(index)
            l_index = index;
        end
    end

    k = polyfit([b_data(l_index), b_data(l_index + 1)], [q_data(l_index), q_data(l_index + 1)], 1);
    q = b .* k(1) + k(2);  
end

% get specific power loss in steel p(b) from induction
function p = p_b_lin_approx(b)
    b_data = [0 0.5 1 1.25 1.5 1.6 1.7 1.8 1.9 2];
    p_data = [0 0.1 0.5 0.9 1.4 1.7 2.0 2.45 3 4];

    if b < 0
        p = 0;
        return
    end
    
    l_index = 0;
    for index = 1:1:length(b_data)
        if b > b_data(index)
            l_index = index;
        end
    end

    k = polyfit([b_data(l_index), b_data(l_index + 1)], [p_data(l_index), p_data(l_index + 1)], 1);
    p = b .* k(1) + k(2);  
end
