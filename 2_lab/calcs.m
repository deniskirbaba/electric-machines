%% load data
data = readtable('values.csv');
n = round(data{:, "Var1"});
ML = data{:, "Var2"};
Pmech = data{:, "Var3"};
U = data{:, "Var4"} / sqrt(3);
I = data{:, "Var5"};
S = data{:, "Var6"};
P1 = data{:, "Var7"};
Q = data{:, "Var8"};
cos = data{:, "Var9"};

%% initial data
f1 = 50;
zp = 1;
m1 = 3;
n0 = 60 * f1 / zp;
omega1 = 2 * pi * f1;
s = (n0 - n) / n0;
eta = Pmech ./ P1;

%% plot ML
plot(ML, n, LineWidth=2)
xlabel('M, Н\cdotм')
ylabel('n, об/мин')
grid on
legend('n(M_L)')

%% create modeling vars
n_m = -500:3200;
s_m = (n0 - n_m) / n0;
idle_speed_idx = find(n_m == n0);
start_idx = find(n_m == 0);

n_is_s = n_m(start_idx:idle_speed_idx);
s_is_s = s_m(start_idx:idle_speed_idx);
%% split the machine on generator and brake modes
% experiment
n_gen = n(n > 0);
ML_gen = ML(n > 0);
n_br = n(n < 0);
ML_br = ML(n < 0);

% modeling
n_gen_m = n_m(n_m > 0);
ML_gen_m = interp1(n_gen, ML_gen, n_gen_m, 'makima');
n_br_m = n_m(n_m < 0);
ML_br_m = linspace(ML_br(end), ML_br(1), length(n_br_m));
%% calculate dry and viscous friction
Mtf_s = abs(ML_gen_m(1) - ML_br_m(end)) / 2; % starting dry moment
M_se = ML_gen_m(1) + Mtf_s; % starting moment

ML_m = [ML_br_m M_se ML_gen_m];

Mtf = -Mtf_s * sign(n_m);
Mvf = (ML_m(idle_speed_idx) - Mtf(idle_speed_idx)) * n_m / n0;
MIM = ML_m - Mtf - Mvf;

%% plot n(M_IM), n(M_tf), n(M_L), n(M_vf)
grid on
hold on
plot(ML, n, LineWidth=2)
plot(Mtf, n_m, MIM, n_m, Mvf, n_m, LineWidth=2)
xlabel('M, Н\cdotм')
ylabel('n, об/мин')
legend('n(M_L)', 'n(M_{tf})', 'n(M_{IM})', 'n(M_{vf})')

%% MIM in range
MIM_is_s = MIM(idle_speed_idx:start_idx);

%% calculate mech power and eta
Pmech_true = MIM .* (n_m / 60 * 2 * pi);
P1 = interp1(n, P1, n_m, 'makima');
eta_motor_ture = abs(Pmech_true(s_m >= 0) ./ P1(s_m >= 0));
eta_gen_true = abs(P1(s_m < 0) ./ Pmech_true(s_m < 0));
eta_true = [eta_motor_ture eta_gen_true];
for idx = find(abs(eta_true) > 0.8)
    eta_true(idx) = 0;
end

%% plot mech power
plot(s, Pmech, s_m, Pmech_true, LineWidth=2)
grid on
xlabel('s')
ylabel('P, Вт')
legend('experiment', 'true')

%% plot eta
plot(s, eta, s_m, eta_true, LineWidth=2)
grid on
xlabel('s')
ylabel('\eta')
legend('experiment', 'true')
%% parameters calculation
P2_data_real = (MIM .* 2 .* pi .* n) ./ 60;
P1n0 = 89.3782;
I1n0 = 0.19377;
r0 = P1n0 / I1n0^2;

U1n0 = 229.186;
cosn0 = 0.621281;
x0 = U1n0 * sqrt(1 - cosn0^2) / I1n0; 

I10 = I(30);
P10 = P1(30);
x1s = (U1n0 * sqrt(1 - cos_data(30))) / (2 * I10);
xk = 2 * x1s;
x2ss = x1s;

x_m = x0 - x1s;
c1 = 1 + x1s / x_m;

%% plot r2s(s)

s_m = 0.410667;
r2s_data = MIM .* sss_data * (2 * pi * f1) ./ ((I) .^ 2 * m1 * zp);

plot(sss_data, r2s_data, LineWidth=2)
xlim([0 1])
xlabel('s')
ylabel("r^'_2, Ом")
grid on
legend("r^'_2(s)")

%% plot M_L(s)
plot(sss_data, MIM, LineWidth=2)
%xlim([0 1])
xlabel('s')
ylabel("M_L, Н\cdotм")
grid on
legend("M_L(s)")

%% more calcs
r2s = 47.982;
r1 = sqrt((r2s / s_m)^2 - xk^2);
rm = r0 - r1;

%% plot Ms(h)

m_st = 1.12;


h_data = linspace(0, 5, 10^6);

uuu = 395 / sqrt(3);
kr = h_data .* (sinh(2 .* h_data) + sin(2 .* h_data)) ./ (cosh(2 .* h_data) - cos(2 .* h_data));
kx = 3 ./ (2 .* h_data) .* (sinh(2 .* h_data) - sin(2 .*  h_data)) ./ (cosh(2 .* h_data) - cos(2 .* h_data));
ms_data = zp .* m1 .* (uuu)^2 .* r2s .* kr ./ (2 .* pi .* f1 .* ((r1 + c1 .* r2s .* kr).^2 + (x1s + c1 .* x2ss .* kx).^2));

plot(h_data, ms_data, LineWidth=2)
yline(m_st, LineWidth=2)
xlabel('h')
ylabel("M, Н\cdotм")
grid on
hold on
legend("M_s(h)", "M_{st}")

%% h calculation

plot(h_data, ms_data - m_st, LineWidth=2)
yline(0, LineWidth=2)
xlabel('h')
ylabel("M, Н\cdotм")
grid on
hold on
legend("M_s(h)", "M = 0")

h = 1.90417;

%% characteristics calcs
beta = 1;
krs = (h .* abs(sss_data).^beta) .* (sinh(2 .* (h .* abs(sss_data).^beta)) + sin(2 .* (h .* abs(sss_data).^beta))) ./ (cosh(2 .* (h .* abs(sss_data).^beta)) - cos(2 .* (h .* abs(sss_data).^beta)));
kxs = 3 ./ (2 .* (h .* abs(sss_data).^beta)) .* (sinh(2 .* (h .* abs(sss_data).^beta)) - sin(2 .*  (h .* abs(sss_data).^beta))) ./ (cosh(2 .* (h .* abs(sss_data).^beta)) - cos(2 .* (h .* abs(sss_data).^beta)));

Z_1 = r1 + 1i * x1s;
Z_2 = krs .* r2s ./ sss_data + 1i .* kxs .* x2ss;
Z_m = rm + 1i * x_m;
Z_in = Z_1 + (Z_m .* Z_2) ./ (Z_m + Z_2);
imzin = imag(Z_in);
reimzin = real(Z_in);
cosphis = atan(imzin./reimzin);

%% cosphi

plot(sss_data, cos_data, sss_data, cosphis, LineWidth=2)
xlabel('s')
ylabel("cosphi")
grid on
legend("experiment", "modeling")

%% I2(s)
I1s = U1_data ./ Z_in;
I2ss = (I1s .* Z_m) ./ (Z_m + Z_2);
plot(sss_data, I, sss_data, abs(I1s), LineWidth=2)
xlabel('s')
ylabel("I, A")
grid on
legend("experiment", "modeling")

%% P1(s)

P1_mod_data = m1 * U1_data .* abs(I1s) .* cosphis;

plot(sss_data, P1, sss_data, P1_mod_data, LineWidth=2)
xlabel('s')
ylabel("P_1, Вт")
grid on
legend("experiment", "modeling")

%% 