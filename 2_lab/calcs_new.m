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
cos_phi = data{:, "Var9"};

%% initial data
f1 = 50;
zp = 1;
m1 = 3;
n0 = 60 * f1 / zp;
omega1 = 2 * pi * f1;
s = (n0 - n) / n0;

n0eta = 2918;
s_eta = (n0eta - n) / n0eta;

%% eta exp calc
eta_exppp_motor = abs(Pmech(s_eta >= 0) ./ P1(s_eta >= 0));
eta_exppp_gen = abs(P1(s_eta < 0) ./ Pmech(s_eta < 0));
eta_exppp = cat(1,eta_exppp_gen, eta_exppp_motor);
for i=find(abs(eta_exppp)>0.8)
    eta_exppp(i) = 0;
end

%% plot ML
plot(ML, n, LineWidth=2)
xlabel('M, Н\cdotм')
ylabel('n, об/мин')
grid on
legend('n(M_L)')

%% create modeling vars
n_modeling = -500:3200;
s_modeling = (n0 - n_modeling) / n0;
seta_modeling = (n0eta - n_modeling) / n0eta;
idle_speed_idx = find(n_modeling == n0);
start_idx = find(n_modeling == 0);

n_is_s = n_modeling(start_idx:idle_speed_idx);
s_is_s = s_modeling(start_idx:idle_speed_idx);

%% split the machine on generator and brake modes
% experiment
n_gen = n(n > 0);
ML_gen = ML(n > 0);
n_br = n(n < 0);
ML_br = ML(n < 0);

% modeling
n_gen_m = n_modeling(n_modeling > 0);
ML_gen_m = interp1(n_gen, ML_gen, n_gen_m, 'makima');
n_br_m = n_modeling(n_modeling < 0);
ML_br_m = linspace(ML_br(end), ML_br(1), length(n_br_m));
%% calculate dry and viscous friction
Mtf_s = abs(ML_gen_m(1) - ML_br_m(end)) / 2; % starting dry moment
M_se = ML_gen_m(1) + Mtf_s; % starting moment

ML_m = [ML_br_m M_se ML_gen_m];

Mtf = -Mtf_s * sign(n_modeling);
Mvf = (ML_m(idle_speed_idx) - Mtf(idle_speed_idx)) * n_modeling / n0;
MIM = ML_m - Mtf - Mvf;

M_motor = MIM(start_idx:idle_speed_idx);

%% plot n(M_IM), n(M_tf), n(M_L), n(M_vf)
grid on
hold on
plot(ML, n, LineWidth=2)
plot(Mtf, n_modeling, MIM, n_modeling, Mvf, n_modeling, LineWidth=2)
xlabel('M, Н\cdotм')
ylabel('n, об/мин')
legend('n(M_L)', 'n(M_{tf})', 'n(M_{IM})', 'n(M_{vf})')

%% MIM in range
MIM_is_s = MIM(idle_speed_idx:start_idx);

%% plot MIM(s)
plot(s_modeling, MIM, LineWidth=2)
xlabel('s')
ylabel("M, Н\cdotм")
grid on
legend("M_{IM}(s)")

%% params calcs
s_m = s_modeling(MIM == max(MIM));
U_calcs = mean(U);

I_motor = interp1(n, I, n_is_s, 'makima');
I_work = interp1(n, I, n_modeling, 'makima');

% входная мощность (считаем для одной фазы)
P1_motor = interp1(n, P1, n_is_s,'makima') / m1;
P1_work = interp1(n, P1, n_modeling,'makima') / m1;
cos_fi_work = P1_work ./ (I_work * U_calcs);
sin_fi_work = sqrt(1 - cos_fi_work.^2);

% холостой ход
r0 = P1_motor(end) / I_motor(end)^2;
x0 = U_calcs * sin_fi_work(end) / I_motor(end);

% режим короткого замыкания [пуск]
x_k = U_calcs * sin_fi_work(1) / I_motor(1);
x_s1 = x_k / 2;
x_s2 =x_k / 2;
x_m = x0 - x_s1;
c1 = 1 + x_s1 / x_m;
r20 = omega1 * M_motor .* s_is_s ./ m1 ./ I_motor.^2 / zp;
delta = (s_is_s(1) - s_is_s(2)) / 2;
expr = abs(s_is_s - s_m);
r2 = r20(expr < delta);
r1 = sqrt((r2 / s_m)^2 - x_k^2);
r_m = r0 - r1;

%% 
P_mech_exp = MIM.*(n_modeling/60*2*pi);
P1_exp = interp1(n,P1,n_modeling,'makima');

eta_motor_exp = abs(P_mech_exp(s_modeling>=0)./P1_exp(s_modeling>=0));
eta_gen_exp = abs(P1_exp(s_modeling<0)./P_mech_exp(s_modeling<0));
eta_exp = cat(2, eta_motor_exp, eta_gen_exp);
for i=find(abs(eta_exp)>0.8)
    eta_exp(i) = 0;
end

%% calcs 2

h0 = 1;
Ms = @(h) Ms_func(h,m1,zp,c1,r2,omega1,r1,x_s1,x_s2,U_calcs,M_se);
h = fzero(Ms,h0)

%% plot
hhh = linspace(0, 5, 10^3);
out = linspace(0, 5, 10^3);
for index = 1:10^3
   out(index) = Ms(hhh(index));
end
plot(hhh, out, LineWidth=2)
yline(0, LineWidth=2)
grid on
xlabel('h')
ylabel("M_s, Н\cdotм")
legend("M_s(h)", "M_s=0")

%% calcs 3
betta = 1.4; %1.3

ksi = h * abs(s_modeling).^betta;
kr = ksi .* (sinh(2*ksi) + sin(2*ksi)) ./ (cosh(2*ksi) - cos(2*ksi));
kx = 3./(2.*ksi).*(sinh(2.*ksi)-sin(2.*ksi))./(cosh(2.*ksi)-cos(2.*ksi));
Z1 = r1 + 1i * x_s1;
Z2 = kr.*c1.*r2./s_modeling + 1i.*c1.*kx.*x_s2;
Zm = r_m +1i*x_m;
Z_in = Z1+Zm.*Z2./(Zm+Z2);
I1_calc = U_calcs ./ Z_in;
I2_calc = I1_calc.*Zm./(Zm+Z2);
cos_fi_calc = cos(angle(Z_in));
P1_calc = m1 * U_calcs .*abs(I1_calc).*cos_fi_calc;
P_mech_calc = m1 *abs(I2_calc).^2.*kr.*c1*r2.*(1-s_modeling) ./ s_modeling;
M_calc = zp * P_mech_calc./ omega1 ./(1-s_modeling);
eta_motor_calc = P_mech_calc./P1_calc;
eta_gen_calc = P1_calc./P_mech_calc;
eta_calc = [eta_motor_calc eta_gen_calc];
for i=find(abs(eta_calc)>1)
    eta_calc(i) = 0;
end 

%% plots
plot(s_modeling,cos_fi_calc,'Color','b','LineWidth',1.5); % экспериментальная
grid on; 
xlabel('s','FontSize',16,'FontAngle','italic');
ylabel("cos\phi",'FontSize',16,'FontAngle','italic');
set(gca,'FontSize',14);

%% plots 2
plot(s_modeling,P1_calc,'Color','b','LineWidth',1.5); % экспериментальная
grid on; 
xlabel('s','FontSize',16,'FontAngle','italic');
ylabel("P_{1}",'FontSize',16,'FontAngle','italic');
set(gca,'FontSize',14);

%% plots 3
plot(ML,n,'Color','b','LineWidth',1.5); % экспериментальная
grid on; 
% title('экспериментальная характеристика','FontSize',20,'FontWeight','bold');
hold on;
plot(MIM,n_modeling, 'Color','g','LineWidth',1.5); % идеальная (без потерь)
plot(M_calc,n_modeling, 'Color', 'r', 'LineStyle','--','LineWidth',1.5);
hold off;
xlabel('M, Н\cdotм','FontSize',16,'FontAngle','italic');
ylabel("n, об/мин",'FontSize',16,'FontAngle','italic');
set(gca,'FontSize',14);
legend('экспериментальная','"истинная"','рассчитанная','Location','southwest');

%% asdf
sin_fi_work(1)
I_motor(1)
P1_motor(1)

%% asdadasad

plot(s_is_s, r20, 'LineWidth', 2); % экспериментальная
grid on; 
xlabel('s');
ylabel("r^'_2");
legend("r^'_2(s)")

%% выфа
plot(s,cos_phi,'LineWidth',2); % экспериментальная
grid on;
hold on;
plot(s_modeling,cos_fi_work,'LineWidth',2); % идеальная (без потерь)
plot(s_modeling,cos_fi_calc, 'LineWidth',2);
hold off;
xlabel('s');
ylabel("cos phi");
legend('experiment','corrected','modeling');

%% pogreshnost
new = interp1(s, cos_phi, s_modeling, 'makima');
plot(s_modeling, cos_fi_calc-new,'LineWidth',2);
grid on;
xlabel('s');
ylabel("error")
legend("modeling - experiment")
%% выфа
plot(s,I,'LineWidth',2); % экспериментальная
grid on;
hold on;
plot(s_modeling,I_work,'LineWidth',2); % идеальная (без потерь)
plot(s_modeling,abs(I2_calc), 'LineWidth',2);
hold off;
xlabel('s');
ylabel("I^'_2");
legend('experiment','corrected','modeling');

%% pogreshnost
new = interp1(s, I, s_modeling, 'makima');
plot(s_modeling, abs(I2_calc)-new,'LineWidth',2);
grid on;
xlabel('s');
ylabel("error")
legend("modeling - experiment")

%% выфа
plot(s,P1,'LineWidth',2); % экспериментальная
grid on;
hold on;
plot(s_modeling,P1_exp,'LineWidth',2); % идеальная (без потерь)
plot(s_modeling,P1_calc, 'LineWidth',2);
hold off;
xlabel('s');
ylabel("P_1");
legend('experiment','corrected','modeling');

%% pogreshnost
new = interp1(s, P1, s_modeling, 'makima');
plot(s_modeling, P1_calc-new,'LineWidth',2);
grid on;
xlabel('s');
ylabel("error")
legend("modeling - experiment")

%% фываыфвафыва
plot(s,Pmech,'LineWidth',2); % экспериментальная
grid on;
hold on;
plot(s_modeling,P_mech_exp,'LineWidth',2); % идеальная (без потерь)
plot(s_modeling,P_mech_calc, 'LineWidth',2);
hold off;
xlabel('s');
ylabel("P_{mech}");
legend('experiment','corrected','modeling');
%% pogreshnost
new = interp1(s, Pmech, s_modeling, 'makima');
plot(s_modeling, P_mech_calc-new,'LineWidth',2);
grid on;
xlabel('s');
ylabel("error")
legend("modeling - experiment")

%% фываыфвафыва
plot(s,ML,'LineWidth',2); % экспериментальная
grid on;
hold on;
plot(s_modeling,MIM,'LineWidth',2); % идеальная (без потерь)
plot(s_modeling,M_calc, 'LineWidth',2);
hold off;
xlabel('s');
ylabel("M, Н\cdotм");
legend('experiment','corrected','modeling');

%% pogreshnost
new = interp1(s, ML, s_modeling, 'makima');
plot(s_modeling, M_calc-new,'LineWidth',2);
grid on;
xlabel('s');
ylabel("error")
legend("modeling - experiment")

%% sdaffasd
qqqqq = interp1(s, eta_exppp, s_modeling, 'makima');

eta_exxx = 0.6 * qqqqq + 0.4 .* eta_motor_calc;

%% фываыфвафывафыв
plot(s,eta_exppp,'LineWidth',2); % экспериментальная
grid on;
hold on;
%plot(s_modeling,MIM,'LineWidth',2); % идеальная (без потерь)
%plot(s_modeling - 0.045,eta_exxx, 'LineWidth',2);
plot(s_modeling,eta_motor_calc, 'LineWidth',2);
ylim([0, 1])
hold off;
xlabel('s');
ylabel("\eta");
legend('experiment','modeling');

%% pogreshnost
new = interp1(P1, Pmech, P1_calc, 'makima');
plot(P_mech_calc-new, P1_calc,'LineWidth',2);
grid on;
xlabel('error');
ylabel("P_1")
legend("modeling - experiment")

%% sdflk
plot(Pmech, eta_exppp,'LineWidth',2); % экспериментальная
grid on;
hold on;
ylim(([0 1]))
%plot(P_mech_exp, eta_exp,'LineWidth',2); % идеальная (без потерь)
plot(P_mech_calc, eta_motor_calc, 'LineWidth',2);
hold off;
xlabel('P_{mech}');
ylabel("\eta");
legend('experiment','modeling');