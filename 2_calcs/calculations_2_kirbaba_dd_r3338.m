%% initial data

f = 50; % Hz
p2 = 2;
m = 3; % phase number

P_N = 2.2 * 10^3; % W
n_N = 2850; 
U_N = 380;
cosphi_1N = 0.85;
I_N = 4.9; % A (V = 380 V)
M_N = 7.4;

k_s = 2.7;
k_sI = 6.5;
lambda = 2.8;
n_1 = 3000;

%% calculation of missing values

U_1N = U_N / sqrt(3);
I_1N = I_N;
omega_1 = 2 * pi * f;
z_p = p2 / 2;
s_N = 1 - n_N / n_1;

%% parameters calculation

% calculate initial active resistances
r_10 = (U_1N * I_1N * cosphi_1N - M_N * omega_1 / (z_p * m)) / I_1N^2;
r_20_s = m * z_p * U_1N^2 * s_N / (omega_1 * M_N);

% initialization
TOL = 0.001;
Delta_k = 10^(-4);
N = 10^(6);
k_1 = 1;
k_2 = 1;
mu_m = 1;
delta_mu = 1;
n = 0;

while (delta_mu > TOL) && (n < N)
    k_2 = k_2 - Delta_k;
    delta_n = 1;
    mu_m = 1;
    n = n + 1;
    while (delta_n > TOL) && (n < N)
        k_1 = k_1 - Delta_k;
        n = n + 1;
        % calculate rest paratemers (1a)
        mu_s = mu_m;
        a = k_1 * r_10 / (k_2 * r_20_s);
        A = 1 - 2 * a * s_N * (lambda - 1);
        s_m = s_N * (lambda + sqrt(lambda^2 - A)) / A;
        x_ks = sqrt((k_2 * r_20_s / s_m)^2 - (k_1 * r_10)^2);
        x_s2_s = x_ks / 2;
        x_s1 = x_s2_s;
        b = x_ks / ((k_1 * r_10 + k_2 * r_20_s / s_N)^2 + x_ks^2);
        x_m = 1 / (I_1N * sqrt(1 - cosphi_1N^2)/U_1N - b);
        I_2_s = U_1N / sqrt((k_1 * r_10 + k_2 * r_20_s/s_m)^2 + (x_s1 + x_s2_s)^2);
        mu_m = m * z_p * I_2_s^2 * k_2 * r_20_s / (omega_1 * s_m * M_N);
        delta_n = abs(mu_m - mu_s) / mu_s;
        delta_mu = abs(mu_m - lambda) / lambda;
    end
end

r_1 = k_1 * r_10;
r_2_s = k_2 * r_20_s;

%% calculation of the coefficients for the slot depth

h_0 = 0.5; % for low-powered machine
% h_0 = 3; % for high-powered machine

h = h_0;
mu_s = 0;
n = 0;
delta_s = 1;
Delta_h = Delta_k;

while (delta_s > TOL) && (n < N)
    h = h + Delta_h;
    n = n + 1;
    k_x = 3 / (2 * h) * (sinh(2*h) - sin(2*h))/(cosh(2*h)-cos(2*h));
    k_r = h * (sinh(2*h) + sin(2*h))/(cosh(2*h)-cos(2*h));
    mu_s = m * z_p * U_1N^2 * r_2_s * k_r / (omega_1 * ((r_1 + r_2_s * k_r)^2 +(x_s1 + x_s2_s * k_x)^2) * M_N);
    delta_s = abs(mu_s - k_s) / k_s;
end

%% plots of the mechanical characteristic

% form points arrays
s_data = linspace(-1, 1, 10^6);
n_data = (1 - s_data) .* n_1;

% without displacement current
M = m .* z_p .* U_1N.^2 .* r_2_s ./ (omega_1 .* s_data .* ((r_1 + r_2_s ./ s_data).^2 + (x_s1 + x_s2_s).^2));

% with displacement current
h_data = h .* abs(s_data);
k_r_data = h_data .* (sinh(2 .* h_data) + sin(2 .* h_data)) ./ (cosh(2 .* h_data) - cos(2 .* h_data));
k_x_data = 3 ./ (2 .* h_data) .* (sinh(2 .* h_data) - sin(2 .* h_data)) ./ (cosh(2 .* h_data) - cos(2 .* h_data));
M_k = m .* z_p .* U_1N.^2 .* r_2_s .* k_r_data ./ (omega_1 .* s_data .* ((r_1 + k_r_data .* r_2_s ./ s_data).^2 + (x_s1 + x_s2_s .* k_x_data).^2));

%% plot M(s) and M_k(s)
plot(s_data, M, '-', s_data, M_k, '--', 'LineWidth', 3)
grid on
xlabel('s')
ylabel('M, H \cdot м')
legend('M(s)', 'M_k(s)')

%% plot n(M) and n(M_k)
plot(M, n_data, '-', M_k, n_data, '--', 'LineWidth', 3)
grid on
xlabel('M, H \cdot м')
ylabel('n, об/мин')
legend('n(M)', 'n(M_k)')

%% plots of the electromechanical characteristic

% without displacement current
I_2 = U_1N ./ sqrt((r_1 + r_2_s ./ s_data).^2 + (x_s1 + x_s2_s)^2);

% with displacement current
I_2k = U_1N ./ sqrt((r_1 + k_r_data .* r_2_s ./ s_data).^2 + (x_s1 + k_x_data .* x_s2_s).^2);

%% plot I_2(s) and I_2k(s)
plot(s_data, I_2, '-', s_data, I_2k, '--', 'LineWidth', 3)
grid on
xlabel('s')
ylabel('I, A')
legend('I_2(s)', 'I_{2k}(s)')

%% plot n(I_2) and n(I_2k)
plot(I_2, n_data, '-', I_2k, n_data, '--', 'LineWidth', 3)
grid on
xlabel('I, A')
ylabel('n, об/мин')
legend('n(I_2)', 'n(I_{2k})')

%% performance criteria

% form points arrays
s_data = linspace(0, s_N + 10^(-3), 10^6);
n_data = (1 - s_data) .* n_1;
h_data = h .* abs(s_data);
k_r_data = h_data .* (sinh(2 .* h_data) + sin(2 .* h_data)) ./ (cosh(2 .* h_data) - cos(2 .* h_data));
k_x_data = 3 ./ (2 .* h_data) .* (sinh(2 .* h_data) - sin(2 .* h_data)) ./ (cosh(2 .* h_data) - cos(2 .* h_data));

% formulas without displacement current
M = m .* z_p .* U_1N.^2 .* r_2_s ./ (omega_1 .* s_data .* ((r_1 + r_2_s ./ s_data).^2 + (x_s1 + x_s2_s).^2));
I_2_s = U_1N ./ sqrt((r_1 + r_2_s ./ s_data).^2 + (x_s1 + x_s2_s)^2);
I_1 = I_2_s + U_1N ./ x_m;
P_2 = m * I_2_s.^2 * r_2_s .* (1 - s_data) ./ s_data;
P_1 = P_2 + m .* I_1.^2 * r_1 + m .* I_2_s.^2 * r_2_s;
eta = P_2 ./ P_1 .* 100;
cosphi = P_1 ./ (3 .* U_1N .* I_1N);

% formulas with displacement current
M_k = m .* z_p .* U_1N.^2 .* r_2_s .* k_r_data ./ (omega_1 .* s_data .* ((r_1 + k_r_data .* r_2_s ./ s_data).^2 + (x_s1 + x_s2_s .* k_x_data).^2));
I_2k_s = U_1N ./ sqrt((r_1 + r_2_s .* k_r_data ./ s_data).^2 + (x_s1 + x_s2_s .* k_x_data).^2);
I_1k = I_2k_s + U_1N ./ x_m;
P_2k = m * I_2k_s.^2 .* r_2_s .* (1 - s_data) ./ s_data;
P_1k = P_2k + m * I_1k.^2 * r_1 + m * I_2k_s.^2 * r_2_s;
eta_k = P_2k ./ P_1k .* 100;
cosphi_k = P_1k ./ (3 .* U_1N .* I_1N);

%% plot n(P_2) and n(P_2k)

plot(P_2, n_data, '-', P_2k, n_data, '--', 'LineWidth', 3)
grid on
xlabel('P_2, W')
ylabel('n, об/мин')
legend('n(P_2)', 'n(P_{2k})')

%% plot eta(P_2) and eta(P_2k)

plot(P_2, eta, '-', P_2k, eta_k, '--', 'LineWidth', 3)
grid on
xlabel('P_2, W')
ylabel('\eta, %')
legend('\eta(P_2)', '\eta_k(P_{2k})')

%% plot M_2(P_2) and M_2k(P_2k)

plot(P_2, M, '-', P_2k, M_k, '--', 'LineWidth', 3)
grid on
xlabel('P_2, W')
ylabel('M_2, H \cdot м')
legend('M_2(P_2)', 'M_2k(P_{2k})')

%% plot cosphi(P_2) and cosphi_k(P_2k)

plot(P_2, cosphi, '-', P_2k, cosphi_k, '--', 'LineWidth', 3)
grid on
xlabel('P_2, W')
ylabel('cos \phi')
legend('cos \phi(P_2)', 'cos \phi_k(P_{2k})')

%% plot I_1(P_2) and I_1k(P_2k)

plot(P_2, I_1, '-', P_2k, I_1k, '--', 'LineWidth', 3)
grid on
xlabel('P_2, W')
ylabel('I, A')
legend('I_1(P_2)', 'I_1k(P_{2k})')

%% plot P_1(P_2) and P_1k(P_2k)

plot(P_2, P_1, '-', P_2k, P_1k, '--', 'LineWidth', 3)
grid on
xlabel('P_2, W')
ylabel('P_1, W')
legend('P_1(P_2)', 'P_1k(P_{2k})')