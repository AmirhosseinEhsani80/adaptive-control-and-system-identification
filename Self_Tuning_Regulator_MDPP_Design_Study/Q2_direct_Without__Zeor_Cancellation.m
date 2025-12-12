%% Amirhossein Ehsani - 810602159 - Q2_b
clc
clear all 
close all

%% Discretization of the Model

s = tf('s');

% Continuous-time system transfer function
sysC = (s+0.5)*(s+7)/((s+1)*(s+4.1)*(s-2));

% Calculate the bandwidth of the continuous-time system
BW = bandwidth(sysC); 

% Discretization ratio and sampling frequency
disc_ratio = 10;
fs = BW * disc_ratio/(2*pi);
Ts = 1/fs; % Sampling time

% Discretize the continuous-time system using zero-order hold
sysD = c2d(sysC, Ts, 'zoh');
[numD, denD] = tfdata(sysD, 'v'); % Get numerator and denominator

B = numD;
A = denD;

% Remove leading zero if present in B
if B(1) == 0
    B = B(2:end);
end

%% Desired System Specifications

% Desired overshoot and settling time
overshoot = 10;
settling_time = 3;

% Calculate damping ratio and natural frequency
zeta = cos(atan2(pi,-log(overshoot/100)));
wn = 4/(zeta*settling_time); 
zeta = 0.999;

% Desired zeros and pole
z1 = -1;
z2 = -20;
p3 = -70;

% Gain calculation
k2 = -p3/(z1*z2);

% Desired continuous-time system
G1 = tf([wn^2],[1, 2*zeta*wn, wn^2]);
G2 = zpk([z1,z2],[p3],k2);

sys_desC = G1*G2;
sys_desD = c2d(sys_desC, Ts, 'zoh'); % Discretize the desired system
[num_desD, den_desD] = tfdata(sys_desD, 'v');

Bm_prim2 = num_desD; 
Am = den_desD;

%% Direct Without Zero Cancellation

syms q;

% Input properties
num_samples = 1000;

%num_samples = 10000; %for long time simulation
step_mag = 1;
t = 0:num_samples-1;
uc = (-1).^-ceil(t/200);
uc(1) = -1;
uc = uc;

% Degrees of desired A and B
deg_desA = 4;
deg_desB = 3;

% Initial conditions
A_estimated = [1,0,0,0]; % Initial A coefficients
B_estimated = [-0.3, 0.5, -0.2]; % Initial B coefficients

% Initial B
Bplus = [1];
Bminus = B_estimated;
Bm = conv(Bplus, Bminus);
BmPrime = (sum(Am)/sum(Bm)); 
BR_estimated = conv([1.0000   0, 0.1], B);
BS_estimated = conv([1.0000   0.2, 0.1], B);

% Initial R, S, T, and Ao
R_estimated = [1, -2.6253, 0.6498];
S_estimated = [18.2185,-23.7952,6.7085];
T_estimated = [3.3711,0,0];
Ao = [1,0,0];

% Real R, S, and T
R_real = [-0.2996, 0.8489 ,-0.9232];
S_real = [-0.2402 ,0.3583 ,-0.1272];
T_real = [-0.1473, 0.3021 ,-0.1550];

d0 = deg_desA - deg_desB;
deg_R = length(R_estimated);
deg_S = deg_R;
deg_T = deg_R; 

% Degree of BR and BS
deg_BR = length(BR_estimated)-1;
deg_BS = deg_BR;

% Filtering
f_filter = conv(Am, Ao);

% Initial conditions for y, ym, u, uf, yf, and ucf
skip_instances = max(deg_desA, deg_BR) + d0 + 1;
y = zeros(1, num_samples);
ym = zeros(1, num_samples);
u = zeros(1, num_samples);
uf = zeros(1, num_samples);
yf = zeros(1, num_samples);
ucf = zeros(1, num_samples);

total_parameters = deg_BR + 1 + deg_BS + 1 + deg_T;
R_estimated_toPlot = zeros([num_samples, deg_R]);
S_estimated_toPlot = zeros([num_samples, deg_S]);
T_toPlot = zeros([num_samples, deg_T]);
theta_hat_toPlot = zeros([num_samples, total_parameters]);

theta_real = [Am(2:end).'; Bm.'];
theta_desired = [Am(2:end).'; Bm.'];

precision = 0.001;

theta_epsilon_zero = [1,0,0].';
els_solver = ELSClass(100 * eye(total_parameters + length(theta_epsilon_zero)), 0.01 * ones([total_parameters,1]), theta_epsilon_zero);

for i = skip_instances:length(uc)
    % Calculate phi_t
    phi_t = [-y(i-1:-1:i-(deg_desA - 1)), u(i-1:-1:i-deg_desB)].'; 
    
    % Calculate y and u
    y(i) = phi_t.' * theta_real;
    u(i) = T_estimated * [uc(i:-1:i-(length(T_estimated)-1))].' + S_estimated * [-y(i:-1:i-(length(S_estimated)-1))].' - R_estimated(2:end) * [u(i-1:-1:i-(length(R_estimated)-1))].';
    u(i) = u(i)./R_estimated(1);

    % Filtering
    uf(i) = u(i) - f_filter(2:end) * uf(i-1:-1:i-(length(f_filter)-1)).';
    yf(i) = y(i) - f_filter(2:end) * yf(i-1:-1:i-(length(f_filter)-1)).';
    ucf(i) = uc(i) - f_filter(2:end) * ucf(i-1:-1:i-(length(f_filter)-1)).';

    phi_d0_filtered = [uf(i-d0:-1:i-d0-deg_BR), yf(i-d0:-1:i-d0-deg_BS), -ucf(i-d0:-1:i-d0-deg_T+1)].';
    phi_t_m = [-y(i-1:-1:i-(deg_desA - 1)), uc(i-1:-1:i-deg_desB)].';
    ym(i) = phi_t_m.' * theta_desired;

    error = y(i) - ym(i);
    theta_hat_new = els_solver.update_ELS(error, phi_d0_filtered);
   
    BR_estimated = theta_hat_new(1:deg_BR+1).';
    BS_estimated = theta_hat_new(deg_BR+2:deg_BR + deg_BS+2).';
    T_estimated = theta_hat_new(deg_BR + deg_BS+3:deg_BR + deg_BS+deg_T+2).';

    [R_estimated, S_estimated] = remove_common_roots(BR_estimated, BS_estimated, precision);

    theta_hat_toPlot(i, :) = theta_hat_new(1:total_parameters).';

    % Plotting only first 3 coefficients of R and S
    if length(R_estimated) >= 3
        R_estimated_toPlot(i, :) = R_estimated(1:3);
        S_estimated_toPlot(i, :) = S_estimated(1:3);
    else
        R_estimated_toPlot(i, :) = [R_estimated, zeros(1,3-length(R_estimated))];
        S_estimated_toPlot(i, :) = [S_estimated, zeros(1,3-length(S_estimated))];
    end

    if length(T_estimated) >= 3
        T_toPlot(i, :) = T_estimated(1:3);
    else
        T_toPlot(i, :) = [T_estimated, zeros(1,3-length(T_estimated))];
    end
end

%% Plotters

% Plot system response
figure()
plot(t, y, 'b', 'LineWidth', 1.5) % System response in blue
hold on;
plot(t, uc, '--r', 'DisplayName','Desired', 'LineWidth', 1.5) % Reference input in dashed red
xlabel("Time (samples)");
ylabel("Output");
title("System Response to Step Input");
grid on
grid minor
legend('System Output','Desired');
set(gca, 'FontSize', 12)
hold off;

% Plot control input
figure()
plot(t, u, 'g', 'LineWidth', 1.5) % Control input in green
xlabel("Sample Number");
ylabel("Control Effort");
title("Control Effort");
grid on
grid minor
set(gca, 'FontSize', 12)

% R Parameters
total_plot_rows = ceil(length(R_estimated_toPlot(1,:))/2) + ceil(length(S_estimated_toPlot(1,:))/2) + ceil(length(T_toPlot(1,:))/2) - 2;
figure()
for i = 2:(length(R_estimated_toPlot(1,:))) % first parameter of R is 1
    subplot(total_plot_rows, 2, i-1);
    title_text = "R_%d";
    plot(1:num_samples, R_estimated_toPlot(:,i), 'b', 'LineWidth', 1.5, 'DisplayName','Predicted') 
    hold on;
    plot(1:num_samples, ones([num_samples,1]) * R_real(i), 'r--', 'LineWidth', 1.5, 'DisplayName','Real') 
    title(sprintf(title_text, i-1));
    legend('Location','best');
    xlabel("Sample Number")
    grid on
    grid minor
    set(gca, 'FontSize', 12)
end
end_of_r_plot = i-1;

% S Parameters
for i = 1:(length(S_estimated_toPlot(1,:))) % first parameter of R is 1
    subplot(total_plot_rows, 2, i + end_of_r_plot);
    title_text = "S_%d";
    plot(1:num_samples, S_estimated_toPlot(:,i), 'g', 'LineWidth', 1.5, 'DisplayName','Predicted') 
    hold on;
    plot(1:num_samples, ones([num_samples,1]) * S_real(i), 'm--', 'LineWidth', 1.5, 'DisplayName','Real') 
    title(sprintf(title_text, i));
    legend('Location','best');
    xlabel("Sample Number")
    grid on
    grid minor
    set(gca, 'FontSize', 12)
end

end_of_S_plot = i + end_of_r_plot;

% T Parameters
for i = 1:(length(T_toPlot(1,:))) % first parameter of R is 1
    subplot(total_plot_rows, 2, i+end_of_S_plot);
    title_text = "T_%d";
    plot(1:num_samples, T_toPlot(:,i), 'k', 'LineWidth', 1.5, 'DisplayName','Predicted') 
    hold on;
    plot(1:num_samples, ones([num_samples,1]) * T_real(i), 'c--', 'LineWidth', 1.5, 'DisplayName','Real') 
    title(sprintf(title_text, i));
    legend('Location','best');
    xlabel("Sample Number")
    grid on
    grid minor
    set(gca, 'FontSize', 12)
end

% Theta Parameters
f5 = figure();
for i = 1:length(theta_hat_toPlot(1,:))
    title_text = "θ_%d";
    subplot(7,2,i);
    plot(1:num_samples, theta_hat_toPlot(:,i), 'b', 'LineWidth', 1.5, 'DisplayName','Predicted')
    title(sprintf(title_text, i));
    legend('Location','best');
    xlabel("Sample Number")
    grid on
    grid minor
    set(gca, 'FontSize', 12)
end
