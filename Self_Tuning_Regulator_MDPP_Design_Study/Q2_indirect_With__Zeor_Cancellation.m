%% Amirhossein Ehsani - 810602159 - Q2_c
clc
clear all 
close all

%% Discretization of the Model
clc
clear;

% Define the transfer function variable
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

%% Desired System Specifications

% Desired overshoot and settling time
overshoot = 10;
settling_time = 3;

% Calculate damping ratio and natural frequency
zeta = cos(atan2(pi,-log(overshoot/100)));
wn = 4/(zeta*settling_time); 

% Desired zeros and pole
z1 = -3;
z2 = -25;
p3 = -20;

% Gain calculation
k2 = -p3/(z1*z2);

% Desired continuous-time system
G1 = tf([wn^2],[1, 2*zeta*wn, wn^2]);
G2 = zpk([z1,z2],[p3],k2);

sys_desC = G1*G2;
sys_desD = c2d(sys_desC, Ts, 'zoh'); % Discretize the desired system
[num_desD, den_desD] = tfdata(sys_desD, 'v');

BmPrime_main = num_desD; 
Am = den_desD;

% Remove leading zero if present in B
if B(1) == 0
    B = B(2:end);
end

%% InDirect with Zero Cancellation

% Input properties
num_samples = 1000;
step_mag = 1;
t = 0:num_samples-1;
uc = (-1).^-ceil(t/200);
uc(1) = -1;

% Degrees of desired A and B
deg_desA = 4;
deg_desB = 3;

% Initial conditions
A_estimated = [1,0,0,0]; % Initial A coefficients
B_estimated = [-0.3, 0.5, -0.2]; % Initial B coefficients

syms q;

% Initial B
Bplus = [B_estimated/B_estimated(1)];
Bminus = B_estimated(1);
BmPrime = BmPrime_main/B_estimated(1);

% Initial R, S, T, and Ao
R_estimated = [1,-2.5,-0.5]; 
S_estimated = [-76,13.4,-5.6];
T = [1.1,0,0];
Ao = [1];

% Real R, S, and T
R_real = [1.0000, -1.1888, 0.2396];
S_real = [5.3248, -7.3002, 2.5526];
T_real = [1.2053, -0.7324, -0.0144];

% Initial conditions for y
skip_instances = max(deg_desA, deg_desB);
total_parameters = deg_desA - 1 + deg_desB;
y = [];
y(1:skip_instances) = 0;
u(1:skip_instances) = 0;

theta_real = [A(2:end).'; B.'];

% RLS initial parameters
rls_solver = RLSClass(100 * eye(total_parameters), 0.1 * ones([total_parameters,1]));
d0 = deg_desA - deg_desB;
deg_R = length(R_estimated);
deg_S = deg_R;
deg_T = deg_R; 

R_estimated_toPlot = zeros([num_samples, deg_R]);
S_estimated_toPlot = zeros([num_samples, deg_S]);
T_toPlot = zeros([num_samples, deg_T]);
theta_hat_toPlot = zeros([num_samples, total_parameters]);

for i = skip_instances:length(uc)
    phi_t = [-y(i-1:-1:i-(deg_desA - 1)), u(i-1:-1:i-deg_desB)].';     
    y(i) = phi_t.' * theta_real;
    u(i) = T * [uc(i:-1:i-(length(T)-1))].' + S_estimated * [-y(i:-1:i-(length(S_estimated)-1))].' - R_estimated(2:end) * [u(i-1:-1:i-(length(R_estimated)-1))].';
    theta_hat_new = rls_solver.update_RLS(y(i), phi_t);

    A_estimated = theta_hat_new(1:(deg_desA - 1)).';
    B_estimated = theta_hat_new(deg_desA:total_parameters).';
    
    Bplus = [B_estimated/B_estimated(1)];
    Bminus = B_estimated(1);
    BmPrime = BmPrime_main/B_estimated(1);
    Ao = [1];
    
    [R_estimated, S_estimated, T] = Diophantine(Am, [1,A_estimated], Ao, Bminus, Bplus, BmPrime);
    
    theta_hat_toPlot(i, :) = theta_hat_new(1:total_parameters).';
    R_estimated_toPlot(i, :) = R_estimated;
    S_estimated_toPlot(i, :) = S_estimated;
    T_toPlot(i, :) = T;
end

%% Plotters

% Plot system response
figure()
plot(t, y, 'b', 'LineWidth', 1.5) % System response in blue
hold on;
plot(t, uc, '--r', 'LineWidth', 1.5) % Reference input in dashed red
xlabel("Time");
ylabel("Output");
title("System Response to Step Input");
grid on
grid minor
legend('System Output', 'Desired');
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
f2 = figure();
f2.Position = [0 0 400 900];
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
    subplot(3,2,i);
    plot(1:num_samples, theta_hat_toPlot(:,i), 'b', 'LineWidth', 1.5, 'DisplayName','Predicted')
    title(sprintf(title_text, i));
    legend('Location','best');
    xlabel("Sample Number")
    grid on
    grid minor
    set(gca, 'FontSize', 12)
end
