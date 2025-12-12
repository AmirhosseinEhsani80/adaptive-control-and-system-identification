%% Amirhossein Ehsani - 810602159 - Q2_Under Parameter
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

% Desired zero and pole
z1 = -3;
p3 = -70;

% Gain calculation
k2 = -p3/(z1);

% Desired continuous-time system
G1 = tf([wn^2],[1, 2*zeta*wn, wn^2]);
G2 = zpk([z1],[],k2);

sys_desC = G1*G2;
sys_desD = c2d(sys_desC, Ts, 'zoh'); % Discretize the desired system
[num_desD, den_desD] = tfdata(sys_desD, 'v');

BmPrime_main = num_desD; 
Am = den_desD;

% Remove leading zero if present in B
if B(1) == 0
    B = B(2:end);
end

%% Direct with Zero Cancellation

% Input properties
num_samples = 1000;
step_mag = 1;
t = 0:num_samples-1;
uc = (-1).^-ceil(t/200);
uc(1) = -1;
uc = uc;

% Degrees of desired A and B
deg_desA = 3;
deg_desB = 2;

% Initial conditions
A_estimated = [1,0,0]; % Initial A coefficients
B_estimated = [-0.3, 0.5]; % Initial B coefficients

syms q;

% Initial B
Bplus = [B_estimated/B_estimated(1)];
Bminus = B_estimated(1);
BmPrime = BmPrime_main/B_estimated(1);

% Initial R, S, T, and Ao
R_solved = [1,-1.6]; 
S_solved = [-2.5,4.3];
T = [-.125,0];
Ao = [1];

% Real R, S, and T
R_real = [1.0000, -0.9710];
S_real = [3.5332, -3.0324];
T_real = [-83.4381, 52.2274];

% Initial conditions for y
skip_instances = max(deg_desA, deg_desB) + 1;
total_parameters = deg_desA - 1 + deg_desB;
y = zeros(1, num_samples);
u = zeros(1, num_samples);

theta_real = [A(2:end).'; B.'];

% For plotting
u_toPlot = zeros(num_samples, 1);
R_solved_toPlot = zeros(num_samples, length(R_solved));
S_solved_toPlot = zeros(num_samples, length(S_solved));
T_toPlot = zeros(num_samples, length(T));
theta_hat_toPlot = zeros(num_samples, total_parameters);

% RLS initial parameters
theta_epsilon_zero = [1,0,0].';
els_solver = ELSClass(100 * eye(total_parameters + length(theta_epsilon_zero)), 0.1 * ones([total_parameters,1]), theta_epsilon_zero);
for i = skip_instances:length(uc)
    % Calculate phi_t_real and phi_t
    phi_t_real = [-y(i-1:-1:i-(deg_desA)), u(i-1:-1:i-deg_desB-1)].';     
    phi_t = [-y(i-1:-1:i-(deg_desA - 1)), u(i-1:-1:i-deg_desB)].';
    
    % Calculate y
    y(i) = phi_t_real.' * theta_real;
    
    % Calculate u
    u(i) = T * [uc(i:-1:i-(length(T)-1))].' + S_solved * [-y(i:-1:i-(length(S_solved)-1))].' - R_solved(2:end) * [u(i-1:-1:i-(length(R_solved)-1))].';
    u_toPlot(i) = u(i);
    
    % Update parameters using RLS
    theta_hat_new = els_solver.update_ELS(y(i), phi_t);
    A_estimated = theta_hat_new(1:(deg_desA - 1)).';
    B_estimated = theta_hat_new(deg_desA:total_parameters).';
    
    Bplus = [B_estimated/B_estimated(1)];
    Bminus = B_estimated(1);
    BmPrime = BmPrime_main/B_estimated(1);
    Ao = [1];

    % Solve Diophantine equation
    [R_solved,S_solved,T] = Diophantine(Am, [1,A_estimated], Ao, Bminus, Bplus, BmPrime);

    % Store values for plotting
    theta_hat_toPlot(i, :) = theta_hat_new(1:total_parameters).';
    R_solved_toPlot(i, :) = R_solved;
    S_solved_toPlot(i, :) = S_solved;
    T_toPlot(i, :) = T;
end

%% Plotters
close all

% Plot system response
figure()
plot(t, y, 'b', 'LineWidth', 1.5) % System response in blue
hold on;
plot(t, uc, '--r', 'DisplayName', 'Desired', 'LineWidth', 1.5) % Reference input in dashed red
grid on
grid minor
title("Output");
legend('Obtained','Desired');
set(gca, 'FontSize', 12)
hold off;

% Plot control input
figure()
plot(t, u, 'g', 'LineWidth', 1.5) % Control input in green
xlabel("Sample Number");
ylabel("Control Effort");
grid on
grid minor
title("Control Effort");
set(gca, 'FontSize', 12)

% R Parameters
total_plot_rows = ceil(length(R_solved)/2) + ceil(length(S_solved)/2) + ceil(length(T)/2);
f2 = figure();
f2.Position = [0 0 400 900];
for i = 2:(length(R_solved)) % first parameter of R is 1
    subplot(total_plot_rows, 2, i-1);
    title_text = "R_%d";
    plot(1:num_samples, R_solved_toPlot(:,i), 'b', 'LineWidth', 1.5, 'DisplayName','Predicted') 
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
for i = 1:(length(S_solved)) % first parameter of R is 1
    subplot(total_plot_rows, 2, i + end_of_r_plot);
    title_text = "S_%d";
    plot(1:num_samples, S_solved_toPlot(:,i), 'g', 'LineWidth', 1.5, 'DisplayName','Predicted') 
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
    subplot(2,2,i);
    plot(1:num_samples, theta_hat_toPlot(:,i), 'b', 'LineWidth', 1.5, 'DisplayName','Predicted')
    hold on
    plot(1:num_samples, ones([num_samples,1]) * theta_real(i) , 'DisplayName','Real')
    title(sprintf(title_text, i));
    legend('Location','best');
    xlabel("Sample Number")
    grid on
    grid minor
    set(gca, 'FontSize', 12)
end
