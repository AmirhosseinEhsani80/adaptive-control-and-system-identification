%% Amirhossein Ehsani - 810602159 - Q2_Disturbance with Integrator

% Clear workspace
clc
clear all 
close all

%% Discretization of the model

% Define transfer function variable
s = tf('s');

% Define continuous-time system
sysC  = (s+0.5)*(s+7)/((s+1)*(s+4.1)*(s-2));
BW = bandwidth(sysC);  

% Define discretization parameters
disc_ratio = 10;
fs = BW * disc_ratio/(2*pi);
Ts = 1/fs;

% Convert continuous-time system to discrete-time
sysD = c2d(sysC, Ts, 'zoh');
[numD, denD] = tfdata(sysD, 'v');

% Define numerator and denominator of the discrete-time system
B = numD;
A = denD;

%% Desired system

% Define desired system specifications
overshoot = 10;
settling_time = 3;

% Calculate desired damping ratio and natural frequency
zeta = cos(atan2(pi,-log(overshoot/100)));
wn = 4/(zeta*settling_time); 

% Define desired poles and zeros
z1 = -3;
z2 = -1;
p3 = -50;

% Calculate gain for desired system
k2 = -p3/(z1*z2);

% Create transfer functions for desired system components
G1 = tf([wn^2],[1, 2*zeta*wn, wn^2]);
G2 = zpk([z1,z2],[p3],k2);

% Create continuous-time and discrete-time desired systems
sys_desC = G1*G2;
sys_desD = c2d(sys_desC, Ts, 'zoh');
[num_desD, den_desD] = tfdata(sys_desD, 'v');

% Define numerator and denominator of the discrete-time desired system
BmPrime_main = num_desD; 
Am = den_desD;

% Remove leading zeros from numerator if present
if B(1) == 0
    B = B(2:end);
end

%% Direct with zero cancellation 

% Input properties
num_samples = 1000;
step_mag = 1;
t = 0:num_samples-1;
uc = (-1).^-ceil(t/200);
uc(1) = -1;
uc = uc;

% Must be changed for every problem
deg_desA = 4;
deg_desB = 3;
A_estimated = [1,0,0,0]; % Initial conditions
B_estimated = [-0.3, 0.5, -0.1];

% Define symbolic variable
syms q;

% Initial B
Bplus = [B_estimated/B_estimated(1)];
Bminus = B_estimated(1);
BmPrime = BmPrime_main/B_estimated(1);

% Initial R S T and Ao
R_solved = [0 0 0]; 
S_solved = [0 0 0];
T = [-0,0,0];
Ao = [1];

% R S and T which were calculated in the last part
R_real = [1.0000   -1.2062    0.2378];
S_real = [5.5540   -7.5999    2.5986];
T_real = [9.3684  -14.2210    5.3327];

% Initial conditions for y
skip_instances = max(deg_desA, deg_desB);
total_parameters = deg_desA - 1 + deg_desB;
y = [];
y(1:skip_instances) = 0;
u(1:skip_instances) = 0;

% True system parameters
theta_real = [A(2:end).'; B.'];

% Disturbance
v = [zeros([1,ceil(num_samples/2)]), 1*ones([1,ceil(num_samples/2)])];

% Integrator
X = [1, 0];

% For plotting
u_toPlot = zeros([num_samples, length(uc)]);
R_solved_toPlot = zeros([num_samples, length(R_solved)]);
S_solved_toPlot = zeros([num_samples, length(S_solved)]);
T_toPlot = zeros([num_samples, length(T)]);
theta_hat_toPlot = zeros([num_samples, total_parameters]);

% RLS initial parameters
theta_epsilon_zero = [1,0,0].';
els_solver = ELSClass(100 * eye(total_parameters + length(theta_epsilon_zero)), 0.1 * ones([total_parameters,1]), theta_epsilon_zero);

for i = skip_instances:length(uc)
    % Calculate Y
    Y = -((1+X(end)) * polyval(R_solved,1))/polyval(B_estimated, 1);
    
    % Estimate R and S with integrator
    R_estimated_withI = conv(X, R_solved) + [0,conv(Y, B_estimated)];
    S_estimated_withI = conv(X, S_solved) - conv(Y, A_estimated);

    % Define phi_t_real and phi_t
    phi_t_real = [-y(i-1:-1:i-(deg_desA - 1)), u(i-1:-1:i-deg_desB)].';  
    phi_t = [-y(i-1:-1:i-(deg_desA - 1)), u_toPlot(i-1:-1:i-deg_desB)].'; 
    
    % Calculate y(i)
    y(i) = phi_t_real.' * theta_real;
    
    % Calculate u(i)
    u(i) = T * [uc(i-1:-1:i-(length(T)))].' + S_estimated_withI * [-y(i:-1:i-(length(S_estimated_withI)-1))].' ...
        - R_estimated_withI(2:end) * [u_toPlot(i-1:-1:i-(length(R_estimated_withI)-1))].' + v(i); 
    
    % Pure input
    u_toPlot(i) = u(i) - v(i); 
    
    % Update parameter estimates using ELS
    theta_hat_new = els_solver.update_ELS(y(i), phi_t);

    % Update estimated A and B
    A_estimated = theta_hat_new(1:(deg_desA - 1)).';
    A_estimated = [1, A_estimated];
    B_estimated = theta_hat_new(deg_desA:total_parameters).';
    
    % Update Bplus, Bminus, BmPrime, and Ao
    Bplus = [B_estimated/B_estimated(1)];
    Bminus = B_estimated(1);
    BmPrime = BmPrime_main/B_estimated(1);
    Ao = [1];

    % Solve Diophantine equation
    [R_solved,S_solved,T] = Diophantine(Am, A_estimated, Ao, Bminus, Bplus, BmPrime);

    % Store data for plotting
    theta_hat_toPlot(i, :) = theta_hat_new(1:total_parameters).';
    R_solved_toPlot(i, :) = R_solved;
    S_solved_toPlot(i, :) = S_solved;
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
total_plot_rows = ceil(length(R_solved_toPlot(1,:))/2) + ceil(length(S_solved_toPlot(1,:))/2) + ceil(length(T_toPlot(1,:))/2) - 2;
f2 = figure();
f2.Position = [0 0 400 900];
for i = 2:(length(R_solved_toPlot(1,:))) % first parameter of R is 1
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
for i = 1:(length(S_solved_toPlot(1,:))) % first parameter of R is 1
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
    subplot(3,2,i);
    plot(1:num_samples, theta_hat_toPlot(:,i), 'b', 'LineWidth', 1.5, 'DisplayName','Predicted')
    hold on
    plot(1:num_samples, ones([num_samples,1]) * theta_real(i) , '--r', 'LineWidth', 1.5,'DisplayName','Real')
    title(sprintf(title_text, i));
    legend('Location','best');
    xlabel("Sample Number")
    grid on
    grid minor
    set(gca, 'FontSize', 12)
end
