%% Amirhossein Ehsani - 810602159 - Q2_c

% Clear workspace
clc
clear all 
close all

%% Discretization of the model

% transfer function variable
s = tf('s');

% Define continuous-time system
sysC = (s-0.5)*(s+7)/((s+1)*(s+4.1)*(s-2));
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

% Remove leading zeros from numerator if present
if B(1) == 0
    B = B(2:end);
end

%% Desired system

% Define desired system specifications
overshoot = 10;
settling_time = 3;

% Calculate desired damping ratio and natural frequency
zeta = cos(atan2(pi,-log(overshoot/100)));
wn = 4/(zeta*settling_time); 

% Define desired poles and zeros
z1 = -3;
z2 = -20;
p3 = -70;

k2 = -p3;

% Create transfer functions for desired system components
G1 = tf([wn^2],[1, 2*zeta*wn, wn^2]);
G2 = zpk([z1,z2],[p3],k2);

% Create continuous-time and discrete-time desired systems
sys_desC = G1*G2;
sys_desD = c2d(sys_desC, Ts, 'zoh');
[num_desD, den_desD] = tfdata(sys_desD, 'v');

% Define numerator and denominator of the discrete-time desired system
Bm_prim2 = num_desD; 
Am = den_desD;

%% Indirect without zero cancellation

% Input properties
num_samples = 1000;
step_mag = 1;
t = 0:num_samples-1;
uc = (-1).^-ceil(t/200);
uc(1) = -1;

% Must be changed for every problem
deg_desA = 4;
deg_desB = 3;
A_estimated = [1,0,0,0]; % Initial conditions
B_estimated = [-0.3, 0.5, -0.2];

% Initial B
Bplus = [1];
Bminus = B_estimated;
Bm = conv(Bplus, Bminus);
BmPrime = (sum(Am)/sum(Bm)); 

% Initial R S T and Ao
R_solved = [1,-25,21.2];
S_solved = [-76,134,-55];
T = [-1,0,0];
Ao = [1,0,0];


% Initial conditions for y
skip_instances = max(deg_desA, deg_desB);
total_parameters = deg_desA - 1 + deg_desB;
y = [];
y(1:skip_instances) = 0;
u(1:skip_instances) = 0;

% True system parameters
theta_real = [A(2:end).'; B.'];

% Define arrays for plotting
R_solved_toPlot = zeros([num_samples, length(R_solved)]);
S_solved_toPlot = zeros([num_samples, length(S_solved)]);
T_toPlot = zeros([num_samples, length(T)]);
theta_hat_toPlot = zeros([num_samples, total_parameters]);

% RLS initial parameters
theta_epsilon_zero = [1,0,0].';
els_solver = ELSClass(100 * eye(total_parameters + length(theta_epsilon_zero)), 0.1 * ones([total_parameters,1]), theta_epsilon_zero);

for i = skip_instances:length(uc)
    % Calculate phi_t
    phi_t = [-y(i-1:-1:i-(deg_desA - 1)), u(i-1:-1:i-deg_desB)].';   

    % Calculate y(i)
    y(i) = phi_t.' * theta_real ;

    % Calculate u(i)
    u(i) = T * [uc(i:-1:i-(length(T)-1))].' + S_solved * [-y(i:-1:i-(length(S_solved)-1))].' - R_solved(2:end) * [u(i-1:-1:i-(length(R_solved)-1))].';

    % Update parameter estimates using ELS
    theta_hat_new = els_solver.update_ELS(y(i), phi_t);

    % Update estimated A and B
    A_estimated = theta_hat_new(1:(deg_desA - 1)).';
    B_estimated = theta_hat_new(deg_desA:total_parameters).';
    
    % Update Bplus, Bminus, Bm, and BmPrime
    Bplus = [1];
    Bminus = B_estimated;
    Bm = conv(Bplus, Bminus);
    BmPrime = (sum(Am)/sum(Bm)); 
    
    % Solve Diophantine equation
    [R_solved,S_solved,T] = Diophantine(Am, [1,A_estimated], Ao, Bminus, Bplus, BmPrime);

    % Store data for plotting
    theta_hat_toPlot(i, :) = theta_hat_new(1:total_parameters).';
    R_solved_toPlot(i, :) = R_solved;
    S_solved_toPlot(i, :) = S_solved;
    T_toPlot(i, :) = T;
end
% R S and T which were calculated in the last part
R_real = R_solved;
S_real = S_solved;
T_real = T;

%% Plotters

% Plot system response
figure()
plot(t, y, 'b', 'LineWidth', 1.5) % System response in blue
hold on;
plot(t, uc, '--r', 'DisplayName', 'Desired', 'LineWidth', 1.5) % Reference input in dashed red
xlabel("Time (samples)");
ylabel("Output");
title("System Response to Step Input");
grid on

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

set(gca, 'FontSize', 12)

% R Parameters
total_plot_rows = ceil(length(R_solved_toPlot(1,:))/2) + ceil(length(S_solved_toPlot(1,:))/2) + ceil(length(T_toPlot(1,:))/2) - 2;
f2 = figure();
f2.Position = [0 0 400 900];
for i = 2:(length(R_solved_toPlot(1,:))) % first parameter of R is 1
    subplot(total_plot_rows,2,i-1);
    title_text = "R_%d";
    plot(1:num_samples, R_solved_toPlot(:,i), 'b', 'LineWidth', 1.5, 'DisplayName','Predicted') 
    hold on;
    plot(1:num_samples, ones([num_samples,1]) * R_real(i), 'r--', 'LineWidth', 1.5, 'DisplayName','Real') 
    title(sprintf(title_text, i-1));
    legend('Obtained','Desired');
    xlabel("Sample Number")
    grid on
    
    set(gca, 'FontSize', 12)
end
end_of_r_plot = i-1;

% S Parameters
for i = 1:(length(S_solved_toPlot(1,:))) % first parameter of R is 1
    subplot(total_plot_rows,2,i + end_of_r_plot);
    title_text = "S_%d";
    plot(1:num_samples, S_solved_toPlot(:,i), 'g', 'LineWidth', 1.5, 'DisplayName','Predicted') 
    hold on;
    plot(1:num_samples, ones([num_samples,1]) * S_real(i), 'm--', 'LineWidth', 1.5, 'DisplayName','Real') 
    title(sprintf(title_text, i));
    legend('Obtained','Desired');
    xlabel("Sample Number")
    grid on
    
    set(gca, 'FontSize', 12)
end

end_of_S_plot = i + end_of_r_plot;

% T Parameters
for i = 1:(length(T_toPlot(1,:))) % first parameter of R is 1
    subplot(total_plot_rows,2,i+end_of_S_plot);
    title_text = "T_%d";
    plot(1:num_samples, T_toPlot(:,i), 'k', 'LineWidth', 1.5, 'DisplayName','Predicted') 
    hold on;
    plot(1:num_samples, ones([num_samples,1]) * T_real(i), 'c--', 'LineWidth', 1.5, 'DisplayName','Real') 
    title(sprintf(title_text, i));
    legend('Obtained','Desired');
    xlabel("Sample Number")
    grid on
    
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
    set(gca, 'FontSize', 12)
end
