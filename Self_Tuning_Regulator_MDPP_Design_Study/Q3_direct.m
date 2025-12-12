%% Amirhossein Ehsani - 810602159 - Q2_d

% Clear workspace
clc
clear all 
close all

%% Discretization of the model
clc
clear;
s = tf('s');

% Define continuous-time system
sysC = (s-0.5)*(s+7)/((s+1)*(s+4.1)*(s-2));
BW = bandwidth(sysC);

% Define discretization parameters
disc_ratio = 1000;
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

zeta = 0.9999;

% Define desired poles and zeros
z1 = -3;
z2 = -20;
p3 = -70;

k2 = -p3/(z1*z2);

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

%% Direct without zero cancellation

% Initialize title for the main plot
main_title = "direct without ZC";

% Input properties
num_samples = 100;
step_mag = 1;
t = 0:num_samples-1;
uc = (-1).^-ceil(t/200);
uc(1) = -1;
uc = uc;

% Must be changed for every problem
deg_desA = 4;
deg_desB = 3;
A_estimated = [1,0,0,0]; % Initial conditions
B_estimated = [-0.3, 0.5, -0.2];

% Initial conditions and parameters
Bplus = [1];
Bminus = B_estimated;
Bm = conv(Bplus, Bminus);
BmPrime = (sum(Am)/sum(Bm)); 

BR_estimated = conv([1,-0.7386,-0.5508], B);
BS_estimated = conv([6.7933,-7.7807,1.9653], B);

R_estimated = [1,-25,21.2];
S_estimated = [-76,134,-55];
T_estimated = [-1,0,0];
Ao = [1,0,0];

R_real = [1,-25.2149,21.2587];
S_real = [-76.5055,134.0505,-55.9795];
T_real = [-2.8686,0,0];

d0 = deg_desA - deg_desB;

% Find degrees of R and Ao
deg_R = length(R_estimated);
deg_S = deg_R;
deg_T = deg_R; 

deg_BR = length(BR_estimated)-1;
deg_BS = deg_BR;

f_filter = conv(Am, Ao);

% Initialize variables for plotting
skip_instances = max(deg_desA, deg_BR) + d0  + 1 ;
y = [];
y(1:skip_instances) = 0;
ym(1:skip_instances) = 0;
u(1:skip_instances) = 0;
uf(1:skip_instances) = 0;
yf(1:skip_instances) = 0;
ucf(1:skip_instances) = 0;
total_parameters = deg_BR+1 + deg_BS+1 + deg_T;
R_estimated_toPlot = zeros([num_samples, deg_R]);
S_estimated_toPlot = zeros([num_samples, deg_S]);
T_toPlot = zeros([num_samples, deg_T]);
theta_hat_toPlot = zeros([num_samples, total_parameters]);

theta_real = [A(2:end).'; B.'];
theta_desired = [Am(2:end).'; Bm.'];
precision = 0.01;

theta_epsilon_zero = [1,0,0].';
els_solver = ELSClass(100 * eye(total_parameters+ length(theta_epsilon_zero)), 0.01 * ones([total_parameters,1]), theta_epsilon_zero);
rls_solver = RLSClass(100 * eye(total_parameters), 0.01 * ones([total_parameters,1]));

for i = skip_instances:length(uc)
    % Calculate phi_t
    phi_t = [-y(i-1:-1:i-(deg_desA - 1)), u(i-1:-1:i-deg_desB)].';
    
    % Calculate system output
    y(i) = phi_t.' * theta_real;
    
    % Calculate control signal
    u(i) = T_estimated * [uc(i:-1:i-(length(T_estimated)-1))].' + S_estimated * [-y(i:-1:i-(length(S_estimated)-1))].' - R_estimated(2:end) * [u(i-1:-1:i-(length(R_estimated)-1))].';
    u(i) = u(i)/R_estimated(1);

    % Filter input, output, and reference signals
    uf(i) = u(i) - f_filter(2:end) * uf(i-1:-1:i-(length(f_filter)-1)).';
    yf(i) = y(i) - f_filter(2:end) * yf(i-1:-1:i-(length(f_filter)-1)).';
    ucf(i) = uc(i) - f_filter(2:end) * ucf(i-1:-1:i-(length(f_filter)-1)).';

    % Filtered signal for calculating phi_d0
    phi_d0_filtered = [uf(i-d0:-1:i-d0-deg_BR), yf(i-d0:-1:i-d0-deg_BS), -ucf(i-d0:-1:i-d0-deg_T+1)].';
    
    % Calculate desired model output
    phi_t_m = [-y(i-1:-1:i-(deg_desA - 1)), uc(i-1:-1:i-deg_desB)].';
    ym(i) = phi_t_m.' * theta_desired;
    error = y(i) - ym(i);

    % Update parameter estimates using Extended Least Squares (ELS)
    theta_hat_new = els_solver.update_ELS(error, phi_d0_filtered);
   
    % Separate estimates for R, S, and T
    BR_estimated = theta_hat_new(1:deg_BR+1).';
    BS_estimated = theta_hat_new(deg_BR+2:deg_BR + deg_BS+2).';
    T_estimated = theta_hat_new(deg_BR + deg_BS+3:deg_BR + deg_BS+deg_T+2).';

    % Remove common roots from R and S
    [R_estimated, S_estimated] = remove_common_roots(BR_estimated,BS_estimated,precision);

    % Store estimates for plotting
    theta_hat_toPlot(i, :) = theta_hat_new(1:total_parameters).';

    % Store deg_R of R for plotting
    if length(R_estimated) >= 3
        R_estimated_toPlot(i, :) = R_estimated(1:3);
        S_estimated_toPlot(i, :) = S_estimated(1:3);
    else
        R_estimated_toPlot(i, :) = [R_estimated, zeros(1,3-length(R_estimated))];
        S_estimated_toPlot(i, :) = [S_estimated, zeros(1,3-length(S_estimated))];   
    end
    % Store deg_T of T for plotting
    if length(T_estimated) >= 3
        T_toPlot(i, :) = T_estimated(1:3);
    else
        T_toPlot(i, :) = [T_estimated, zeros(1,3-length(T_estimated))];
    end
end

R_real = R_estimated_toPlot(100,:);
S_real = S_estimated_toPlot(100,:);
T_real = T_toPlot(100,:);
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
total_plot_rows = ceil(length(R_estimated_toPlot(1,:))/2) + ceil(length(S_estimated_toPlot(1,:))/2) + ceil(length(T_toPlot(1,:))/2) - 2;
f2 = figure();
f2.Position = [0 0 400 900];
for i = 2:(length(R_estimated_toPlot(1,:))) % first parameter of R is 1
    subplot(total_plot_rows,2,i-1);
    title_text = "R_%d";
    plot(1:num_samples, R_estimated_toPlot(:,i), 'b', 'LineWidth', 1.5, 'DisplayName','Predicted') 
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
for i = 1:(length(S_estimated_toPlot(1,:))) % first parameter of R is 1
    subplot(total_plot_rows,2,i + end_of_r_plot);
    title_text = "S_%d";
    plot(1:num_samples, S_estimated_toPlot(:,i), 'g', 'LineWidth', 1.5, 'DisplayName','Predicted') 
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
for i = 1:6
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
