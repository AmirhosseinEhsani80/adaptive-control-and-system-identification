%% Amirhossein Ehsani Simulation 4 _ Q3 _ Adaptive MPC J2
clc
clear all
close all

%% System

% Sampling time
Ts = 0.5;

% Define transfer function variable z
z = tf('z', Ts);
delay = 2;
% Define the discrete-time system
A = (z-0.1)*(z-0.65);
B = (z-0.46);
sysD = z^-delay * B / A;
[numD, denD] = tfdata(sysD, 'v');

numeratorDiscrete = numD;
denominatorDiscrete = denD;

% Removing zeros of B at the beginning
numeratorPrime = remove_first_zeros(numeratorDiscrete);

d = delay;

%% Initial conditions
A_estimated = [1, 0.1, 0.1, 0.1, 0.1];
B_estimated = [0, 10, 10]; % First element must be zero

len_desA = length(A_estimated);
len_desB = length(B_estimated);

%% Run

% Input properties
num_samples = 1000;
step_mag = 1;
t = 0:num_samples-1;
desired_output = (-1).^ceil(t/200);
desired_output(1) = -1;
input_noise_variance = 0;
input_noise = sqrt(input_noise_variance) * randn(1, num_samples);
desired_output = desired_output + input_noise;

% Noise and its degree
noise_poly = [1, 0, 0].';
len_noise = length(noise_poly);
noise_variance = 0.001; % System noise variance
system_noise = sqrt(noise_variance) * randn(1, num_samples);

% Disturbance
disturbance = zeros([1, num_samples]);

% Initial conditions for y
skip_instances = d;
output_y = [];
output_y(1:skip_instances) = 0;
control_u(1:skip_instances) = 0;

% Other initial conditions
total_parameters = len_desA - 1 + len_desB;
theta_real = [denominatorDiscrete(2:end).'; numeratorPrime.'];
theta_real_toPlot = [denominatorDiscrete(2:end), zeros(1, delay), numeratorPrime];
theta_hat_toPlot = zeros([num_samples, total_parameters]);

% RLS initial parameters
rls_solver = RLSClass(100 * eye(total_parameters), 0.01 * ones([total_parameters,1]));

instances_to_estimate_d = 30;
threshold = 0.1; % Threshold to determine d
window = 10;

numeratorPrime_estimated = B_estimated(2:end);
d_estimated = length(B_estimated) - length(numeratorPrime_estimated);

% RLS initial parameters
total_parameters = len_desA - 1 + len_desB - d_estimated; 
rls_solver_new = RLSClass(100 * eye(total_parameters), 0.1 * ones([total_parameters,1]));

% Controller parameters
beta_thres = 0.01;
u_thres = 5; 
landa = 0.6; 

for i = skip_instances:num_samples-d
    phi_t = [[-output_y(i-1:-1:max(i-(length(denominatorDiscrete) - 1), 1)), zeros(1, 1-(i-(length(denominatorDiscrete) - 1)))], ...
             [control_u(i-delay:-1:max(i-length(numeratorPrime)-delay+1, 1)), zeros(1, 1-(i-length(numeratorPrime)-delay+1))]].'; % For simulations

    phi_t_d = [[-output_y(i-1:-1:max(i-(len_desA - 1), 1)), zeros(1, 1-(i-(len_desA - 1)))], ...
               [control_u(i-d_estimated:-1:max(i-length(numeratorPrime_estimated)-d_estimated+1, 1)), zeros(1, 1-(i-length(numeratorPrime_estimated)-d_estimated+1))]].'; % For estimation
    
    noise_t = [system_noise(i:-1:max(i-(len_noise-1), 1)), zeros(1, 1-(i-(len_noise-1)))] * noise_poly;
    output_y(i) = phi_t.' * theta_real + noise_t + numeratorDiscrete * [disturbance(i:-1:max(i-(length(numeratorDiscrete)-1), 1)), zeros(1, 1-(i-(length(numeratorDiscrete) - 1)))].';

    theta_hat_new = rls_solver_new.update_RLS(output_y(i), phi_t_d);

    A_estimated = theta_hat_new(1:(len_desA - 1)).';
    numeratorPrime_estimated = theta_hat_new(len_desA:total_parameters).';
    theta_hat_toPlot(i, :) = [A_estimated, zeros(1, d_estimated), numeratorPrime_estimated];

    % Check if d is estimated correctly
    if i > instances_to_estimate_d
        if abs(mean(theta_hat_toPlot(i:-1:i-window, end - length(numeratorPrime_estimated)+1))) < threshold
            numeratorPrime_estimated = numeratorPrime_estimated(2:end);
            d_estimated = d_estimated + 1;
            total_parameters = total_parameters - 1;
            rls_solver_new = RLSClass(100 * eye(total_parameters), [A_estimated, numeratorPrime_estimated].');
        end
    end
        
    % Controller design
    n = length(A_estimated);
    n1 = length(numeratorPrime_estimated) - 1;    

    a1 = [1, A_estimated]; 
    b1 = [1];
    D = [1, zeros(1, n + d_estimated - 1)]; % *****here we have q^-1 instead of q^1****
    
    [F, G] = solve_diophantine_general(a1, b1, D);
    F = remove_first_zeros(F);

    alpha = G;
    beta = conv(F, numeratorPrime);
    beta_prime = beta(2:end); 

    if abs(beta(1)) < beta_thres
        beta(1) = beta_thres;
    end

    beta_0 = beta(1);

    control_u(i) = (beta_0 * (desired_output(i+d) - filter_signal(alpha, output_y, i) - filter_signal(beta_prime, control_u, i-1))/(beta_0^2 + landa)) / beta_0;
    
    if abs(control_u(i)) > u_thres
        control_u(i) = u_thres;
    end
end

%% Plotters
close all
figure()
plot(t(1:end-d), output_y, 'LineWidth', 1.5, 'Color', [0 0.5 1], 'DisplayName','Real Output')
hold on; 
plot(t, desired_output - input_noise, "--r", 'LineWidth', 1.5, 'DisplayName','Desired Output')
xlabel("Sample Time", 'FontSize', 12);
title("Output Response", 'FontSize', 14);
legend({'Real Output', 'Desired Output'}, 'FontSize', 12);
grid on
set(gca, 'FontSize', 12);

figure
plot(t(1:end-d), control_u, 'LineWidth', 1.5, 'Color', [0.5 0 0.5], 'DisplayName','Control Input')
xlabel("Sample Time", 'FontSize', 12);
title("Control Effort", 'FontSize', 14);
legend({'Control Input'}, 'FontSize', 12);
grid on
set(gca, 'FontSize', 12);

% Theta plot with new colors
colors = lines(length(theta_hat_toPlot(1,:)));
figure()
for j = 1:length(theta_hat_toPlot(1,:))
    title_text = sprintf("θ_%d", j);
    subplot(ceil(length(theta_hat_toPlot(1,:))/2), 2, j);
    plot(1:num_samples, theta_hat_toPlot(:, j), 'LineWidth', 1.5, 'Color', colors(j, :), 'DisplayName','Predicted')
    hold on;
    plot(1:num_samples, ones([num_samples, 1]) * theta_real_toPlot(j), "--", 'LineWidth', 1.5, 'Color', 'k', 'DisplayName','Real')
    title(title_text, 'FontSize', 14);
    legend('Location', 'best', 'FontSize', 12);
    xlabel("Sample Number", 'FontSize', 12);
    grid on
    set(gca, 'FontSize', 12);
end
