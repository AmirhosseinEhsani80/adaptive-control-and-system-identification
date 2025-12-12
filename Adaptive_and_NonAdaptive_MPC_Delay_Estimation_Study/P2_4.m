%% Amirhossein Ehsani Simulation 4 _ Q2 _ MPC Constant Future
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
[numDiscrete, denDiscrete] = tfdata(sysD, 'v');

numeratorDiscrete = numDiscrete;
denominatorDiscrete = denDiscrete;

% Removing zeros of B at the beginning
numeratorPrime = remove_first_zeros(numeratorDiscrete);

%% One step ahead
n = length(denominatorDiscrete) - 1;
n1 = length(numeratorPrime) - 1;

%% Finding F and G
d = delay * 2;
a1 = [denominatorDiscrete(1:end-delay), zeros(1, d-1)]; 
na1 = length(a1) - 1;
b1 = [1];
D = [1, zeros(1, na1 + d - 1)]; 

[F, G] = solve_diophantine_general(a1, b1, D);
F = remove_first_zeros(F);

%% Finding R and R_bar

LHS = conv(numeratorPrime, F);
LHS = remove_first_zeros(LHS);
D1 = conv(LHS, [1, zeros(1, d - delay)]);

R = D1(1:d - delay + 1);
R_bar = D1(d - delay + 2:end);

%% Simulation

question1 = input('Do we have noise?\n 1. YES\n 2. NO\n Input:  ');
question2 = input('Do we have disturbance?\n 1. YES\n 2. NO\nInput:  ');

% Input properties
num_samples = 1000;
step_mag = 1;
t = 0:num_samples-1;
desired_output = (-1).^ceil(t/200);
desired_output(1) = -1;
input_noise_variance = 0.001;
input_noise = sqrt(input_noise_variance) * randn(1, num_samples);

if question1 == 1
    desired_output = desired_output + input_noise;
else
    input_noise = 0;
end

% Noise and its degree
noise_poly = [1].';
len_noise = length(noise_poly);
noise_variance = 0.001; % System noise variance
noise_variance = 0.00; % System noise variance
system_noise = sqrt(noise_variance) * randn(1, num_samples);

% Disturbance
if question2 == 1
    disturbance = [zeros([1, ceil(num_samples/2)]), 8*ones([1, ceil(num_samples/2)])];
else
    disturbance = zeros([1, num_samples]);
end

len_desA = 5;
len_desB = 2;

% Initial conditions for y
skip_instances = d;
total_parameters = len_desA - 1 + len_desB;
output_y = [];

output_y(1:skip_instances) = 0;
control_u(1:skip_instances) = 0;
u_bar(1:skip_instances) = 0;

% Other initial conditions:
theta_real = [denominatorDiscrete(2:end).'; numeratorPrime.'];

% Controller parameters
for i = skip_instances:num_samples-d     
    phi_t = [[-output_y(i-1:-1:max(i-(length(denominatorDiscrete) - 1), 1)), zeros(1, 1-(i-(length(denominatorDiscrete) - 1)))], ...
             [control_u(i-delay:-1:max(i-length(numeratorPrime)-delay+1, 1)), zeros(1, 1-(i-length(numeratorPrime)-delay+1))]].';  % Added zeros for initial conditions
    noise_t = [system_noise(i:-1:max(i-(len_noise-1), 1)), zeros(1, 1-(i-(len_noise-1)))] * noise_poly;
    output_y(i) = phi_t.' * theta_real + noise_t + numeratorDiscrete * [disturbance(i:-1:max(i-(length(numeratorDiscrete)-1), 1)), zeros(1, 1-(i-(length(numeratorDiscrete) - 1)))].';
    control_u(i) = (desired_output(i+d) - filter_signal(G, output_y, i) - filter_signal(R_bar(1:end), control_u, i-1)) / sum(R);
end

%% Plotters

% Output response plot
figure()
plot(t(1:end-d), output_y, 'LineWidth', 1.5, 'Color', [0 0.5 1], 'DisplayName','Real Output')
hold on; 
plot(t, desired_output - input_noise, "--r", 'LineWidth', 1.5, 'DisplayName','Desired Output')
xlabel("Sample Time", 'FontSize', 12);
title("Output Response", 'FontSize', 14);
legend({'Real Output', 'Desired Output'}, 'FontSize', 12);
grid on
set(gca, 'FontSize', 12);

% Control effort plot
figure
plot(t(1:end-d), control_u, 'LineWidth', 1.5, 'Color', [0.5 0 0.5], 'DisplayName','Control Input')
xlabel("Sample Time", 'FontSize', 12);
title("Control Effort", 'FontSize', 14);
legend({'Control Input'}, 'FontSize', 12);
grid on
set(gca, 'FontSize', 12);
