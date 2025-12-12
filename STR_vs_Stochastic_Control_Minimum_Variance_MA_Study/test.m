%% Amirhossein Ehsani_810602159_ Sim 3_ Q1_DirSTR _Without Noise
clc;
clear all;
close all;

%% Discretization of the Model
z = tf('z', 0.1); % Discrete-time transfer function with sampling time 0.1 seconds
sysDiscrete = (5.401*10^(-7)*z-3.7807*10^(-7))/((2.141*10^(-6)*z^2- 2.0554*10^(-6)*z+ 2.7405*10^(-7)));
[numDiscrete, denDiscrete] = tfdata(sysDiscrete, 'v'); % Extract numerator and denominator
samplingTime = 0.1; % Sampling time

% Normalize B and A polynomials
B_discrete = numDiscrete/(5.401*10^(-7));
A_discrete = denDiscrete/(2.141*10^(-6));

if B_discrete(1) == 0
    B_discrete = B_discrete(2:end);
end

%% Desired System Parameters
desiredOvershoot = 0.15;
desiredRiseTime = 0.1;

dampingRatio = cos(atan2(pi, -log(desiredOvershoot))); % Damping ratio
naturalFrequency = 1.8/(desiredRiseTime); % Natural frequency

desiredZero = -20; % Zero of the desired system
gainFactor = -1/desiredZero;

G1_desired = tf([naturalFrequency^2], [1, 2*dampingRatio*naturalFrequency, naturalFrequency^2]);
G2_desired = zpk([desiredZero], [], gainFactor);
sysDesiredC = G1_desired * G2_desired; % Desired continuous system
sysDesiredD = c2d(sysDesiredC, samplingTime, 'zoh'); % Discretize desired system
[numDesiredD, denDesiredD] = tfdata(sysDesiredD, 'v');

B_prime_main = numDesiredD; 
A_m = denDesiredD;

if B_discrete(1) == 0
    B_discrete = B_discrete(2:end);
end
if B_prime_main(1) == 0
    B_prime_main = B_prime_main(2:end);
end

%% Direct STR with Zero Cancellation

% Input properties
numberOfSamples = 1000;
stepMagnitude = 1;
timeVector = 0:numberOfSamples-1;
controlSignal = (-1).^-ceil(timeVector/200);
controlSignal(1) = -1;

% Initial conditions for estimated polynomials
degA = 3;
degB = 2;
A_estimate = [1, -1, 0.15]; % Initial A polynomial
B_estimate = [1, -0.69]; % Initial B polynomial

% Desired B
B_plus = [B_prime_main/B_prime_main(1)];
B_minus = B_prime_main(1);
B_m = conv(B_plus, B_minus);
B_prime = B_prime_main/B_prime_main(1);

% Initial R, S, T, and A0 polynomials
R_estimate = [1, -0.1]; 
S_estimate = [1, 0.1];
T_estimate = [1.1, 0];
A_0 = [1];

% R, S, and T calculated in the last part
R_real = [0.9874   -0.7105];
S_real = [0.9904    0.0551];
T_real = [1.1472    0.0533];

d0 = degA - degB;
[degR, degAo] = find_degrees(A_m, A_discrete, A_0, B_minus, B_plus, B_prime);
degS = degR;
degT = degR; 

f_filter = conv(A_m, A_0);

assert(f_filter(1) == 1, "AmAo is not monic");

% Initial conditions for y
skipInstances = max(degA, degB);

output = zeros(1, numberOfSamples);
modelOutput = zeros(1, numberOfSamples);
input = zeros(1, numberOfSamples);
uf = zeros(1, numberOfSamples);
yf = zeros(1, numberOfSamples);
ucf = zeros(1, numberOfSamples);

totalParams = degR + 1 + degS + 1 + degT + 1;

% For plotting
R_estimate_toPlot = zeros(numberOfSamples, length(R_estimate));
S_estimate_toPlot = zeros(numberOfSamples, length(S_estimate));
T_toPlot = zeros(numberOfSamples, length(T_estimate));
theta_hat_toPlot = zeros(numberOfSamples, totalParams);

theta_real = [A_discrete(2:end).'; B_discrete.'];
theta_desired = [A_m(2:end).'; B_m.'];

% RLS initial parameters
theta_init = [0.5, 0.5].';
els_solver = ELSClass(100 * eye(totalParams + length(theta_init)), 0.01 * ones(totalParams, 1), theta_init);

% Main loop for parameter estimation and control
for i = skipInstances:numberOfSamples
    phi = [-output(i-1:-1:i-(degA - 1)), input(i-1:-1:i-degB)].'; % Assumed degB = degB
    output(i) = phi.' * theta_real;
    input(i) = T_estimate * [controlSignal(i:-1:i-(length(T_estimate)-1))].' + S_estimate * [-output(i:-1:i-(length(S_estimate)-1))].' - R_estimate(2:end) * [input(i-1:-1:i-(length(R_estimate)-1))].';
    input(i) = input(i) / R_estimate(1);

    uf(i) = input(i) - f_filter(2:end) * uf(i-1:-1:i-(length(f_filter)-1)).';
    yf(i) = output(i) - f_filter(2:end) * yf(i-1:-1:i-(length(f_filter)-1)).';
    ucf(i) = controlSignal(i) - f_filter(2:end) * ucf(i-1:-1:i-(length(f_filter)-1)).';

    phi_d0_filtered = [uf(i-d0:-1:i-d0-degR), yf(i-d0:-1:i-d0-degS), -ucf(i-d0:-1:i-d0-degT)].';
    
    phi_t_m = [-output(i-1:-1:i-(degA - 1)), controlSignal(i-1:-1:i-degB)].';
    modelOutput(i) = phi_t_m.' * theta_desired;
    error = output(i) - modelOutput(i);

    theta_hat_new = els_solver.update_ELS(error, phi_d0_filtered);
   
    R_estimate = theta_hat_new(1:degR+1).';
    S_estimate = theta_hat_new(degR+2:degR + degS+2).';
    T_estimate = theta_hat_new(degR + degS+3:degR + degS+degT+3).';

    theta_hat_toPlot(i, :) = theta_hat_new(1:totalParams).';
    R_estimate_toPlot(i, :) = R_estimate;
    S_estimate_toPlot(i, :) = S_estimate;
    T_toPlot(i, :) = T_estimate;
end

%% Metrics Calculation
varianceY = var(output);
varianceU = var(input);
varianceAbsY = var(abs(output));
varianceAbsU = var(abs(input));
meanY = mean(output);
meanU = mean(input);

accumulatedLossY = cumsum(output.^2);
accumulatedLossU = cumsum(input.^2);

%% Plotting

% Output Plot
figure();
plot(timeVector, output, 'LineWidth', 1.5);
hold on;
plot(timeVector, controlSignal, '--r', 'LineWidth', 1.5);
grid on;
title('Output Response', 'FontSize', 14);
legend('Obtained', 'Desired');
xlabel('Time (samples)');
ylabel('Output');
set(gca, 'FontSize', 12);

% Control Effort Plot
figure();
plot(timeVector, input, 'LineWidth', 1.5);
grid on;
title('Control Effort', 'FontSize', 14);
xlabel('Time (samples)');
ylabel('Control Signal');
set(gca, 'FontSize', 12);

% Accumulated Loss for Output and Control Signal
figure();

% Subplot for Accumulated Loss for Output
subplot(1, 2, 1);
plot(timeVector, accumulatedLossY, 'LineWidth', 1.5);
title('Accumulated Loss for Output (Y)', 'FontSize', 14);
xlabel('Time (samples)');
ylabel('Accumulated Loss');
grid on;
set(gca, 'FontSize', 12);

% Subplot for Accumulated Loss for Control Signal
subplot(1, 2, 2);
plot(timeVector, accumulatedLossU, 'LineWidth', 1.5);
title('Accumulated Loss for Control Signal (U)', 'FontSize', 14);
xlabel('Time (samples)');
ylabel('Accumulated Loss');
grid on;
set(gca, 'FontSize', 12);

% R parameters Plot
totalPlotRows = 3;
figure();
set(gcf, 'Position', [0 0 600 900]);
for i = 2:length(R_estimate) % First parameter of R is 1
    subplot(totalPlotRows, 2, i-1);
    plot(1:numberOfSamples, R_estimate_toPlot(:, i), 'LineWidth', 1.5, 'Color', [0.1 0.6 0.8]);
    hold on;
    plot(1:numberOfSamples, ones(numberOfSamples, 1) * R_real(i), '--', 'LineWidth', 1.5, 'Color', [0.8 0.1 0.2]);
    title(sprintf('R_%d', i-1), 'FontSize', 12);
    legend('Estimated', 'Real', 'Location', 'best');
    xlabel('Time (samples)');
    grid on;
    set(gca, 'FontSize', 10);
end

% S parameters Plot
end_of_r_plot = i-1;
for i = 1:length(S_estimate)
    subplot(totalPlotRows, 2, i + end_of_r_plot);
    plot(1:numberOfSamples, S_estimate_toPlot(:, i), 'LineWidth', 1.5, 'Color', [0.1 0.8 0.1]);
    hold on;
    plot(1:numberOfSamples, ones(numberOfSamples, 1) * S_real(i), '--', 'LineWidth', 1.5, 'Color', [0.8 0.1 0.1]);
    title(sprintf('S_%d', i), 'FontSize', 12);
    legend('Estimated', 'Real', 'Location', 'best');
    xlabel('Time (samples)');
    grid on;
    set(gca, 'FontSize', 10);
end

% T parameters Plot
end_of_S_plot = i + end_of_r_plot;
for i = 1:length(T_estimate)
    subplot(totalPlotRows, 2, i + end_of_S_plot);
    plot(1:numberOfSamples, T_toPlot(:, i), 'LineWidth', 1.5, 'Color', [0.6 0.3 0.9]);
    hold on;
    plot(1:numberOfSamples, ones(numberOfSamples, 1) * T_real(i), '--', 'LineWidth', 1.5, 'Color', [0.9 0.6 0.2]);
    title(sprintf('T_%d', i), 'FontSize', 12);
    legend('Estimated', 'Real', 'Location', 'best');
    xlabel('Time (samples)');
    grid on;
    set(gca, 'FontSize', 10);
end

% Theta parameters Plot
figure();
set(gcf, 'Position', [0 0 1200 600]);
for i = 1:4
    subplot(ceil(length(theta_hat_toPlot(1, :))/2), 2, i);
    plot(1:numberOfSamples, theta_hat_toPlot(:, i), 'LineWidth', 1.5, 'Color', [0.2 0.7 0.9]);
    hold on;
    plot(1:numberOfSamples, ones(numberOfSamples, 1) * theta_real(i), '--', 'LineWidth', 1.5, 'Color', [0.9 0.2 0.2]);
    title(sprintf('θ_%d', i), 'FontSize', 12);
    legend('Estimated', 'Real', 'Location', 'best');
    xlabel('Time (samples)');
    grid on;
    set(gca, 'FontSize', 10);
end
