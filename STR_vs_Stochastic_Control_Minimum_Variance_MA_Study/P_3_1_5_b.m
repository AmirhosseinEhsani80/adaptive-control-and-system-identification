%% Amirhossein Ehsani_810602159_ Sim 3_ Q1_IndirSTR _Colored Noise_Adaptive_Minimum Variance
clc;
clear all;
close all;

%% Discretization of the Model
s = tf('s');
m = 10;
g = 9.81;
c = 5;
k = 100;
Gc = (m*g)/(2*(s+c/k)*(3*m/2*s^2+c*s+2*k)); % Continuous-time transfer function

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

%% Initial Conditions
rng(50);
numberOfSamples = 1000;
stepMagnitude = 1;
timeVector = 0:numberOfSamples-1;
noisePolynomial = [1, 0.3, 0.3].'; % Noise polynomial
degreeNoise = length(noisePolynomial);
systemNoiseVariance = 0.001; % System noise variance
systemNoise = sqrt(systemNoiseVariance) * randn(1, numberOfSamples); % Generate system noise

% Initial conditions for parameter estimation
len_desA = 3;
len_desB = 2;
degA = len_desA - 1;
degB = len_desB;
A_estimate = [1, -1, 0.1]; % Initial A polynomial
B_estimate = [1, -0.6]; % Initial B polynomial
theta_epsilon_zero = [0.1, 0.1].'; % Initial conditions for noise polynomial
C_estimate = [1, theta_epsilon_zero.'];
d0 = degA - degB;

% Initial R and S polynomials
R_calculated = [1, -0.1]; 
S_calculated = [1, 0.1];

% Real R and S polynomials (calculated in the last part)
R_real = [1, -0.7];
S_real = [1.36, -0.028];

% Skip instances for initial conditions
skipInstances = max(len_desA, len_desB);
totalParams = len_desA - 1 + len_desB;

% Initialize output and input
output = [];
output(1:skipInstances) = 0.1;
input = zeros(1, numberOfSamples);

%% Adaptive Minimum Variance Control
theta_real = [A_discrete(2:end).'; B_discrete.'];
theta_real_toPlot = [A_discrete(2:end).'; B_discrete.'; noisePolynomial(2:end)];

% Initialize plotting matrices
R_calculated_toPlot = zeros(numberOfSamples, length(R_calculated));
S_calculated_toPlot = zeros(numberOfSamples, length(S_calculated));
theta_hat_toPlot = zeros(numberOfSamples, length(theta_real_toPlot));

% Initialize ELS solver
els_solver = ELSClass(100 * eye(totalParams + length(theta_epsilon_zero)), 0.1 * ones(totalParams, 1), theta_epsilon_zero);

% Main loop for parameter estimation and control
for i = skipInstances:numberOfSamples
    phi = [-output(i-1:-1:i-(len_desA - 1)), input(i-1:-1:i-len_desB)].';  
    noise_t = [systemNoise(i:-1:i-(degreeNoise-1))] * noisePolynomial;
    output(i) = phi.' * theta_real + noise_t;
    input(i) = (S_calculated * [-output(i:-1:i-(length(S_calculated)-1))].' - R_calculated(2:end) * [input(i-1:-1:i-(length(R_calculated)-1))].') / R_calculated(1);
    
    % Update parameter estimates using ELS
    theta_hat_new = els_solver.update_ELS(output(i), phi);
    A_estimate = [1, theta_hat_new(1:(len_desA - 1)).'];
    B_estimate = theta_hat_new(len_desA:totalParams).';
    C_not_modified = [1, theta_hat_new(totalParams+1:end).'];
    C_estimate = modify_polynomial(C_not_modified); % Modify polynomial if necessary
    
    % Solve Diophantine equation
    A_dio = A_estimate;
    B_dio = B_estimate; 
    D_dio = conv([1, zeros(1, d0-1)], C_estimate);
    D_dio = conv(D_dio, B_estimate);
    [alpha, beta] = solve_diophantine_general(A_dio, B_dio, D_dio, 0);
    R_calculated = alpha; 
    S_calculated = beta; 
    if R_calculated(1) == 0
        R_calculated = R_calculated(2:end);
    end
    if S_calculated(1) == 0
        S_calculated = S_calculated(2:end);
    end
    
    % Store estimates for plotting
    theta_hat_toPlot(i, :) = theta_hat_new(1:end).';
    R_calculated_toPlot(i, :) = R_calculated;
    S_calculated_toPlot(i, :) = S_calculated;
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

%% Plots and Metrics
Mat = [meanU, meanY, varianceU, varianceY];
desired_output = zeros(1000, 1);

% Output Plot
figure();
plot(timeVector, output, 'LineWidth', 1.5);
hold on;
plot(timeVector, desired_output, '--r', 'LineWidth', 1.5);
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
subplot(1, 2, 1);
plot(timeVector, accumulatedLossY, 'LineWidth', 1.5);
title('Y Accumulated Loss', 'FontSize', 14);
xlabel('Time (samples)');
ylabel('Accumulated Loss');
grid on;
set(gca, 'FontSize', 12);

subplot(1, 2, 2);
plot(timeVector, accumulatedLossU, 'LineWidth', 1.5);
title('U Accumulated Loss', 'FontSize', 14);
xlabel('Time (samples)');
ylabel('Accumulated Loss');
grid on;
set(gca, 'FontSize', 12);

% R and S parameters Plot
totalPlotRows = ceil(length(R_calculated)/2) + ceil(length(S_calculated)/2);
figure();
set(gcf, 'Position', [0 0 800 900]);

% R parameters Plot
for i = 1:length(R_calculated)
    subplot(totalPlotRows, 2, i);
    plot(1:numberOfSamples, R_calculated_toPlot(:, i), 'LineWidth', 1.5, 'Color', [0.1 0.6 0.8]);
    hold on;
    plot(1:numberOfSamples, ones(numberOfSamples, 1) * R_real(i), '--k', 'LineWidth', 1.5);
    title(sprintf('R_%d', i), 'FontSize', 12);
    legend('Predicted', 'Real', 'Location', 'best');
    xlabel('Time (samples)');
    grid on;
    set(gca, 'FontSize', 10);
end

% S parameters Plot
end_of_r_plot = i;
for i = 1:length(S_calculated)
    subplot(totalPlotRows, 2, i + end_of_r_plot);
    plot(1:numberOfSamples, S_calculated_toPlot(:, i), 'LineWidth', 1.5, 'Color', [0.2 0.7 0.3]);
    hold on;
    plot(1:numberOfSamples, ones(numberOfSamples, 1) * S_real(i), '--k', 'LineWidth', 1.5);
    title(sprintf('S_%d', i), 'FontSize', 12);
    legend('Predicted', 'Real', 'Location', 'best');
    xlabel('Time (samples)');
    grid on;
    set(gca, 'FontSize', 10);
end

% Theta parameters Plot
figure();
set(gcf, 'Position', [0 0 1200 600]);
for i = 1:length(theta_hat_toPlot(1, :))
    subplot(ceil(length(theta_hat_toPlot(1, :))/2), 2, i);
    plot(1:numberOfSamples, theta_hat_toPlot(:, i), 'LineWidth', 1.5, 'Color', [0.4 0.4 0.8]);
    hold on;
    plot(1:numberOfSamples, ones(numberOfSamples, 1) * theta_real_toPlot(i), '--k', 'LineWidth', 1.5);
    title(sprintf('θ_%d', i), 'FontSize', 12);
    legend('Predicted', 'Real', 'Location', 'best');
    xlabel('Time (samples)');
    grid on;
    set(gca, 'FontSize', 10);
end