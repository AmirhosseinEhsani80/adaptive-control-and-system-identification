%% Amirhossein Ehsani_810602159_ Sim 3_ Q1_IndirSTR _Colored Noise_Adaptive_Moving Average_Non Minimum Phase System
clc;
clear all;
close all;

%% Discretization of the Model
z = tf('z', 0.1); % Discrete-time transfer function with sampling time 0.1 seconds
sysDiscrete = (z + 2.4)/(z^2 - 0.8*z + 0.1);
[numDiscrete, denDiscrete] = tfdata(sysDiscrete, 'v'); % Extract numerator and denominator
samplingTime = 0.1; % Sampling time

% Normalize B and A polynomials
B_discrete = numDiscrete;
A_discrete = denDiscrete;
if B_discrete(1) == 0
    B_discrete = B_discrete(2:end);
end

%% Initial Conditions
rng(50);
numSamples = 1000;
timeVector = 0:numSamples-1;
noisePolynomial = [1, 0.4, 0.1].'; % Noise polynomial
degreeNoise = length(noisePolynomial);
noiseVariance = 0.001; % System noise variance
systemNoise = sqrt(noiseVariance) * randn(1, numSamples); % Generate system noise

% Initial control input
controlInput = zeros(1, numSamples);

% Desired system parameters
lenA = 3;
lenB = 2;
degA = lenA - 1;
degB = lenB;
A_estimated = [1, -1, 0.1]; % Initial estimated A polynomial
B_estimated = [1, -0.6]; % Initial estimated B polynomial
theta_epsilon_zero = [0.5, 0.5].'; % Initial conditions for noise polynomial
C_estimated = [1, theta_epsilon_zero.'];


d0 = degA - degB;

% Initial R and S polynomials
R_estimated = [1, -0.1]; 
S_estimated = [1, 0.1];

% Real R and S polynomials (calculated previously)
R_real = [1.0000    0.8929];
S_real = [0.3203   -0.0349];

% Skip instances for initial conditions
skipInstances = max(lenA, lenB);
totalParams = lenA - 1 + lenB;

% Initialize output and input
output = [];
output(1:skipInstances) = 0.1;
input = zeros(1, numSamples);

%% Adaptive Moving Average NMP
theta_real = [A_discrete(2:end).'; B_discrete.'];
theta_real_toPlot = [A_discrete(2:end).'; B_discrete.'; noisePolynomial(2:end)];

% Initialize plotting matrices
R_estimated_toPlot = zeros(numSamples, length(R_estimated));
S_estimated_toPlot = zeros(numSamples, length(S_estimated));
theta_hat_toPlot = zeros(numSamples, length(theta_real_toPlot));

% Initialize ELS solver
els_solver = ELSClass(100 * eye(totalParams + length(theta_epsilon_zero)), 0.1 * ones(totalParams, 1), theta_epsilon_zero);

% Main loop for parameter estimation and control
for i = skipInstances:numSamples 
    phi = [-output(i-1:-1:i-(lenA - 1)), input(i-1:-1:i-lenB)].';  
    noise_t = [systemNoise(i:-1:i-(degreeNoise-1))] * noisePolynomial;
    output(i) = phi.' * theta_real + noise_t + B_discrete * controlInput(i:-1:i-(length(B_discrete)-1)).';
    input(i) = (S_estimated * [-output(i:-1:i-(length(S_estimated)-1))].' - R_estimated(2:end) * [input(i-1:-1:i-(length(R_estimated)-1))].') / R_estimated(1);
    
    % Update parameter estimates using ELS
    theta_hat_new = els_solver.update_ELS(output(i), phi);
    A_estimated = [1, theta_hat_new(1:(lenA - 1)).'];
    B_estimated = theta_hat_new(lenA:totalParams).';
    C_not_modified = [1, theta_hat_new(totalParams+1:end).'];
    C_estimated = modify_polynomial(C_not_modified); % Modify polynomial if necessary
    
    % Factorize B polynomial
    [Bplus, Bminus] = factor_polynomial(B_discrete);
    d = degA - (length(Bplus) - 1);
    A_dio = A_estimated;
    B_dio = B_estimated; 
    D_dio = conv([1, zeros(1, d-1)], C_estimated);
    D_dio = conv(D_dio, Bplus);
    
    % Solve Diophantine equation
    [alpha, beta] = solve_diophantine_general(A_dio, B_dio, D_dio, 0);
    R_estimated = alpha; 
    S_estimated = beta; 
    if R_estimated(1) == 0
        R_estimated = R_estimated(2:end);
    end
    if S_estimated(1) == 0
        S_estimated = S_estimated(2:end);
    end
    
    % Store estimates for plotting
    theta_hat_toPlot(i, :) = theta_hat_new(1:end).';
    R_estimated_toPlot(i, :) = R_estimated;
    S_estimated_toPlot(i, :) = S_estimated;
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

Mat = [meanU, meanY, varianceU, varianceY];

%% Plots and Metrics
desired_output = zeros(numSamples, 1);
% Output Plot
figure();
plot(timeVector, output, 'LineWidth', 1.5, 'Color', [0.1 0.6 0.8]);
grid on;
hold on;
plot(timeVector, desired_output, '--r', 'LineWidth', 1.5);
grid on;
title('Output Response', 'FontSize', 14);
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
totalPlotRows = ceil(length(R_estimated)/2) + ceil(length(S_estimated)/2);
figure();
set(gcf, 'Position', [0 0 800 900]);

% R parameters Plot
for i = 1:length(R_estimated)
    subplot(totalPlotRows, 2, i);
    plot(1:numSamples, R_estimated_toPlot(:, i), 'LineWidth', 1.5, 'Color', [0.1 0.6 0.8]);
    hold on;
    plot(1:numSamples, ones(numSamples, 1) * R_real(i), '--k', 'LineWidth', 1.5);
    title(sprintf('R_%d', i), 'FontSize', 12);
    legend('Predicted', 'Real', 'Location', 'best');
    xlabel('Time (samples)');
    grid on;
    set(gca, 'FontSize', 10);
end

% S parameters Plot
end_of_r_plot = i;
for i = 1:length(S_estimated)
    subplot(totalPlotRows, 2, i + end_of_r_plot);
    plot(1:numSamples, S_estimated_toPlot(:, i), 'LineWidth', 1.5, 'Color', [0.2 0.7 0.3]);
    hold on;
    plot(1:numSamples, ones(numSamples, 1) * S_real(i), '--k', 'LineWidth', 1.5);
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
    plot(1:numSamples, theta_hat_toPlot(:, i), 'LineWidth', 1.5, 'Color', [0.4 0.4 0.8]);
    hold on;
    plot(1:numSamples, ones(numSamples, 1) * theta_real_toPlot(i), '--k', 'LineWidth', 1.5);
    title(sprintf('θ_%d', i), 'FontSize', 12);
    legend('Predicted', 'Real', 'Location', 'best');
    xlabel('Time (samples)');
    grid on;
    set(gca, 'FontSize', 10);
end
