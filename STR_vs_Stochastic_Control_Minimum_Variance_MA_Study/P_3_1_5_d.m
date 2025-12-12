%% Amirhossein Ehsani_810602159_ Sim 3_ Q1_IndirSTR _Colored Noise_Adaptive_Moving Average
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
numSamples = 1000;
timeVector = 0:numSamples-1;
noisePolynomial = [1, 0.3, 0.3].'; % Noise polynomial
degreeNoise = length(noisePolynomial);
noiseVariance = 0.1; % System noise variance
systemNoise = sqrt(noiseVariance) * randn(1, numSamples); % Generate system noise

% Initial conditions for parameter estimation
lenA = 3;
lenB = 2;
degA = lenA - 1;
degB = lenB;
A_estimate = [1, -1, 0.1]; % Initial A polynomial
B_estimate = [1, -0.6]; % Initial B polynomial
thetaNoise = [0.5, 0.5].'; % Initial conditions for noise polynomial
C_estimate = [1, thetaNoise.'];
d0 = degA - degB;

% Initial R and S polynomials
R_estimate = [1, -0.1]; 
S_estimate = [1, 0.1];

% Real R and S polynomials (calculated in the last part)
R_real = [1  -11.6339];
S_real = [13.0247   -2.1049];

% Skip instances for initial conditions
skipInstances = max(lenA, lenB);
totalParams = lenA - 1 + lenB;

% Initialize output and input
output = [];
output(1:skipInstances) = 0.1;
input = zeros(1, numSamples);

%% Adaptive Moving Average Control
theta_real = [A_discrete(2:end).'; B_discrete.'];
theta_real_toPlot = [A_discrete(2:end).'; B_discrete.'; noisePolynomial(2:end)];

% Initialize plotting matrices
R_estimate_toPlot = zeros(numSamples, length(R_estimate));
S_estimate_toPlot = zeros(numSamples, length(S_estimate));
theta_hat_toPlot = zeros(numSamples, length(theta_real_toPlot));

% Initialize ELS solver
els_solver = ELSClass(100 * eye(totalParams + length(thetaNoise)), 0.1 * ones(totalParams, 1), thetaNoise);

% Main loop for parameter estimation and control
for i = skipInstances:numSamples
    phi = [-output(i-1:-1:i-(lenA - 1)), input(i-1:-1:i-lenB)].';  
    noise_t = [systemNoise(i:-1:i-(degreeNoise-1))] * noisePolynomial;
    output(i) = phi.' * theta_real + noise_t;
    input(i) = (S_estimate * [-output(i:-1:i-(length(S_estimate)-1))].' - R_estimate(2:end) * [input(i-1:-1:i-(length(R_estimate)-1))].') / R_estimate(1);
    
    % Update parameter estimates using ELS
    theta_hat_new = els_solver.update_ELS(output(i), phi);
    A_estimate = [1, theta_hat_new(1:(lenA - 1)).'];
    B_estimate = theta_hat_new(lenA:totalParams).';
    C_not_modified = [1, theta_hat_new(totalParams+1:end).'];
    C_estimate = modify_polynomial(C_not_modified); % Modify polynomial if necessary
    
    % Solve Diophantine equation
    d = 2;
    Bplus = 1;
    A_dio = A_estimate;
    B_dio = B_estimate; 
    D_dio = conv([1, zeros(1, d-1)], C_estimate);
    D_dio = conv(D_dio, Bplus);
    [alpha, beta] = solve_diophantine_general(A_dio, B_dio, D_dio, 0);
    R_estimate = alpha; 
    S_estimate = beta; 
    if R_estimate(1) == 0
        R_estimate = R_estimate(2:end);
    end
    if S_estimate(1) == 0
        S_estimate = S_estimate(2:end);
    end
    
    % Store estimates for plotting
    theta_hat_toPlot(i, :) = theta_hat_new(1:end).';
    R_estimate_toPlot(i, :) = R_estimate;
    S_estimate_toPlot(i, :) = S_estimate;
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
desired_output = zeros(numSamples, 1);

% Output Plot
figure();
plot(timeVector, output, 'LineWidth', 1.5, 'Color', [0.1 0.6 0.8]);
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
totalPlotRows = ceil(length(R_estimate)/2) + ceil(length(S_estimate)/2);
figure();
set(gcf, 'Position', [0 0 800 900]);

% R parameters Plot
for i = 1:length(R_estimate)
    subplot(totalPlotRows, 2, i);
    plot(1:numSamples, R_estimate_toPlot(:, i), 'LineWidth', 1.5, 'Color', [0.1 0.6 0.8]);
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
for i = 1:length(S_estimate)
    subplot(totalPlotRows, 2, i + end_of_r_plot);
    plot(1:numSamples, S_estimate_toPlot(:, i), 'LineWidth', 1.5, 'Color', [0.2 0.7 0.3]);
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