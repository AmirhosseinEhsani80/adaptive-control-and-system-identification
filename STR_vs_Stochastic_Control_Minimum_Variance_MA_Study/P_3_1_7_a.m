%% Amirhossein Ehsani_810602159_ Sim 3_ Q1_IndirSTR _Colored Noise_Non Adaptive_Minimum Variance_Non Minimum Phase System
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
C = noisePolynomial.';
degreeNoise = length(noisePolynomial);
noiseVariance = 0.001; % System noise variance
systemNoise = sqrt(noiseVariance) * randn(1, numSamples); % Generate system noise

% Initial control input
controlInput = zeros(1, numSamples);

% Desired system parameters
lenA = 3;
lenB = 2;
degA = lenA - 1;
degB = lenB - 1;
d0 = degA - degB;

% Solve Diophantine equation
A_dio = A_discrete;
B_dio = B_discrete; 
D_dio = conv([1, zeros(1, d0-1)], C);
D_dio = conv(D_dio, modify_polynomial(B_discrete));
[alpha, beta] = solve_diophantine_general(A_dio, B_dio, D_dio, 0);
R = alpha; 
S = beta; 

if R(1) == 0
    R = R(2:end);
end
if S(1) == 0
    S = S(2:end);
end

% Skip instances for initial conditions
skipInstances = max(lenA, lenB);
totalParams = lenA - 1 + lenB;

% Initialize output and input
output = [];
output(1:skipInstances) = 0.1;
input = zeros(1, numSamples);

%% Non-Adaptive Minimum Variance NMP
theta_real = [A_discrete(2:end).'; B_discrete.'];
theta_hat_toPlot = zeros(numSamples, totalParams);

for i = skipInstances:numSamples
    phi = [-output(i-1:-1:i-(lenA - 1)), input(i-1:-1:i-lenB)].';  
    noise_t = [systemNoise(i:-1:i-(degreeNoise-1))] * noisePolynomial;
    output(i) = phi.' * theta_real + noise_t + B_discrete * controlInput(i:-1:i-(length(B_discrete)-1)).';  
    input(i) = (S * [-output(i:-1:i-(length(S)-1))].' - R(2:end) * [input(i-1:-1:i-(length(R)-1))].') / R(1);
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
