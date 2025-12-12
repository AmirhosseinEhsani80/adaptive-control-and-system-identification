%% Amirhossein Ehsani_810602159_ Sim 3_ Q1_IndirSTR _Colored Noise_Non Adaptive_Moving Average
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
C = noisePolynomial.';
degreeNoise = length(noisePolynomial);
systemNoiseVariance = 0.001; % System noise variance
systemNoise = sqrt(systemNoiseVariance) * randn(1, numberOfSamples); % Generate system noise

% Initial conditions for parameter estimation
len_desA = 3;
len_desB = 2;
degA = len_desA - 1;
degB = len_desB - 1;
d0 = degA - degB;
d = 2; % Degree of moving average
Bplus = 1; % B+ is set to 1

% Solve Diophantine equation
A_dio = A_discrete;
B_dio = B_discrete; 
D_dio = conv([1, zeros(1, d-1)], C);
D_dio = conv(D_dio, Bplus);
[alpha, beta] = solve_diophantine_general(A_dio, B_dio, D_dio, 0);
R = alpha; 
S = beta; 

% Skip instances for initial conditions
skipInstances = max(len_desA, len_desB);
totalParams = len_desA - 1 + len_desB;

%% Non-Adaptive Moving Average Control
output = [];
output(1:skipInstances) = 0.1;
input = zeros(1, numberOfSamples);
theta_real = [A_discrete(2:end).'; B_discrete.'];
theta_hat_toPlot = zeros(numberOfSamples, totalParams);

for i = skipInstances:numberOfSamples 
    phi = [-output(i-1:-1:i-(len_desA - 1)), input(i-1:-1:i-len_desB)].';  
    noise_t = [systemNoise(i:-1:i-(degreeNoise-1))] * noisePolynomial;
    output(i) = phi.' * theta_real + noise_t;  
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

%% Plots and Metrics
Mat = [meanU, meanY, varianceU, varianceY];
desired_output = zeros(numberOfSamples, 1);

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