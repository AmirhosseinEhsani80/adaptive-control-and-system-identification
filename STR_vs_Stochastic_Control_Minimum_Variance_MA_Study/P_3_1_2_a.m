%% Amirhossein Ehsani_810602159_ Sim 3_ Q1_IndirSTR _Without Noise
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

% Discrete system
z = tf('z', 0.1); % Sampling time of 0.1 seconds
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

cRatio = cos(atan2(pi,-log(desiredOvershoot))); % c ratio
naturalFrequency = 1.8/(desiredRiseTime); % Natural frequency

desiredZero = -20; % Zero of the desired system
gainFactor = -1/desiredZero;

G1_desired = tf([naturalFrequency^2],[1, 2*cRatio*naturalFrequency, naturalFrequency^2]);
G2_desired = zpk([desiredZero],[],gainFactor);
sysDesiredC = G1_desired*G2_desired; % Desired continuous system
sysDesiredD = c2d(sysDesiredC, samplingTime, 'zoh'); % Discretize desired system
[numDesiredD, denDesiredD] = tfdata(sysDesiredD, 'v');

B_prime_main = numDesiredD; 
A_m = denDesiredD;

%% Indirect STR with Zero Cancellation 

% Input properties
numberOfSamples = 1000;
stepMagnitude = 1;
timeVector = 0:numberOfSamples-1;
controlSignal = (-1).^-ceil(timeVector/200);
controlSignal(1) = -1;

% Initial conditions for estimated polynomials
lengthA = 3;
lengthB = 2;
A_estimate = [1,-1,0.15]; % Initial A polynomial
B_estimate = [.9,-0.69]; % Initial B polynomial

% Initial parameters for RLS
theta_init = [0.1,0.1].';
C_estimate = [1, theta_init.']; % C polynomial must be monic

% Normalize B polynomial
B_plus = [B_estimate/B_estimate(1)];
B_minus = B_estimate(1);
B_prime = B_prime_main/B_estimate(1);

% Initial R, S, T, and A0 polynomials
R_estimate = [1,-0.1]; 
S_estimate = [1,0.1];
T_estimate = [1.1,0];
A_0 = [1];

% Initial conditions for output and input
initialSamples = max(lengthA, lengthB);
totalParams = lengthA - 1 + lengthB;
output = zeros(1, numberOfSamples);
input = zeros(1, numberOfSamples);

% Actual system parameters
theta_actual = [A_discrete(2:end).'; B_discrete.'];

% For plotting
R_plot = zeros(numberOfSamples, length(R_estimate));
S_plot = zeros(numberOfSamples, length(S_estimate));
T_plot = zeros(numberOfSamples, length(T_estimate));
theta_plot = zeros(numberOfSamples, totalParams);

% RLS initial parameters
theta_init = [0.5, 0.5].';
els_solver = ELSClass(100 * eye(totalParams + length(theta_init)), 0.1 * ones(totalParams,1), theta_init);

% Main loop for parameter estimation and control
for i = initialSamples:numberOfSamples
    phi = [-output(i-1:-1:i-(lengthA - 1)), input(i-1:-1:i-lengthB)].'; 
    output(i) = phi.' * theta_actual; % Output based on real system parameters
    
    input(i) = T_estimate * [controlSignal(i:-1:i-(length(T_estimate)-1))].' + S_estimate * [-output(i:-1:i-(length(S_estimate)-1))].' - R_estimate(2:end) * [input(i-1:-1:i-(length(R_estimate)-1))].';

    theta_new = els_solver.update_ELS(output(i), phi); % Update estimated parameters

    A_estimate = theta_new(1:(lengthA - 1)).';
    B_estimate = theta_new(lengthA:totalParams).';
    
    B_plus = [B_estimate/B_estimate(1)];
    B_minus = B_estimate(1);
    B_prime = B_prime_main/B_estimate(1);
    A_0 = [1];

    [R_estimate, S_estimate, T_estimate] = Diophantine(A_m, [1,A_estimate], A_0, B_minus, B_plus, B_prime); % Solve Diophantine equations

    theta_plot(i, :) = theta_new(1:totalParams).';
    R_plot(i, :) = R_estimate;
    S_plot(i, :) = S_estimate;
    T_plot(i, :) = T_estimate;
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

accumulatedLossY_indirect = accumulatedLossY;

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
for i = 1:length(R_estimate)
    subplot(totalPlotRows,2,i);
    plot(1:numberOfSamples, R_plot(:,i), 'LineWidth', 1.5, 'Color', [0.1 0.6 0.8]);
    hold on;
    plot(1:numberOfSamples, ones(numberOfSamples,1) * R_estimate(i), '--', 'LineWidth', 1.5, 'Color', [0.8 0.1 0.2]);
    title(sprintf('R_%d', i), 'FontSize', 12);
    legend('Estimated', 'Real', 'Location', 'best');
    xlabel('Time (samples)');
    grid on;
    set(gca, 'FontSize', 10);
end

% S parameters Plot
for i = 1:length(S_estimate)
    subplot(totalPlotRows,2,i + 2);
    plot(1:numberOfSamples, S_plot(:,i), 'LineWidth', 1.5, 'Color', [0.1 0.8 0.1]);
    hold on;
    plot(1:numberOfSamples, ones(numberOfSamples,1) * S_estimate(i), '--', 'LineWidth', 1.5, 'Color', [0.8 0.1 0.1]);
    title(sprintf('S_%d', i), 'FontSize', 12);
    legend('Estimated', 'Real', 'Location', 'best');
    xlabel('Time (samples)');
    grid on;
    set(gca, 'FontSize', 10);
end

% T parameters Plot
for i = 1:length(T_estimate)
    subplot(totalPlotRows,2,i + 4);
    plot(1:numberOfSamples, T_plot(:,i), 'LineWidth', 1.5, 'Color', [0.6 0.3 0.9]);
    hold on;
    plot(1:numberOfSamples, ones(numberOfSamples,1) * T_estimate(i), '--', 'LineWidth', 1.5, 'Color', [0.9 0.6 0.2]);
    title(sprintf('T_%d', i), 'FontSize', 12);
    legend('Estimated', 'Real', 'Location', 'best');
    xlabel('Time (samples)');
    grid on;
    set(gca, 'FontSize', 10);
end

% Theta parameters Plot
figure();
set(gcf, 'Position', [0 0 1200 600]);
for i = 1:length(theta_plot(1,:))
    subplot(ceil(length(theta_plot(1,:))/2), 2, i);
    plot(1:numberOfSamples, theta_plot(:,i), 'LineWidth', 1.5, 'Color', [0.2 0.7 0.9]);
    hold on;
    plot(1:numberOfSamples, ones(numberOfSamples,1) * theta_actual(i), '--', 'LineWidth', 1.5, 'Color', [0.9 0.2 0.2]);
    title(sprintf('θ_%d', i), 'FontSize', 12);
    legend('Estimated', 'Real', 'Location', 'best');
    xlabel('Time (samples)');
    grid on;
    set(gca, 'FontSize', 10);
end
