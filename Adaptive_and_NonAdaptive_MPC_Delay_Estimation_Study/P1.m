%% Amirhossein Ehsani Simulation 4 _Q1
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
B = (z-0.46)
sysD = z^-delay * B / A;
[numDiscrete, denDiscrete] = tfdata(sysD, 'v');

numeratorDiscrete = numDiscrete;
denominatorDiscrete = denDiscrete;

% Adjust the numerator if leading coefficients are zero
if numeratorDiscrete(1:3) == 0
    numeratorDiscrete = numeratorDiscrete(4:end);
end

%% Desired system specifications
overshootPercentage = 5;
riseTime = 1;

% Calculate damping ratio and natural frequency
zeta = cos(atan2(pi,-log(overshootPercentage/100)));
wn = 1.8/(riseTime);  

% Desired pole and zeros
desiredZero = -25;
desiredPole1 = -35;
desiredPole2 = -40;

% Gain calculation
gainK2 = -(desiredPole1 * desiredPole2) / (desiredZero);

% Continuous-time transfer functions
G1 = tf([wn^2], [1, 2*zeta*wn, wn^2]);
G2 = zpk([desiredZero], [desiredPole1, desiredPole2], gainK2);

% Desired discrete-time system
sys_desC = G1 * G2;
sys_desD = c2d(sys_desC, Ts, 'zoh');
[numDesiredDiscrete, denDesiredDiscrete] = tfdata(sys_desD, 'v');

numeratorDesiredPrime = numDesiredDiscrete; 
denominatorDesired = denDesiredDiscrete;

%% Without pole-zero cancellation

syms q;

% Adjust the coefficients for the Diophantine equation
Bplus = [numeratorDiscrete / numeratorDiscrete(1)];
Bminus = numeratorDiscrete(1);
numeratorDesiredPrime = numeratorDesiredPrime / numeratorDiscrete(1);
Ao = [1, 0, 0];

% Solve the Diophantine equation
[R_solved, S_solved, T] =  Diophantine(denominatorDesired, denominatorDiscrete, Ao, Bminus, Bplus, numeratorDesiredPrime);
Ac_solved = conv(denominatorDiscrete, R_solved) + [zeros(1, length(Ao)), conv(numeratorDiscrete, S_solved)]; % Depends on degree

%% Simulation
num_samples = 1000;
step_mag = 1;
t = 0:num_samples-1;
uc = (-1).^ceil(t/200);
uc(1) = -1;

% System output and control effort
y_output = filter(conv(numeratorDiscrete, T), Ac_solved, uc);
u_effort = filter(conv(denominatorDiscrete, T), Ac_solved, uc);

%% Plotting Results

% Output response plot
figure()
plot(t, y_output, 'LineWidth', 1.5, 'Color', [0 0.5 1])
hold on;
plot(t, uc, "--r", 'LineWidth', 1.5)
title("System Output Response", 'FontSize', 14)
grid on
xlabel('Sample time', 'FontSize', 12)
ylabel('Response', 'FontSize', 12)
legend({'System Output', 'Control Input'}, 'FontSize', 12)
set(gca, 'FontSize', 12)

% Control input plot
figure()
plot(t, u_effort, 'LineWidth', 1.5, 'Color', [0.5 0 0.5])
title("Control Effort", 'FontSize', 14)
grid on
xlabel('Sample time', 'FontSize', 12)
ylabel('Effort', 'FontSize', 12)
legend({'Control Effort'}, 'FontSize', 12)
set(gca, 'FontSize', 12)
