%% Amirhossein Ehsani - 810602159 - Q1_2
clc
clear all 
close all
%% Discretization of the Model

% Define the transfer function variable
s = tf('s');

% Continuous-time system transfer function
sysC = (s+0.5)*(s+7)/((s+1)*(s+4.1)*(s-2));

% Calculate the bandwidth of the continuous-time system
BW = bandwidth(sysC); 

% Discretization ratio and sampling frequency
disc_ratio = 10;
fs = BW * disc_ratio/(2*pi);
Ts = 1/fs; % Sampling time

% Discretize the continuous-time system using zero-order hold
sysD = c2d(sysC, Ts, 'zoh');
[numD, denD] = tfdata(sysD, 'v'); % Get numerator and denominator

B = numD;
A = denD;

%% Desired System Specifications

% Desired overshoot and settling time
overshoot = 10;
settling_time = 3;

% Calculate damping ratio and natural frequency
zeta = cos(atan2(pi,-log(overshoot/100)));
wn = 4/(zeta*settling_time); 

% Adjust damping ratio and natural frequency
zeta = 0.8;
wn = 4/(zeta*settling_time); 

% Desired pole and gain
p3 = -70;
k2 = -p3;

% Desired continuous-time system
G1 = tf([wn^2],[1, 2*zeta*wn, wn^2]);
G2 = zpk([-1,-3],[p3],k2);

sys_desC = G1*G2;
sys_desD = c2d(sys_desC, Ts, 'zoh'); % Discretize the desired system
[num_desD, den_desD] = tfdata(sys_desD, 'v');

Bm_prim2 = num_desD; 
Am = den_desD;

%% Without Pole-Zero Cancellation

syms q;

% Define polynomials
Bplus = [1];
Bminus = B;
Ao = [1,0,0];
Bm = conv(Bplus, Bminus);
BmPrime = (sum(Am)/sum(Bm)); 

% Solve the Diophantine equation
[R_solved, S_solved, T] = Diophantine(Am, A, Ao, Bminus, Bplus, BmPrime);
Ac_solved = conv(A, R_solved) + conv(B, S_solved);

%% Simulation

% Number of samples and input signal
num_samples = 1000;
step_mag = 1;
t = 0:num_samples-1;
uc = (-1).^ceil(t/200);

% Compute the system response
y = filter(conv(B, T), Ac_solved, uc);
u = filter(conv(A, T), Ac_solved, uc);

% Plot system response
figure()
plot(t, y, 'b', 'LineWidth', 1.5) % System response in blue
hold on;
plot(t, uc, '--r', 'LineWidth', 1.5) % Reference input in dashed red
grid on
title('System Response to Step Input')
xlabel('Time (samples)')
ylabel('Output')
legend('System Output', 'Reference Input')
set(gca, 'FontSize', 12)
hold off;

% Plot control input
figure()
plot(t, u, 'B', 'LineWidth', 1.5) % Control input in green
grid on
title('Control Effort')
xlabel('Time (samples)')
ylabel('Control Signal')
set(gca, 'FontSize', 12)
