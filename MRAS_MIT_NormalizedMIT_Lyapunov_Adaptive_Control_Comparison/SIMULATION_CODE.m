%% Amirhossein Ehsani_ Simulation 5
clc
clear all
close all
%% System Parameters

ab = 59;
R2 = (ab+24) * 0.01;

tau1 = 63.85;
tau2 = 1048.2575;

num_real = [R2/(tau1*tau2)];
den_real = [1, (tau1 + tau2 + R2) / (tau1*tau2), 1/ (tau1*tau2)];
sys = tf(num_real, den_real);

%% MIT Rule
theta = sym('theta', [3 1]);
syms b a1 a2 p t
y_closedloop(t) = (b * theta(1)) / (p^2 + (b*theta(3) + a1) * p  + (b*theta(2)+a2));

diff_theta1 = diff(y_closedloop, theta(1));
diff_theta2 = diff(y_closedloop, theta(2));
diff_theta3 = diff(y_closedloop, theta(3));


%% Desired System
zeta = 0.6;
wn = 0.04; 

sys_des = tf([wn^2],[1, 2*zeta*wn, wn^2]);
[num_des, den_des] = tfdata(sys_des, 'v');
