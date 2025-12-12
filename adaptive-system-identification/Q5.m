%% Simulation_1 Q_5 Amirhossein Ehsani_810602159
clc;
clear all;
close all;

%% System
a12 = 2.2;
a21 = -4;
a22 = -0.34;
N = 200;
x1(1) = 0;
x2(1) = 0;
sigma = 0.5;
e = sigma * randn(1, N);

q1 = input('What is the first value of P:');
q2 = input('What is the first value of theta:');

for i=2:N
    pt(1,i)=a12;
    pt(2,i)=a21;
    pt(3,i)=a22;
    x1(i)=a12*x2(i-1)+e(i);
    x2(i)=a21*sin(x1(i-1))+a22*x2(i-1)+e(i);
end

%% RLS_1
% Phi 
Phi = zeros(N,1);
Phi(1,:)=[0];
for i=2:N
    Phi(i)=x2(i-1);
end
phi = Phi';

% Initialize P
P(:,:,1) = q1*eye(1);

T_RLS = q2*ones(1,N);
T_RLS(:,1) = 0;
K = zeros(1,N);
K(:,1) = 0;
y = x1;
for t=2:N
    p_inv = inv(P(:,:,(t-1))) + Phi(t,:)' * Phi(t,:);
    P(:,:,t) = inv(p_inv);
    K(:,t) = P(:,:,t) * Phi(t,:)';
    T_RLS(:,t) = T_RLS(:,t-1) + K(:,t) * (y(t) - Phi(t,:) * T_RLS(:,t-1));
end

true_parameters = [a12];
% Error calculations
error_RLS = norm(T_RLS(:,N) - true_parameters) / 12 * 100;
RMSD_RLS = sqrt(sum((T_RLS(:,N) - true_parameters).^2) / 12);
% Display the final error 
fprintf('Final RLS Error: %f%%\n', error_RLS);
fprintf('RMSD for RLS: %f\n', RMSD_RLS);

%% Plotting Results for RLS 1
colors = lines(14); % This will generate 7 distinct colors

% RLS Plot
figure(1);
plot(1:N, T_RLS(1,:), 'Color', colors(1,:), 'LineWidth', 1.5); hold on;
plot(pt(1,:)', 'r:', 'LineWidth', 1.5);
hold off;
title('Estimated Parameters from RLS', 'FontSize', 14);
legend('a12','real parameter', 'Location', 'best');
xlabel('Sample');
ylabel('Parameter Value');
grid on;

%% RLS 2
clear all
% Sys
a12 = 2.2;
a21 = -4;
a22 = -0.34;
N = 1000;
x1(1) = 0;
x2(1) = 0;
sigma = 0.2;
e = sigma * randn(1, N);

q1 = input('What is the first value of P:');
q2 = input('What is the first value of theta:');

for i=2:N
    pt(1,i)=a12;
    pt(2,i)=a21;
    pt(3,i)=a22;
    x1(i)=a12*x2(i-1)+e(i);
    x2(i)=a21*sin(x1(i-1))+a22*x2(i-1)+e(i);
end

% Phi 
Phi = zeros(N,2);
Phi(1,:)=[0 0];
Phi(2,:) = [sin(x1(1)) x2(1)];
for i=3:N
    Phi(i,:)=[sin(x1(i-1)) x2(i-1)];
end
phi = Phi';

% Initialize P

P(:,:,1) = q1*eye(2);
P(:,:,2) = P(:,:,1);

T_RLS = q2*ones(2,N);
T_RLS(:,1) = 0;
T_RLS(:,2) = 0;

K = zeros(2,N);
K(:,1) = 0;
K(:,2) = 0;
y = x2;
for t=2:N
    p_inv = inv(P(:,:,(t-1))) + Phi(t,:)' * Phi(t,:);
    P(:,:,t) = inv(p_inv);
    K(:,t) = P(:,:,t) * Phi(t,:)';
    T_RLS(:,t) = T_RLS(:,t-1) + K(:,t) * (y(t) - Phi(t,:) * T_RLS(:,t-1));
end

true_parameters = [a21, a22];
% Error calculations
error_RLS = norm(T_RLS(:,N) - true_parameters) / 12 * 100;
RMSD_RLS = sqrt(sum((T_RLS(:,N) - true_parameters).^2) / 12);
% Display the final error 
fprintf('Final RLS Error: %f%%\n', error_RLS);
fprintf('RMSD for RLS: %f\n', RMSD_RLS);
T_RLS_2=T_RLS;
%% Plotting Results for RLS 2
colors = lines(14); % This will generate 7 distinct colors

% RLS Plot
figure(2);
plot(1:N, T_RLS(1,:), 'Color', colors(1,:), 'LineWidth', 1.5); hold on;
plot(1:N, T_RLS(2,:), 'Color', colors(1,:), 'LineWidth', 1.5); hold on;

plot(pt(2,:)', 'r:', 'LineWidth', 1.5);
plot(pt(3,:)', 'r:', 'LineWidth', 1.5);
hold off;
title('Estimated Parameters from RLS', 'FontSize', 14);
legend('a21','a22','real parameter', 'Location', 'best');
xlabel('Sample');
ylabel('Parameter Value');
grid on;
N1=N;
N2=N;

