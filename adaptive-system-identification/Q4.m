%% Simulation_1 Q_4 Amirhossein Ehsani_810602159
clear all
close all
clc
%% System
% Parameters
a1 = 0.25;
a2 = -0.02;
a3 = 0.01;
a4 = 0.02;

b1 = 0.6;
b2 = 0.8;
b3 = 0.5;
b4 = 0.4;

c0 = 1;
c1 = 0.25;

sigma = 0.07;
N = 5000;
pt=ones(10,N);

% Noise:
e = sigma * randn(1, N);

% Input Signal:
q1 = input('What is the Input Signal for the system?\n\1)Impulse \n\2)Step \n\3)Sum of Two Sin \n\4)White Noise\n\');
q2 = input('What type of noise do we have?\n\1)White Noise \n\2)Colored Noise\n\');
q3 = input('What is the first Value for P?');
q4 = input('What is the first value for Theta?');
if q1 == 1
    u=zeros(N,1);
    u(1)=1;
    signal = "Impulse"
elseif q1 == 2
    u=ones(N,1);
    signal = "Step"
elseif q1 == 3
    u1 = sin(2*pi*pt);
    u2=sin(2*pi*5*pt);
    u=u1+u2;
    signal = "Sine1 + Sine2"
elseif q1 == 4
    u=sqrt(1.7)*randn(N,1);
    u=u-mean(u);
    signal = "White Noise"
end
if q2 == 1
    c1 = 0;
else 
    c1 = 0.3;
end

%% System Response

y(1) = e(1);
y(2) = a1*y(1)+b1*u(1)+e(2);
y(3) = a1*y(2)+a2*y(1)+b1*u(2)+b2*u(1)+e(3);
y(4) = a1*y(3)+a2*y(2)+a3*y(1)+b1*u(3)+b2*u(2)+b3*u(1)+e(4);
for t=5:N
    pt(1,t)=a1;
    pt(2,t)=a2;
    pt(3,t)=a3;
    pt(4,t)=a4;
    pt(5,t)=b1;
    pt(6,t)=b2;
    pt(7,t)=b3;
    pt(8,t)=b4;
    
    y(t)=a1*y(t-1)+a2*y(t-2)+a3*y(t-3)+a4*y(t-4)+b1*u(t-1)+b2*u(t-2)+b3*u(t-3)+b4*u(t-4)+e(t)+c1*e(t-1);
end

%% Estimation Using KF

phi = zeros(N,8);
phi(2,:) = [y(1) 0 0 0 u(1) 0 0 0];
phi(3,:) = [y(2) y(1) 0 0 u(2) u(1) 0 0];
phi(4,:) = [y(3) y(2) y(1) 0 u(3) u(2) u(1) 0];
for i=5:N
    phi(i,:) = [y(i-1) y(i-2) y(i-3) y(i-4) u(i-1) u(i-2) u(i-3) u(i-4)];
end

phi = phi';

p_kf(:,:,1)=q3*eye(8);
p_kf(:,:,2)=p_kf(:,:,1);
p_kf(:,:,3)=p_kf(:,:,1);
p_kf(:,:,4)=p_kf(:,:,1);

Teta_kf(:,1)=q4*[1;1;1;1;
                 1;1;1;1];
Teta_kf(:,2)=Teta_kf(:,1);
Teta_kf(:,3)=Teta_kf(:,1);
Teta_kf(:,4)=Teta_kf(:,1);

K_KF(:,1)=zeros(8,1);
K_KF(:,2)=zeros(8,1);
K_KF(:,3)=zeros(8,1);
K_KF(:,4)=zeros(8,1);

for t=5:N
    K_KF(:,t)=p_kf(:,:,t-1)*phi(:,t)*inv(1+phi(:,t)'*p_kf(:,:,t-1)*phi(:,t));
    p_kf(:,:,t)=p_kf(:,:,t-1)-p_kf(:,:,t-1)*phi(:,t)*inv(1+phi(:,t)'*p_kf(:,:,t-1)*phi(:,t))*phi(:,t)'*p_kf(:,:,t-1)+0.05;
    epsilon(t)=y(t)-phi(:,t)'*Teta_kf(:,t-1);
    Teta_kf(:,t)=Teta_kf(:,t-1)+K_KF(:,t)*epsilon(t);
end

theta0=[a1 a2 a3 a4 b1 b2 b3 b4]';
Teta_kf(:,N);
error_KF=norm(Teta_kf(:,N)-theta0)/5*100;
RMSD_KF=sqrt(sum((Teta_kf(:,N)-theta0).^2)/5);


%% Kalman Filter (KF) Plot

%PLOTS% Define a color palette for better visuals
colors = lines(10); % This will generate 7 distinct colors

figure(1);
plot(1:N, Teta_kf(1,:), 'Color', colors(1,:), 'LineWidth', 1.5); hold on;
plot(1:N, Teta_kf(2,:), 'Color', colors(2,:), 'LineWidth', 1.5);
plot(1:N, Teta_kf(3,:), 'Color', colors(3,:), 'LineWidth', 1.5);
plot(1:N, Teta_kf(4,:), 'Color', colors(4,:), 'LineWidth', 1.5);
plot(1:N, Teta_kf(5,:), 'Color', colors(5,:), 'LineWidth', 1.5);
plot(1:N, Teta_kf(6,:), 'Color', colors(6,:), 'LineWidth', 1.5);
plot(1:N, Teta_kf(7,:), 'Color', colors(7,:), 'LineWidth', 1.5);
plot(1:N, Teta_kf(8,:), 'Color', colors(8,:), 'LineWidth', 1.5);

plot(pt', 'b:', 'LineWidth', 1.5);
hold off;
title('Estimated Parameters from Kalman Filter', 'FontSize', 14);
legend('a1','a2','a3','a4','b1','b2','b3','b4');
xlabel('Sample');
ylabel('Parameter Value');
grid on;

%% RLS

P(:,:,1)=q3*eye(8);
P(:,:,2)=P(:,:,1);
P(:,:,3)=P(:,:,1);
P(:,:,4)=P(:,:,1);

Teta_rls=q4*ones(8,N);
Teta_rls(:,1:4)= 0;
K = zeros(8,N);
K(:,1:4)= 0;
for t=5:N
    p_inv = inv(P(:,:,t-1))+phi(:,t)*phi(:,t)';
    P(:,:,t) = inv(p_inv);
    K(:,t)=P(:,:,t)*phi(:,t);
    Teta_rls(:,t)=Teta_rls(:,t-1)+K(:,t)*(y(t)-phi(:,t)'*Teta_rls(:,t-1));
end

theta0=[a1 a2 a3 a4 b1 b2 b3 b4]';
Teta_rls(:,N);
error_RLS=norm(Teta_rls(:,N)-theta0)/5*100;
RMSD_RLS=sqrt(sum((Teta_rls(:,N)-theta0).^2)/5);

%% RLS Plot
figure(2);
plot(1:N, Teta_rls(1,:), 'Color', colors(1,:), 'LineWidth', 1.5); hold on;
plot(1:N, Teta_rls(2,:), 'Color', colors(2,:), 'LineWidth', 1.5);
plot(1:N, Teta_rls(3,:), 'Color', colors(3,:), 'LineWidth', 1.5);
plot(1:N, Teta_rls(4,:), 'Color', colors(4,:), 'LineWidth', 1.5);

plot(1:N, Teta_rls(5,:), 'Color', colors(5,:), 'LineWidth', 1.5);
plot(1:N, Teta_rls(6,:), 'Color', colors(6,:), 'LineWidth', 1.5);
plot(1:N, Teta_rls(7,:), 'Color', colors(7,:), 'LineWidth', 1.5);
plot(1:N, Teta_rls(8,:), 'Color', colors(8,:), 'LineWidth', 1.5);

plot(pt', 'b:', 'LineWidth', 1.5);
hold off;
title('Estimated Parameters from RLS', 'FontSize', 14);
legend('a1','a2','a3','a4','b1','b2','b2','b3', 'Location', 'best');
xlabel('Sample');
ylabel('Parameter Value');
grid on;



%% Results:
method.Name = {'Real System';'RLS';'KF'};
method.a1 = [a1;Teta_rls(1,N);Teta_kf(1,N)]; 
method.a2 = [a2;Teta_rls(2,N);Teta_kf(2,N)];
method.a3 = [a3;Teta_rls(3,N);Teta_kf(3,N)];
method.a4 = [a4;Teta_rls(4,N);Teta_kf(4,N)];

method.b1 = [b1;Teta_rls(5,N);Teta_kf(5,N)];
method.b2 = [b2;Teta_rls(6,N);Teta_kf(6,N)];
method.b3 = [b3;Teta_rls(7,N);Teta_kf(7,N)];
method.b4 = [b4;Teta_rls(8,N);Teta_kf(8,N)];

method.RMSD_error = [0;RMSD_RLS;RMSD_KF];
T = struct2table(method)


