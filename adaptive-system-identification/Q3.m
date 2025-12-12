%% Simulation_1 Q_3 Amirhossein Ehsani_810602159
clear all
close all
clc
%% Unstabel System
% a1 = 0.25;
% a2 = -0.02;
% a3 = 0.01;
% a4 = 0.02;
% 
% b1 = -0.6;
% b2 = 0.8;
% b3 = 0.5;
% b4 = 0.4;
%% Stable System 
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

Ts = 0.01;
z = tf('z',Ts);
G_d = (b1*z^3+b2*z^2+b3*z+b4)/(z^4-a1*z^3-a2*z^2-a3*z-a4);
% Proportional Controller 
K=0.2;
H = feedback(G_d,K);
step(H)

a1 = H.den{1}(2);
a2 = H.den{1}(3);
a3 = H.den{1}(4);
a4 = H.den{1}(5);

b1 = H.num{1}(2);
b2 = H.num{1}(3);
b3 = H.num{1}(4);
b4 = H.num{1}(5);
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
    
    y(t)=a1*y(t-1)+a2*y(t-2)+a3*y(t-1)+b1*u(t-1)+b2*u(t-2)+b3*u(t-4)+e(t)+c1*e(t-1);
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

%% Output KF
y_KF(1) = e(1);
y_KF(2) = Teta_kf(1,N)*y_KF(1)+Teta_kf(5,N)*u(1)+e(2);
y_KF(3) = Teta_kf(1,N)*y_KF(2)+Teta_kf(2,N)*y_KF(1)+Teta_kf(5,N)*u(2)+Teta_kf(6,N)*u(1)+e(3);
y_KF(4) = Teta_kf(1,N)*y_KF(3)+Teta_kf(2,N)*y_KF(2)+Teta_kf(3,N)*y_KF(1)+Teta_kf(5,N)*u(3)+Teta_kf(6,N)*u(2)+Teta_kf(7,N)*u(1)+e(4);
for t=5:N
    y_KF(t)=Teta_kf(1,N)*y_KF(t-1)+Teta_kf(2,N)*y_KF(t-2)+Teta_kf(3,N)*y_KF(t-3)+Teta_kf(4,N)*y_KF(t-4)+Teta_kf(5,N)*u(t-1)+Teta_kf(6,N)*u(t-2)+Teta_kf(7,N)*u(t-3)+Teta_kf(8,N)*u(t-4)+e(t)+c1*e(t-1);
end
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

Teta_RLS=q4*ones(8,N);
Teta_RLS(:,1:4)= 0;
K = zeros(8,N);
K(:,1:4)= 0;
for t=5:N
    p_inv = inv(P(:,:,t-1))+phi(:,t)*phi(:,t)';
    P(:,:,t) = inv(p_inv);
    K(:,t)=P(:,:,t)*phi(:,t);
    Teta_RLS(:,t)=Teta_RLS(:,t-1)+K(:,t)*(y(t)-phi(:,t)'*Teta_RLS(:,t-1));
end

theta0=[a1 a2 a3 a4 b1 b2 b3 b4]';
Teta_RLS(:,N);
error_RLS=norm(Teta_RLS(:,N)-theta0)/5*100;
RMSD_RLS=sqrt(sum((Teta_RLS(:,N)-theta0).^2)/5);
%% Output RLS
y_RLS(1) = e(1);
y_RLS(2) = Teta_RLS(1,N)*y_RLS(1)+Teta_RLS(5,N)*u(1)+e(2);
y_RLS(3) = Teta_RLS(1,N)*y_RLS(2)+Teta_RLS(2,N)*y_RLS(1)+Teta_RLS(5,N)*u(2)+Teta_RLS(6,N)*u(1)+e(3);
y_RLS(4) = Teta_RLS(1,N)*y_RLS(3)+Teta_RLS(2,N)*y_RLS(2)+Teta_RLS(3,N)*y_RLS(1)+Teta_RLS(5,N)*u(3)+Teta_RLS(6,N)*u(2)+Teta_RLS(7,N)*u(1)+e(4);
for t=5:N
    y_RLS(t)=Teta_RLS(1,N)*y_RLS(t-1)+Teta_RLS(2,N)*y_RLS(t-2)+Teta_RLS(3,N)*y_RLS(t-3)+Teta_RLS(4,N)*y_RLS(t-4)+Teta_RLS(5,N)*u(t-1)+Teta_RLS(6,N)*u(t-2)+Teta_RLS(7,N)*u(t-3)+Teta_RLS(8,N)*u(t-4)+e(t)+c1*e(t-1);
end

%% RLS Plot
figure(2);
plot(1:N, Teta_RLS(1,:), 'Color', colors(1,:), 'LineWidth', 1.5); hold on;
plot(1:N, Teta_RLS(2,:), 'Color', colors(2,:), 'LineWidth', 1.5);
plot(1:N, Teta_RLS(3,:), 'Color', colors(3,:), 'LineWidth', 1.5);
plot(1:N, Teta_RLS(4,:), 'Color', colors(4,:), 'LineWidth', 1.5);

plot(1:N, Teta_RLS(5,:), 'Color', colors(5,:), 'LineWidth', 1.5);
plot(1:N, Teta_RLS(6,:), 'Color', colors(6,:), 'LineWidth', 1.5);
plot(1:N, Teta_RLS(7,:), 'Color', colors(7,:), 'LineWidth', 1.5);
plot(1:N, Teta_RLS(8,:), 'Color', colors(8,:), 'LineWidth', 1.5);

plot(pt', 'b:', 'LineWidth', 1.5);
hold off;
title('Estimated Parameters from RLS', 'FontSize', 14);
legend('a1','a2','a3','a4','b1','b2','b2','b3', 'Location', 'best');
xlabel('Sample');
ylabel('Parameter Value');
grid on;



%% Results:
method.Name = {'Real System';'RLS';'KF'};
method.a1 = [a1;Teta_RLS(1,N);Teta_kf(1,N)]; 
method.a2 = [a2;Teta_RLS(2,N);Teta_kf(2,N)];
method.a3 = [a3;Teta_RLS(3,N);Teta_kf(3,N)];
method.a4 = [a4;Teta_RLS(4,N);Teta_kf(4,N)];

method.b1 = [b1;Teta_RLS(5,N);Teta_kf(5,N)];
method.b2 = [b1;Teta_RLS(6,N);Teta_kf(6,N)];
method.b3 = [b1;Teta_RLS(7,N);Teta_kf(7,N)];
method.b4 = [b1;Teta_RLS(8,N);Teta_kf(8,N)];

method.RMSD_error = [0;RMSD_RLS;RMSD_KF];
T = struct2table(method)
%% Plotting Outputs RLS and KF
tf = 0.5;
Ts=0.0001;
t = 0:Ts:tf;
N = length(t);
y(5001)=y(5000);
y_KF(5001)=y_KF(5000);
y_RLS(5001)=y_RLS(5000);

figure(3)
plot(t, y, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]) % Actual system response in blue
hold on;
plot(t, y_KF, 'LineWidth', 1, 'Color', [0.8500, 0.3250, 0.0980]) % Estimated system response in orange
hold off;

title(['Actual System Response and Over Parameter Estimated System Response by KF'], 'Interpreter', 'none')
xlabel('Time (s)')
ylabel('System Output')
legend('Actual System', 'Estimated System', 'Location', 'best')
grid on
set(gca, 'FontSize', 12)
set(gcf, 'Color', 'w') % Sets the background to white
set(gca, 'Box', 'on') % Draws a box around the plot

figure(4)
plot(t, y, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]) % Actual system response in blue
hold on;
plot(t, y_RLS, 'LineWidth', 1, 'Color', [0.8500, 0.3250, 0.0980]) % Estimated system response in orange
hold off;

title(['Actual System Response and Over Parameter Estimated System Response by RLS'], 'Interpreter', 'none')
xlabel('Time (s)')
ylabel('System Output')
legend('Actual System', 'Estimated System', 'Location', 'best')
grid on
set(gca, 'FontSize', 12)
set(gcf, 'Color', 'w') % Sets the background to white
set(gca, 'Box', 'on') % Draws a box around the plot











