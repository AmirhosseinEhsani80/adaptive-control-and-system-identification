%% Simulation_1 Q_2 Amirhossein Ehsani_810602159
% RLS SA PA LMS 

clear all
close all
clc
%%
% Defining Parameters
syms k1 k2 k3 k4 b1 b2 m R r J

%STABLE BUT VERY SLOW!

% k1=0.6
% k2=0.4
% k3=0.4
% k4=0.13
% b1=0.6
% b2=0.4

%STABLE AND GOOD!
k1=600;
k2=400;
k3=400;
k4=133;
b1=25;
b2=20;

m=1;
R=0.2;
r=0.1;
J=1/2*m*R^2;

% State Space matrices
A=[0 1 0 0 0 0;
    (-k1-k2)/m (-b1-b2)/m k2/m b2/m 0 0;
    0 0 0 1 0 0;
    k2/m b2/m (-k2-k3)/m -b2/m k3/m 0;
    0 0 0 0 0 1;
    0 0 (R^2*k4+k3*r^2)/J 0 (-R^2*k4-k3*r^2)/J 0];
B=[0;1/m;0;0;0;0];
C=[0 0 1 0 0 0];
D=[0];
[num, den]=ss2tf(A,B,C,D);
G=tf(num,den);

%% Discrete Transfer function 

% Calculate the frequency response using Bode plot data
[mag, phase, w] = bode(G);

% Convert magnitude from absolute to dB
magdB = 20*log10(squeeze(mag));

% Find the -3 dB point, which is the bandwidth of the system
refdB = magdB(1);
minus3dB = refdB - 3;
bwIndex = find(magdB <= minus3dB, 1, 'first');

% Bandwidth frequency in rad/s and converting to Hz
bwFrequencyRad = w(bwIndex);
bwFrequencyHz = bwFrequencyRad / (2*pi);
% Calculating Ts using Shanon theoreom
Ts = 1/(10*bwFrequencyHz);
disp(['Sampling time (s): ', num2str(Ts)]); 

% Discrete transfer function 
G_d = c2d(G,Ts);

% Discrete state space
ss_G = ss(G);
ssd_G = c2d(ss_G, Ts);

%% Estimation Parameters
a1 = -G_d.den{1}(2);
a2 = -G_d.den{1}(3);
a3 = -G_d.den{1}(4);
a4 = -G_d.den{1}(5);
a5 = -G_d.den{1}(6);
a6 = -G_d.den{1}(7);
b1 = G_d.num{1}(2);
b2 = G_d.num{1}(3);
b3 = G_d.num{1}(4);
b4 = G_d.num{1}(5);
b5 = G_d.num{1}(5);
b6 = G_d.num{1}(6);
c0 = 1;
c1 = 0.25;
sigma = 0.07; % variance of white noise
N = 5000;
pt = ones(14,N);

% White Noise 
e = sigma*rand(1,N);

% Creating Input Signal 

% Input Signal:
q1 = input('What is the Input Signal for the system?\n\1)Impulse \n\2)Step \n\3)Sum of Two Sin \n\4)White Noise\n\');
q2 = input('What type of noise do we have?\n\1)White Noise \n\2)Colored Noise\n3)No Noise\n');
q3 = input('What is the first Value for P?');
q4 = input('What is the first value for Theta?');
if q1 == 1
    u=zeros(N,1);
    u(1)=1;
    signal = "Impulse"
elseif q1 == 2
    u=ones(N,1)
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
elseif q2==2 
    c1 = 0.3;
elseif q2==3
    c0=0;
    c1=0;
end

%% System Response

% The first 6 signals must be created manually 
y(1) = e(1);
y(2) = a1*y(1)+b1*u(1)+e(2);
y(3) = a1*y(2)+a2*y(1)+b1*u(2)+b2*u(1)+e(3);
y(4) = a1*y(3)+a2*y(2)+a3*y(1)+b1*u(3)+b2*u(2)+b3*u(1)+e(4);
y(5) = a1*y(4)+a2*y(3)+a3*y(2)+a4*y(1)+b1*u(4)+b2*u(3)+b3*u(2)+b4*u(1)+e(5);
y(6) = a1*y(5)+a2*y(4)+a3*y(3)+a4*y(2)+a5*y(1)+b1*u(5)+b2*u(4)+b3*u(3)+b4*u(2)+b5*u(1)+c0*e(6)+c1*e(5);

for t=7:N
    pt(1,t) = a1;
    pt(2,t) = a2;
    pt(3,t) = a3;
    pt(4,t) = a4;
    pt(5,t) = a5;
    pt(6,t) = a6;
    pt(7,t) = b1;
    pt(8,t) = b2;
    pt(9,t) = b3;
    pt(10,t) = b4;
    pt(11,t) = b5;
    pt(12,t) = b6;
    pt(13,t) = c0;
    pt(14,t) = c1;
    y(t) = a1*y(t-1)+a2*y(t-2)+a3*y(t-3)+a4*y(t-4)+a5*y(t-5)+a6*y(t-6)+...
        b1*u(t-1)+b2*u(t-2)+b3*u(t-3)+b4*u(t-4)+b5*u(t-5)+b6*u(t-6)+c0*e(t)+c1*e(t-1);
end
%% RLS Method 

% Phi Matrix

Phi = zeros(N,12);
Phi(1,:)=[0 0 0 0 0 0 0 0 0 0 0 0];
Phi(2,:)=[y(1) 0 0 0 0 0 u(1) 0 0 0 0 0];
Phi(3,:)=[y(2) y(1) 0 0 0 0 u(2) u(1) 0 0 0 0];
Phi(4,:)=[y(3) y(2) y(1) 0 0 0 u(3) u(2) u(1) 0 0 0];
Phi(5,:)=[y(4) y(3) y(2) y(1) 0 0 u(4) u(3) u(2) u(1) 0 0];
Phi(6,:)=[y(5) y(4) y(3) y(2) y(1) 0 u(5) u(4) u(3) u(2) u(1) 0];

for i=7:N
    Phi(i,:)=[y(i-1) y(i-2) y(i-3) y(i-4) y(i-5) y(i-6) u(i-1) u(i-2) u(i-3) u(i-4) u(i-5) u(i-6)];
end
phi = Phi';
% Initialize P for the first iteration with dimension 12x12 for the 12 parameters

P(:,:,1) = q3 * eye(12); 
P(:,:,2)=P(:,:,1);
P(:,:,3)=P(:,:,1);
P(:,:,4)=P(:,:,1);
P(:,:,5)=P(:,:,1);
P(:,:,6)=P(:,:,1);

T_RLS = q4 * ones(12, N); 
T_RLS(:,1:6) = 0; % The first 6 values are set to zero as they are manually calculated
K = zeros(12, N); % Gain vector initialization
K(:,1:6) = 0; % First 6 values are set to zero

for t = 7:N
    p_inv = inv(P(:,:,(t-1))) + Phi(t,:)' * Phi(t,:);
    P(:,:,t) = inv(p_inv);
    K(:,t) = P(:,:,t) * Phi(t,:)';
    T_RLS(:,t) = T_RLS(:,t-1) + K(:,t) * (y(t) - Phi(t,:) * T_RLS(:,t-1));
end

true_parameters = [a1; a2; a3; a4; a5; a6; b1; b2; b3; b4; b5; b6];

% Error calculations
error_RLS = norm(T_RLS(:,N) - true_parameters) / 12 * 100;
RMSD_RLS = sqrt(sum((T_RLS(:,N) - true_parameters).^2) / 12);

% Display the final error 
fprintf('Final RLS Error: %f%%\n', error_RLS);
fprintf('RMSD for RLS: %f\n', RMSD_RLS);

% Final estimated parameters and true parameters
final_parameters_RLS = T_RLS(:,N);

method.Name = {'Real System'; 'Recurisive Least Square'};
method.a1 = [a1; final_parameters_RLS(1)];
method.a2 = [a2; final_parameters_RLS(2)];
method.a3 = [a3; final_parameters_RLS(3)];
method.a4 = [a4; final_parameters_RLS(4)];
method.a5 = [a5; final_parameters_RLS(5)];
method.a6 = [a6; final_parameters_RLS(6)];


method.b1 = [b1; final_parameters_RLS(7)];
method.b2 = [b2; final_parameters_RLS(8)];
method.b3 = [b3; final_parameters_RLS(9)];
method.b4 = [b4; final_parameters_RLS(10)];
method.b5 = [b5; final_parameters_RLS(11)];
method.b6 = [b6; final_parameters_RLS(12)];

T_RLS_T = struct2table(method)

% Plotting Results for RLS
colors = lines(14); % This will generate 7 distinct colors
% RLS Plot
figure(1);
plot(1:N, T_RLS(1,:), 'Color', colors(1,:), 'LineWidth', 1.5); hold on;
plot(1:N, T_RLS(2,:), 'Color', colors(2,:), 'LineWidth', 1.5);
plot(1:N, T_RLS(3,:), 'Color', colors(3,:), 'LineWidth', 1.5);
plot(1:N, T_RLS(4,:), 'Color', colors(4,:), 'LineWidth', 1.5);
plot(1:N, T_RLS(5,:), 'Color', colors(5,:), 'LineWidth', 1.5);
plot(1:N, T_RLS(6,:), 'Color', colors(6,:), 'LineWidth', 1.5);
plot(1:N, T_RLS(7,:), 'Color', colors(7,:), 'LineWidth', 1.5);
plot(1:N, T_RLS(8,:), 'Color', colors(8,:), 'LineWidth', 1.5);
plot(1:N, T_RLS(9,:), 'Color', colors(9,:), 'LineWidth', 1.5);
plot(1:N, T_RLS(10,:), 'Color', colors(10,:), 'LineWidth', 1.5);
plot(1:N, T_RLS(11,:), 'Color', colors(11,:), 'LineWidth', 1.5);
plot(1:N, T_RLS(12,:), 'Color', colors(12,:), 'LineWidth', 1.5);
plot(pt', 'b:', 'LineWidth', 1.5);
hold off;
title('Estimated Parameters from RLS', 'FontSize', 14);
legend('a1','a2','a3','a4','a5','a6','b1','b2','b3','b4','b5','b6','real parameters', 'Location', 'best');
xlabel('Sample');
ylabel('Parameter Value');
grid on;

%% SA Method

T_SA = q4*ones(12,N);
T_SA(:,1:6)=-2.8;
for t=7:N
    P_sa(6)= phi(:,6)'* phi(:,6)+phi(:,5)'* phi(:,5)+phi(:,4)'* phi(:,4);
    P_sa(t) = P_sa(t-1)+phi(:,t)'*phi(:,t);
    epsilon_sa(t)=(y(t)-phi(:,t)'*T_SA(:,t-1));
    T_SA(:,t)=T_SA(:,t-1)+(P_sa(t)^-1)*phi(:,t)*epsilon_sa(t);
end

n=1:N;
T_SA(:,N);
error_sa=norm(T_SA(:,N)-true_parameters)/5*100
RMSD_SA=sqrt(sum((T_SA(:,N)-true_parameters).^2)/5)

% SA Plot
figure(2);
plot(1:N, T_SA(1,:), 'Color', colors(1,:), 'LineWidth', 1.5); hold on;
plot(1:N, T_SA(2,:), 'Color', colors(2,:), 'LineWidth', 1.5);
plot(1:N, T_SA(3,:), 'Color', colors(3,:), 'LineWidth', 1.5);
plot(1:N, T_SA(4,:), 'Color', colors(4,:), 'LineWidth', 1.5);
plot(1:N, T_SA(5,:), 'Color', colors(5,:), 'LineWidth', 1.5);
plot(1:N, T_SA(6,:), 'Color', colors(6,:), 'LineWidth', 1.5);
plot(1:N, T_SA(7,:), 'Color', colors(7,:), 'LineWidth', 1.5);
plot(1:N, T_SA(8,:), 'Color', colors(8,:), 'LineWidth', 1.5);
plot(1:N, T_SA(9,:), 'Color', colors(9,:), 'LineWidth', 1.5);
plot(1:N, T_SA(10,:), 'Color', colors(10,:), 'LineWidth', 1.5);
plot(1:N, T_SA(11,:), 'Color', colors(11,:), 'LineWidth', 1.5);
plot(1:N, T_SA(12,:), 'Color', colors(12,:), 'LineWidth', 1.5);
plot(pt', 'b:', 'LineWidth', 1.5);
hold off;
title('Estimated Parameters from SA', 'FontSize', 14);
legend('a1','a2','a3','a4','a5','a6','b1','b2','b3','b4','b5','b6','real parameters', 'Location', 'best');
xlabel('Sample');
ylabel('Parameter Value');
grid on;


% Final estimated parameters and true parameters
final_parameters_SA = T_SA(:,N);

method.Name = {'Real System'; 'SA'};
method.a1 = [a1; final_parameters_SA(1)];
method.a2 = [a2; final_parameters_SA(2)];
method.a3 = [a3; final_parameters_SA(3)];
method.a4 = [a4; final_parameters_SA(4)];
method.a5 = [a5; final_parameters_SA(5)];
method.a6 = [a6; final_parameters_SA(6)];


method.b1 = [b1; final_parameters_SA(7)];
method.b2 = [b2; final_parameters_SA(8)];
method.b3 = [b3; final_parameters_SA(9)];
method.b4 = [b4; final_parameters_SA(10)];
method.b5 = [b5; final_parameters_SA(11)];
method.b6 = [b6; final_parameters_SA(12)];

T_SA_T = struct2table(method)

%% PA Method

landa = 1.2;
alpha_PA = 0.001;
T_PA = q4*ones(12,N);
T_PA(:,1) = -2.8;
Alpha(1) = 0.1;
for t=2:N
    epsilon_PA(t)=(y(t)-phi(:,t)'*T_PA(:,t-1));
    K_PA(t)= 1/(phi(:,t)'*phi(:,t)+alpha_PA)*landa;
    T_PA(:,t)=T_PA(:,t-1)+K_PA(t)*phi(:,t)*epsilon_PA(t);
end
n=1:N;
T_PA(:,N);
error_PA=norm(T_PA(:,N)-true_parameters)/5*100;
RMSD_PA=sqrt(sum((T_PA(:,N)-true_parameters).^2)/5);

% PA Plot
figure(3);
plot(1:N, T_PA(1,:), 'Color', colors(1,:), 'LineWidth', 1.5); hold on;
plot(1:N, T_PA(2,:), 'Color', colors(2,:), 'LineWidth', 1.5);
plot(1:N, T_PA(3,:), 'Color', colors(3,:), 'LineWidth', 1.5);
plot(1:N, T_PA(4,:), 'Color', colors(4,:), 'LineWidth', 1.5);
plot(1:N, T_PA(5,:), 'Color', colors(5,:), 'LineWidth', 1.5);
plot(1:N, T_PA(6,:), 'Color', colors(6,:), 'LineWidth', 1.5);
plot(1:N, T_PA(7,:), 'Color', colors(7,:), 'LineWidth', 1.5);
plot(1:N, T_PA(8,:), 'Color', colors(8,:), 'LineWidth', 1.5);
plot(1:N, T_PA(9,:), 'Color', colors(9,:), 'LineWidth', 1.5);
plot(1:N, T_PA(10,:), 'Color', colors(10,:), 'LineWidth', 1.5);
plot(1:N, T_PA(11,:), 'Color', colors(11,:), 'LineWidth', 1.5);
plot(1:N, T_PA(12,:), 'Color', colors(12,:), 'LineWidth', 1.5);
plot(pt', 'b:', 'LineWidth', 1.5);
hold off;
title('Estimated Parameters from PA', 'FontSize', 14);
legend('a1','a2','a3','a4','a5','a6','b1','b2','b3','b4','b5','b6','real parameters', 'Location', 'best');
xlabel('Sample');
ylabel('Parameter Value');
grid on;


% Final estimated parameters and true parameters
final_parameters_PA = T_PA(:,N);

method.Name = {'Real System'; 'PA'};
method.a1 = [a1; final_parameters_PA(1)];
method.a2 = [a2; final_parameters_PA(2)];
method.a3 = [a3; final_parameters_PA(3)];
method.a4 = [a4; final_parameters_PA(4)];
method.a5 = [a5; final_parameters_PA(5)];
method.a6 = [a6; final_parameters_PA(6)];


method.b1 = [b1; final_parameters_PA(7)];
method.b2 = [b2; final_parameters_PA(8)];
method.b3 = [b3; final_parameters_PA(9)];
method.b4 = [b4; final_parameters_PA(10)];
method.b5 = [b5; final_parameters_PA(11)];
method.b6 = [b6; final_parameters_PA(12)];

T_PA_T = struct2table(method)

%% LMS Method

T_LMS = q4*ones(12,N);
T_LMS(:,1) = -2.8;
Alpha(1) = 0.1;
for t=2:N
    epsilon_LMS(t)=(y(t)-phi(:,t)'*T_LMS(:,t-1));
    T_LMS(:,t)=T_LMS(:,t-1)+Alpha(t-1)*phi(:,t)*epsilon_LMS(t);
    Alpha(t) = Alpha(t-1)*0.75;
end
n=1:N;
T_LMS(:,N);
error_LMS=norm(T_LMS(:,N)-true_parameters)/5*100;
RMSD_LMS=sqrt(sum((T_LMS(:,N)-true_parameters).^2)/5);

% LMS Plot
figure(4);
plot(1:N, T_LMS(1,:), 'Color', colors(1,:), 'LineWidth', 1.5); hold on;
plot(1:N, T_LMS(2,:), 'Color', colors(2,:), 'LineWidth', 1.5);
plot(1:N, T_LMS(3,:), 'Color', colors(3,:), 'LineWidth', 1.5);
plot(1:N, T_LMS(4,:), 'Color', colors(4,:), 'LineWidth', 1.5);
plot(1:N, T_LMS(5,:), 'Color', colors(5,:), 'LineWidth', 1.5);
plot(1:N, T_LMS(6,:), 'Color', colors(6,:), 'LineWidth', 1.5);
plot(1:N, T_LMS(7,:), 'Color', colors(7,:), 'LineWidth', 1.5);
plot(1:N, T_LMS(8,:), 'Color', colors(8,:), 'LineWidth', 1.5);
plot(1:N, T_LMS(9,:), 'Color', colors(9,:), 'LineWidth', 1.5);
plot(1:N, T_LMS(10,:), 'Color', colors(10,:), 'LineWidth', 1.5);
plot(1:N, T_LMS(11,:), 'Color', colors(11,:), 'LineWidth', 1.5);
plot(1:N, T_LMS(12,:), 'Color', colors(12,:), 'LineWidth', 1.5);
plot(pt', 'b:', 'LineWidth', 1.5);
hold off;
title('Estimated Parameters from LMS', 'FontSize', 14);
legend('a1','a2','a3','a4','a5','a6','b1','b2','b3','b4','b5','b6','real parameters', 'Location', 'best');
xlabel('Sample');
ylabel('Parameter Value');
grid on;


% Final estimated parameters and true parameters
final_parameters_LMS = T_LMS(:,N);

method.Name = {'Real System'; 'LMS'};
method.a1 = [a1; final_parameters_LMS(1)];
method.a2 = [a2; final_parameters_LMS(2)];
method.a3 = [a3; final_parameters_LMS(3)];
method.a4 = [a4; final_parameters_LMS(4)];
method.a5 = [a5; final_parameters_LMS(5)];
method.a6 = [a6; final_parameters_LMS(6)];


method.b1 = [b1; final_parameters_LMS(7)];
method.b2 = [b2; final_parameters_LMS(8)];
method.b3 = [b3; final_parameters_LMS(9)];
method.b4 = [b4; final_parameters_LMS(10)];
method.b5 = [b5; final_parameters_LMS(11)];
method.b6 = [b6; final_parameters_LMS(12)];

T_LMS_T = struct2table(method)

%% ELS Estimating At Noise Presence
%e=e';
phi_ELS = zeros(N,13);
phi_ELS(1,:)=[0 0 0 0 0 0 0 0 0 0 0 0 0];
phi_ELS(2,:)=[y(1) 0 0 0 0 0 u(1) 0 0 0 0 0 0.25*e(1)];
phi_ELS(3,:)=[y(2) y(1) 0 0 0 0 u(2) u(1) 0 0 0 0 0.25*e(2)];
phi_ELS(4,:)=[y(3) y(2) y(1) 0 0 0 u(3) u(2) u(1) 0 0 0 0.25*e(3)];
phi_ELS(5,:)=[y(4) y(3) y(2) y(1) 0 0 u(4) u(3) u(2) u(1) 0 0 0.25*e(4)];
phi_ELS(6,:)=[y(5) y(4) y(3) y(2) y(1) 0 u(5) u(4) u(3) u(2) u(1) 0 0.25*e(5)];

for i=7:N
    phi_ELS(i,:)=[y(i-1) y(i-2) y(i-3) y(i-4) y(i-5) y(i-6) u(i-1) u(i-2) u(i-3) u(i-4) u(i-5) u(i-6) 0.25*e(i-1)];
end

y=y';
T_ELS = inv(phi_ELS'*phi_ELS)*phi_ELS'*y;
E_ELS = y-phi_ELS*T_ELS;

Ts = 0.04;
tf = 200;
t = 0:Ts:tf; t = t';
y(5001)=y(5000);
% The first 6 signals must be created manually 
y_ELS(1) = e(1);
y_ELS(2) = T_ELS(1)*y(1)+T_ELS(7)*u(1)+e(2);
y_ELS(3) = T_ELS(1)*y(2)+T_ELS(2)*y(1)+T_ELS(7)*u(2)+T_ELS(8)*u(1)+e(3);
y_ELS(4) = T_ELS(1)*y(3)+T_ELS(2)*y(2)+T_ELS(3)*y(1)+T_ELS(7)*u(3)+T_ELS(8)*u(2)+T_ELS(9)*u(1)+e(4);
y_ELS(5) = T_ELS(1)*y(4)+T_ELS(2)*y(3)+T_ELS(3)*y(2)+T_ELS(4)*y(1)+T_ELS(7)*u(4)+T_ELS(8)*u(3)+T_ELS(9)*u(2)+T_ELS(10)*u(1)+e(5);
y_ELS(6) = T_ELS(1)*y(5)+T_ELS(2)*y(4)+T_ELS(3)*y(3)+T_ELS(4)*y(2)+T_ELS(5)*y(1)+T_ELS(7)*u(5)+T_ELS(8)*u(4)+T_ELS(9)*u(3)+T_ELS(10)*u(2)+T_ELS(11)*u(1)+e(5);

for i = 7:N
    y_ELS(i) = T_ELS(1)*y(i-1)+T_ELS(2)*y(i-2)+T_ELS(3)*y(i-3)+T_ELS(4)*y(i-4)+T_ELS(5)*y(i-5)+T_ELS(6)*y(i-6)+...
        T_ELS(7)*u(i-1)+T_ELS(8)*u(i-2)+T_ELS(9)*u(i-3)+T_ELS(10)*u(i-4)+T_ELS(11)*u(i-5)+T_ELS(12)*u(i-6);
end
y_ELS(5001) = y_ELS(5000);

% ELS Plot
figure(5);
plot(t, y, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]) % Actual system response in blue
hold on;
plot(t, y_ELS, 'LineWidth', 2, 'Color', [0.8500, 0.3250, 0.0980]) % Estimated system response in orange
hold off;

title(['Actual System Response and Estimated System Response Using ELS method' ], 'Interpreter', 'none')
xlabel('Time (s)')
ylabel('System Output')
legend('Actual System', 'Estimated System', 'Location', 'best')
grid on
set(gca, 'FontSize', 12)
set(gcf, 'Color', 'w') % Sets the background to white
set(gca, 'Box', 'on') % Draws a box around the plot




% Final estimated parameters and true parameters

method.Name = {'Real System'; 'ELS'};
method.a1 = [a1; T_ELS(1)];
method.a2 = [a2; T_ELS(2)];
method.a3 = [a3; T_ELS(3)];
method.a4 = [a4; T_ELS(4)];
method.a5 = [a5; T_ELS(5)];
method.a6 = [a6; T_ELS(6)];


method.b1 = [b1; T_ELS(7)];
method.b2 = [b2; T_ELS(8)];
method.b3 = [b3; T_ELS(9)];
method.b4 = [b4; T_ELS(10)];
method.b5 = [b5; T_ELS(11)];
method.b6 = [b6; T_ELS(12)];
method.c1 = [c1; T_ELS(13)];

T_ELS_T = struct2table(method)


%% RESULTS FOR COMPARING

% Create a new figure for subplots
figure(6); % Using figure number 6 to avoid overwriting existing figures

% Subplot for RLS method
subplot(2,2,1); % This creates a 2x2 grid, and places the following plot in the first position
plot(1:N, T_RLS, 'LineWidth', 1.5);
hold on
plot(pt', 'b:', 'LineWidth', 1.5);
title('RLS Method');
legend('a1','a2','a3','a4','a5','a6','b1','b2','b3','b4','b5','b6');
xlabel('Sample');
ylabel('Parameter Value');
grid on;

% Subplot for SA method
subplot(2,2,2);
plot(1:N, T_SA, 'LineWidth', 1.5);
hold on
plot(pt', 'b:', 'LineWidth', 1.5);
title('SA Method');
legend('a1','a2','a3','a4','a5','a6','b1','b2','b3','b4','b5','b6');
xlabel('Sample');
ylabel('Parameter Value');
grid on;

% Subplot for PA method
subplot(2,2,3);
plot(1:N, T_PA, 'LineWidth', 1.5);
hold on
plot(pt', 'b:', 'LineWidth', 1.5);
title('PA Method');
legend('a1','a2','a3','a4','a5','a6','b1','b2','b3','b4','b5','b6');
xlabel('Sample');
ylabel('Parameter Value');
grid on;

% Subplot for LMS method
subplot(2,2,4);
plot(1:N, T_LMS, 'LineWidth', 1.5);
hold on
plot(pt', 'b:', 'LineWidth', 1.5);
title('LMS Method');
legend('a1','a2','a3','a4','a5','a6','b1','b2','b3','b4','b5','b6');
xlabel('Sample');
ylabel('Parameter Value');
grid on;

% Adjust layout to make sure labels and titles are not overlapping
sgtitle('Comparison of Parameter Estimation Methods');


% Assuming true_parameters contains the real values for a1 to a6 and b1 to b6
true_parameters = [a1; a2; a3; a4; a5; a6; b1; b2; b3; b4; b5; b6];

% Gather all final parameters and RMSD errors
all_parameters = [true_parameters, final_parameters_RLS, final_parameters_SA, final_parameters_PA, final_parameters_LMS];
all_RMSD = [NaN, RMSD_RLS, RMSD_SA, RMSD_PA, RMSD_LMS]; % No RMSD for true parameters, so use NaN

% Create a table with parameter estimates and RMSD errors for each method
method_names = {'True Parameters', 'RLS', 'SA', 'PA', 'LMS'};
method_table = array2table(all_parameters, 'VariableNames', method_names, ...
    'RowNames', {'a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'b1', 'b2', 'b3', 'b4', 'b5', 'b6'});
RMSD_table = array2table(all_RMSD, 'VariableNames', method_names, 'RowNames', {'RMSD'});

% Combine both tables vertically
combined_table = [method_table; RMSD_table];

% Display the combined table
disp(combined_table);
