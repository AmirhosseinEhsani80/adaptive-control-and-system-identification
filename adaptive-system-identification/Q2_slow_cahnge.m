%% Simulation_1 Q_2_slow_change_Params Amirhossein Ehsani_810602159


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
c1 = 0;
sigma = 0.07;
N = 1000;
%Noise:
e = sigma * randn(1, N);
% Input Signal:
u = sigma / 0.1 * randn(1, N);

q1 = input('What is the value of landa for landa-RLS?\n\')
landa = q1
q2 = input('How many times do you want to reset the covariance error matrice(P)?\n\1)One Time \n\2)Two times \n\4)Four times \n\5)Five times\n\')
TT = N/q2

%% System Response
x = 0;
z=0;
% The first 6 signals must be created manually 
y(1) = e(1);
y(2) = a1*y(1)+b1*u(1)+e(2);
y(3) = a1*y(2)+a2*y(1)+b1*u(2)+b2*u(1)+e(3);
y(4) = a1*y(3)+a2*y(2)+a3*y(1)+b1*u(3)+b2*u(2)+b3*u(1)+e(4);
y(5) = a1*y(4)+a2*y(3)+a3*y(2)+a4*y(1)+b1*u(4)+b2*u(3)+b3*u(2)+b4*u(1)+e(5);
y(6) = a1*y(5)+a2*y(4)+a3*y(3)+a4*y(2)+a5*y(1)+b1*u(5)+b2*u(4)+b3*u(3)+b4*u(2)+b5*u(1)+e(5);

for t=7:N
    if t>=500  && t<700
        x = x+0.002;
        z= z+0.001;
        b4=b4+z;
        b3=b3+x;
    end
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

%% Landa-RLS

P_1(:,:,1) = 10^5 * eye(12);
P_1(:,:,2)=P_1(:,:,1);
P_1(:,:,3)=P_1(:,:,1);
P_1(:,:,4)=P_1(:,:,1);
P_1(:,:,5)=P_1(:,:,1);
P_1(:,:,6)=P_1(:,:,1);
% Assuming landa (lambda) is known and defined somewhere in your code
Teta_rls_landa = zeros(12, N); 
Teta_rls_landa(:, 1:6) = 0;
K = zeros(12, N); 
K(:, 1:6) = 0; 
phi=Phi';

% Adjust the for loop to start from t=7 due to manual calculations
for t = 7:N
    P_1(:,:,t)=(P_1(:,:,t-1)-P_1(:,:,t-1)*phi(:,t)*inv(landa+phi(:,t)'*P_1(:,:,t-1)*phi(:,t))*phi(:,t)'*P_1(:,:,t-1))/landa;
    K(:,t)=P_1(:,:,t)*phi(:,t);
    Teta_rls_landa(:,t)=Teta_rls_landa(:,t-1)+K(:,t)*(y(t)-phi(:,t)'*Teta_rls_landa(:,t-1));
end

% Final estimated parameters with lambda-RLS
final_parameters_lambda = Teta_rls_landa(:, N);

% True parameters vector for error calculation
true_parameters = [a1; a2; a3; a4; a5; a6; b1; b2; b3; b4; b5; b6];

% Error calculations
error_RLS_lambda = norm(Teta_rls_landa(:,N) - true_parameters) / 12 * 100;
RMSD_RLS_lambda = sqrt(sum((Teta_rls_landa(:,N) - true_parameters).^2) / 12);

% Display the final error 
fprintf('Final RLS Error: %f%%\n', error_RLS_lambda);
fprintf('RMSD for RLS: %f\n', RMSD_RLS_lambda);
%% P_reset_RLS
P_2(:,:,1) = 10^5 * eye(12);
P_2(:,:,2)=P_2(:,:,1);
P_2(:,:,3)=P_2(:,:,1);
P_2(:,:,4)=P_2(:,:,1);
P_2(:,:,5)=P_2(:,:,1);
P_2(:,:,6)=P_2(:,:,1);
% Assuming landa (lambda) is known and defined somewhere in your code
Teta_rls_p = zeros(12, N); 
Teta_rls_p(:, 1:6) = 0;
K = zeros(12, N); 
K(:, 1:6) = 0; 
phi=Phi';

for t=7:N
    if mod(t,TT)==0
        P_2(:,:,t-1)=10^5*eye(12);
    end
    p_inv = inv(P_2(:,:,t-1))+phi(:,t)*phi(:,t)';
    P_2(:,:,t) = inv(p_inv);
    K(:,t)=P_2(:,:,t)*phi(:,t);
    Teta_rls_p(:,t)=Teta_rls_p(:,t-1)+K(:,t)*(y(t)-phi(:,t)'*Teta_rls_p(:,t-1));    

end

Teta_rls_p(:,N);
error_RLS_p=norm(Teta_rls_p(:,N)-true_parameters)/5*100
RMSD_RLS_p=sqrt(sum((Teta_rls_landa(:,N)-true_parameters).^2)/5)


%% Plotting Results

colors = [
    0, 0.4470, 0.7410; % MATLAB default blue
    0.8500, 0.3250, 0.0980; % MATLAB default orange
    0.9290, 0.6940, 0.1250; % MATLAB default yellow
    0.4940, 0.1840, 0.5560; % MATLAB default purple
    0.4660, 0.6740, 0.1880; % MATLAB default green
    0.3010, 0.7450, 0.9330; % MATLAB default cyan
    0.6350, 0.0780, 0.1840; % MATLAB default red
    0, 0, 0; % Black
    0.5, 0.5, 0.5; % Gray
    1, 0, 0; % Red
    0, 1, 0; % Lime
    0, 0, 1; % Blue
];

figure(1);
hold on; 

% Plotting the estimated parameters
for i = 1:6
    plot(1:N, Teta_rls_landa(i,:),'Color',colors(i,:) ,'LineWidth', 1.5, 'DisplayName', sprintf('a%d', i));
    plot(1:N,pt(i,:)','Color',colors(i,:), 'LineWidth', 1.5,'DisplayName', sprintf('True a%d', i)); % Plotting pt with a blue dotted line
end
for i = 7:12
    plot(1:N, Teta_rls_landa(i,:), 'Color',colors(i,:),'LineWidth', 1.5, 'DisplayName', sprintf('b%d', i-6));
    plot(1:N,pt(i,:)', 'Color',colors(i,:),'LineWidth', 1.5,'DisplayName', sprintf('True b%d', i))
end

title('Estimated vs. True Parameters from \lambda-RLS');
legend('show', 'Location', 'best'); % Adjust to show all labels correctly
xlabel('Sample');
ylabel('Parameter Value');
grid on;
hold off;

figure(2);
hold on; 

% Plotting the estimated parameters
for i = 1:6
    plot(1:N, Teta_rls_p(i,:),'Color',colors(i,:) ,'LineWidth', 1.5, 'DisplayName', sprintf('a%d', i));
    plot(1:N,pt(i,:)','Color',colors(i,:), 'LineWidth', 1.5,'DisplayName', sprintf('True a%d', i)); % Plotting pt with a blue dotted line
end
for i = 7:12
    plot(1:N, Teta_rls_p(i,:), 'Color',colors(i,:),'LineWidth', 1.5, 'DisplayName', sprintf('b%d', i-6));
    plot(1:N,pt(i,:)', 'Color',colors(i,:),'LineWidth', 1.5,'DisplayName', sprintf('True b%d', i-6))
end

title('Estimated vs. True Parameters from P-Reset-RLS');
legend('show', 'Location', 'best'); % Adjust to show all labels correctly
xlabel('Sample');
ylabel('Parameter Value');
grid on;
hold off;

%% Results in Table
method.Name = {'Real System';'Landa-RLS';'RLS-P-Reset'};
method.a1 =[a1;Teta_rls_landa(1,N);Teta_rls_p(1,N)]; 
method.a2 =[a2;Teta_rls_landa(2,N);Teta_rls_p(2,N)];
method.a3 =[a3;Teta_rls_landa(3,N);Teta_rls_p(3,N)];
method.a4 =[a4;Teta_rls_landa(4,N);Teta_rls_p(4,N)];
method.a5 =[a5;Teta_rls_landa(5,N);Teta_rls_p(5,N)];
method.a6 =[a6;Teta_rls_landa(6,N);Teta_rls_p(6,N)];


method.b1 =[b1;Teta_rls_landa(7,N);Teta_rls_p(7,N)];
method.b2 =[b2;Teta_rls_landa(8,N);Teta_rls_p(8,N)];
method.b3 =[b1;Teta_rls_landa(9,N);Teta_rls_p(9,N)];
method.b4 =[b2;Teta_rls_landa(10,N);Teta_rls_p(10,N)];
method.b5 =[b1;Teta_rls_landa(11,N);Teta_rls_p(11,N)];
method.b6 =[b2;Teta_rls_landa(12,N);Teta_rls_p(12,N)];

method.RMSD_error = [0;RMSD_RLS_lambda;RMSD_RLS_p];

T = struct2table(method)







