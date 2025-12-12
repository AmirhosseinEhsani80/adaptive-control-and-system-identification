%% Simulation_1 Q_4_Slow_Change_Parameters Amirhossein Ehsani_810602159
clc;
clear all;
close all;

%% System
s = tf('s');
G = tf((s-4)*(s-2)*(s-3)/((s+4)*(s+5)*(s+6)*(s+7)));
Ts=0.01;
G_d = c2d(G,Ts);
step(G_d);

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

syms A B
A = a1;
B = b1;

c1 = 0.25;
sigma = 0.07;
N = 2000;
% Noise:
e = sigma * randn(1, N);
% Input Signal:
u = sigma / 0.1 * randn(1, N);
q2 = input('How many times do you want to reset the covariance error matrice(P)?\n\1)One Time \n\2)Two times \n\4)Four times \n\5)Five times\n\')
TT = N/q2

% System Response:
y(1)=e(1);
y(2)=a1*y(1)+b1*u(1)+e(2)+c1*e(2);
y(3)=a1*y(2)+a2*y(1)+b1*u(2)+b2*u(1)+e(3)+c1*e(3);
y(4)=a1*y(3)+a2*y(2)+a3*y(1)+b1*u(3)+b2*u(2)+b3*u(1)+e(4)+c1*e(4);

x = 0;
z=0;
for t=5:N    
if t>=900  && t<1100
    x = x+0.0002;
    b2=b2+x;
end
if t>=400  && t<600
    z= z+0.0001;
    b3=b3+z;
end
y(t)=a1*y(t-1)+a2*y(t-2)+a3*y(t-3)+a4*y(t-4)+b1*u(t-1)+b2*u(t-2)+b3*u(t-3)+b4*u(t-4)+e(t)+c1*e(t-1);
pt(1,t)=a1;
pt(2,t)=a2;
pt(3,t)=a3;
pt(4,t)=a4;

pt(5,t)=b1;
pt(6,t)=b2;
pt(7,t)=b3;
pt(8,t)=b4;

pt(9,t)=c0;
pt(10,t)=c1;

end 
phi = zeros(N,8);
phi(2,:) = [y(1) 0 0 0 u(1) 0 0 0];
phi(3,:) = [y(2) y(1) 0 0 u(2) u(1) 0 0];
phi(4,:) = [y(3) y(2) y(1) 0 u(3) u(2) u(1) 0];
for i=5:N
    phi(i,:) = [y(i-1) y(i-2) y(i-3) y(i-4) u(i-1) u(i-2) u(i-3) u(i-4)];
end

phi = phi';
theta0 = [a1 a2 a3 a4 b1 b2 b3 b4]
% K-P-Reset:
p_kf(:,:,1)=10^5*eye(8);
p_kf(:,:,2)=p_kf(:,:,1);
p_kf(:,:,3)=p_kf(:,:,1);
p_kf(:,:,4)=p_kf(:,:,1);

Teta_kf(:,1)=[0;0;0;0;0;0;0;0];
Teta_kf(:,2)=Teta_kf(:,1);
Teta_kf(:,3)=Teta_kf(:,1);
Teta_kf(:,4)=Teta_kf(:,1);

K_KF(:,1)=zeros(8,1);
K_KF(:,2)=zeros(8,1);
K_KF(:,3)=zeros(8,1);
K_KF(:,4)=zeros(8,1);

for t=5:N
    if mod(t,TT)==0
        p_kf(:,:,t-1)=10^5*eye(8);
    end
    K_KF(:,t)=p_kf(:,:,t-1)*phi(:,t)*inv(1+phi(:,t)'*p_kf(:,:,t-1)*phi(:,t));
    p_kf(:,:,t)=p_kf(:,:,t-1)-p_kf(:,:,t-1)*phi(:,t)*inv(1+phi(:,t)'*p_kf(:,:,t-1)*phi(:,t))*phi(:,t)'*p_kf(:,:,t-1)+0.05;
    epsilon(t)=y(t)-phi(:,t)'*Teta_kf(:,t-1);
    Teta_kf(:,t)=Teta_kf(:,t-1)+K_KF(:,t)*epsilon(t);
end
Teta_kf(:,N)
error_KF=norm(Teta_kf(:,N)-theta0)/5*100
RMSD_KF=sqrt(sum((Teta_kf(:,N)-theta0).^2)/5);
% Results in Table
method.Name = {'Real System';'KF'};
method.a1 =[a1;Teta_kf(1,N)]; 
method.a2 =[a2;Teta_kf(2,N)];
method.a3 =[a3;Teta_kf(3,N)];
method.a4 =[a4;Teta_kf(4,N)];

method.b1 =[b1;Teta_kf(5,N)];
method.b2 =[b2;Teta_kf(6,N)];
method.b3 =[b3;Teta_kf(7,N)];
method.b4 =[b4;Teta_kf(8,N)];
%method.RMSD_error = [0;RMSD_KF];
T = struct2table(method)
% Reults in Plot
% Define a color palette for better visuals
colors = lines(10); % This will generate 7 distinct colors

% KF Plot
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

figure(2);
hold on; 

% Plotting the estimated parameters
for i = 1:4
    plot(1:N, Teta_kf(i,:),'Color',colors(i,:) ,'LineWidth', 1.5, 'DisplayName', sprintf('a%d', i));
    plot(1:N,pt(i,:)','Color',colors(i,:), 'LineWidth', 1.5,'DisplayName', sprintf('True a%d', i)); % Plotting pt with a blue dotted line
end
for i = 5:8
    plot(1:N, Teta_kf(i,:), 'Color',colors(i,:),'LineWidth', 1.5, 'DisplayName', sprintf('b%d', i-4));
    plot(1:N,pt(i,:)', 'Color',colors(i,:),'LineWidth', 1.5,'DisplayName', sprintf('True b%d', i-4))
end

title('Estimated vs. True Parameters from KF');
legend('show', 'Location', 'best'); % Adjust to show all labels correctly
xlabel('Sample');
ylabel('Parameter Value');
grid on;
hold off;

