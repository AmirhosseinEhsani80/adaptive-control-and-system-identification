%% Simulation_1 Q_1 Amirhossein Ehsani_810602159
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
k1=600
k2=400
k3=400
k4=133
b1=25
b2=20

m=1
R=0.2
r=0.1
J=1/2*m*R^2

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
[num, den]=ss2tf(A,B,C,D)
G=tf(num,den)

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
G_d = c2d(G,Ts)

% Discrete state space
ss_G = ss(G);
ssd_G = c2d(ss_G, Ts)

tf = 5;
t = 0:Ts:tf; t = t';
N = length(t);

% figure(1)
% step(G,t)
% hold on;
% step(G_d,'r',t)
% legend('Continious System','Discrete System ')
% set(gca,'fontsize',9)

%% Creating Input Signals
signal = input("Whats your desired signal for System Identification ? \nType the number of your signal\n1 for 'impulse'\n2 for 'step'\n3 for 'sin'\n4 for 'Sum of Two Sine'\n5 for 'Ramp'\n6 for 'white noise'\n7 for 'colored noise':  ")
if signal == 1
    u=zeros(N,1);
    u(1)=1;
    s = "Impulse"
elseif signal == 2
    u=ones(N,1)
    s = "Step"
elseif signal == 3
    u=sin(2*pi*t);
    s = "Sine"
elseif signal == 4
    u1 = sin(2*pi*t);
    u2=sin(2*pi*5*t);
    u=u1+u2;
    s = "Sine1 + Sine2"
elseif signal == 5
    u=t;
    s = "Ramp"
elseif signal == 6
    u=sqrt(1.7)*randn(N,1);
    u=u-mean(u);
    s = "White Noise"
elseif signal == 7
    u=sqrt(1.3)*randn(N,1);
    u=u-mean(u);
    s = "Colored Noise"
    for i=2:N
        u(i)=u(i)+0.35*u(i-1);
    end
end

%% System Response (Difference Equations Same Order Model) 
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
T = [a1;a2;a3;a4;a5;a6;b1;b2;b3;b4;b5;b6];
y = zeros(N,1);

% We will also creat a white noise for the system
var = 0.15;
e = sqrt(var)*randn(N,1);
e = e-mean(e);

% The first 6 signals must be created manually 
y(1) = e(1);
y(2) = a1*y(1)+b1*u(1)+e(2);
y(3) = a1*y(2)+a2*y(1)+b1*u(2)+b2*u(1)+e(3);
y(4) = a1*y(3)+a2*y(2)+a3*y(1)+b1*u(3)+b2*u(2)+b3*u(1)+e(4);
y(5) = a1*y(4)+a2*y(3)+a3*y(2)+a4*y(1)+b1*u(4)+b2*u(3)+b3*u(2)+b4*u(1)+e(5);
y(6) = a1*y(5)+a2*y(4)+a3*y(3)+a4*y(2)+a5*y(1)+b1*u(5)+b2*u(4)+b3*u(3)+b4*u(2)+b5*u(1)+e(5);

for i = 7:N
    y(i) = a1*y(i-1)+a2*y(i-2)+a3*y(i-3)+a4*y(i-4)+a5*y(i-5)+a6*y(i-6)+...
        b1*u(i-1)+b2*u(i-2)+b3*u(i-3)+b4*u(i-4)+b5*u(i-5)+b6*u(i-6);
end

%% Estimation using LS method (the same order as the system)

Phi = [0 0 0 0 0 0 0 0 0 0 0 0;
    y(1) 0 0 0 0 0 u(1) 0 0 0 0 0;
    y(2) y(1) 0 0 0 0 u(2) u(1) 0 0 0 0;
    y(3) y(2) y(1) 0 0 0 u(3) u(2) u(1) 0 0 0;
    y(4) y(3) y(2) y(1) 0 0 u(4) u(3) u(2) u(1) 0 0;
    y(5) y(4) y(3) y(2) y(1) 0 u(5) u(4) u(3) u(2) u(1) 0;
    y(6:N-1) y(5:N-2) y(4:N-3) y(3:N-4) y(2:N-5) y(1:N-6) u(6:N-1) u(5:N-2) u(4:N-3) u(3:N-4) u(2:N-5) u(1:N-6)];

T_LS = inv(Phi'*Phi)*Phi'*y;
y_LS = zeros(N,1);

% The first 6 signals must be created manually 
y_LS(1) = e(1);
y_LS(2) = T_LS(1)*y(1)+T_LS(7)*u(1)+e(2);
y_LS(3) = T_LS(1)*y(2)+T_LS(2)*y(1)+T_LS(7)*u(2)+T_LS(8)*u(1)+e(3);
y_LS(4) = T_LS(1)*y(3)+T_LS(2)*y(2)+T_LS(3)*y(1)+T_LS(7)*u(3)+T_LS(8)*u(2)+T_LS(9)*u(1)+e(4);
y_LS(5) = T_LS(1)*y(4)+T_LS(2)*y(3)+T_LS(3)*y(2)+T_LS(4)*y(1)+T_LS(7)*u(4)+T_LS(8)*u(3)+T_LS(9)*u(2)+T_LS(10)*u(1)+e(5);
y_LS(6) = T_LS(1)*y(5)+T_LS(2)*y(4)+T_LS(3)*y(3)+T_LS(4)*y(2)+T_LS(5)*y(1)+T_LS(7)*u(5)+T_LS(8)*u(4)+T_LS(9)*u(3)+T_LS(10)*u(2)+T_LS(11)*u(1)+e(5);

for i = 7:N
    y_LS(i) = T_LS(1)*y(i-1)+T_LS(2)*y(i-2)+T_LS(3)*y(i-3)+T_LS(4)*y(i-4)+T_LS(5)*y(i-5)+T_LS(6)*y(i-6)+...
        T_LS(7)*u(i-1)+T_LS(8)*u(i-2)+T_LS(9)*u(i-3)+T_LS(10)*u(i-4)+T_LS(11)*u(i-5)+T_LS(12)*u(i-6);
end

E = y-y_LS';
JJ=0.5*(E'*E);

method.Name = {'Real System'; 'Least Square'};
method.a1 = [a1; T_LS(1)];
method.a2 = [a2; T_LS(2)];
method.a3 = [a3; T_LS(3)];
method.a4 = [a4; T_LS(4)];
method.a5 = [a5; T_LS(5)];
method.a6 = [a6; T_LS(6)];


method.b1 = [b1; T_LS(7)];
method.b2 = [b2; T_LS(8)];
method.b3 = [b3; T_LS(9)];
method.b4 = [b4; T_LS(10)];
method.b5 = [b5; T_LS(11)];
method.b6 = [b6; T_LS(12)];

T = struct2table(method)
disp("The table variables belongse to: "+s)

% Calculating RMS error
rms_error = sqrt(mean(E.^2));

figure(2)
plot(t, y, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]) % Actual system response in blue
hold on;
plot(t, y_LS, 'LineWidth', 2, 'Color', [0.8500, 0.3250, 0.0980]) % Estimated system response in orange
hold off;

title(['Actual System Response and Estimated System Response by ' s], 'Interpreter', 'none')
xlabel('Time (s)')
ylabel('System Output')
legend('Actual System', 'Estimated System', 'Location', 'best')
grid on
set(gca, 'FontSize', 12)
set(gcf, 'Color', 'w') % Sets the background to white
set(gca, 'Box', 'on') % Draws a box around the plot

%% Under Parameter order 4

Phi = [0 0 0 0 0 0 0 0;
    y(1) 0 0 0 u(1) 0 0 0;
    y(2) y(1) 0 0 u(2) u(1) 0 0;
    y(3) y(2) y(1) 0 u(3) u(2) u(1) 0;
       y(4:N-1) y(3:N-2) y(2:N-3) y(1:N-4) u(4:N-1) u(3:N-2) u(2:N-3) u(1:N-4)];
T_ls = inv(Phi'*Phi)*Phi'*y;
y_LS = zeros(N,1);
y_LS(1) = e(1);
y_LS(2) = T_ls(1)*y_LS(1)+T_ls(5)*u(1)+e(2);
y_LS(3) = T_ls(1)*y_LS(2)+T_ls(2)*y_LS(1)+T_ls(5)*u(2)+T_ls(6)*u(1)+e(3);
y_LS(4) = T_ls(1)*y_LS(3)+T_ls(2)*y_LS(2)+T_ls(3)*y_LS(1)+T_ls(5)*u(3)+T_ls(6)*u(2)+T_ls(7)*u(1)+e(4);
for i = 5:N
    y_LS(i) = T_ls(1)*y_LS(i-1)+T_ls(2)*y_LS(i-2)+T_ls(3)*y_LS(i-3)+T_ls(4)*y_LS(i-4)+...
              T_ls(5)*u(i-1)+T_ls(6)*u(i-2)+T_ls(7)*u(i-3)+T_ls(8)*u(i-4)+e(i);
end

E = y-y_LS';
JJ=0.5*(E'*E);

figure(3)
plot(t, y, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]) % Actual system response in blue
hold on;
plot(t, y_LS, 'LineWidth', 2, 'Color', [0.8500, 0.3250, 0.0980]) % Estimated system response in orange
hold off;

title(['Actual System Response and Under Parameter Estimated System Response by ' s], 'Interpreter', 'none')
xlabel('Time (s)')
ylabel('System Output')
legend('Actual System', 'Estimated System', 'Location', 'best')
grid on
set(gca, 'FontSize', 12)
set(gcf, 'Color', 'w') % Sets the background to white
set(gca, 'Box', 'on') % Draws a box around the plot
clear method

method.Name = {'Real System'; 'Least Square'};
method.a1 = [a1; T_LS(1)];
method.a2 = [a2; T_LS(2)];
method.a3 = [a3; T_LS(3)];
method.a4 = [a4; T_LS(4)];

method.b1 = [b1; T_LS(7)];
method.b2 = [b2; T_LS(8)];
method.b3 = [b3; T_LS(9)];
method.b4 = [b4; T_LS(10)];

T = struct2table(method)
disp("Under param")

%% Over parameter order 7

Phi = [0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    y(1) 0 0 0 0 0 0 u(1) 0 0 0 0 0 0;
    y(2) y(1) 0 0 0 0 0 u(2) u(1) 0 0 0 0 0;
    y(3) y(2) y(1) 0 0 0 0 u(3) u(2) u(1) 0 0 0 0;
    y(4) y(3) y(2) y(1) 0 0 0 u(4) u(3) u(2) u(1) 0 0 0;
    y(5) y(4) y(3) y(2) y(1) 0 0 u(5) u(4) u(3) u(2) u(1) 0 0;
    y(6) y(5) y(4) y(3) y(2) y(1) 0 u(6) u(5) u(4) u(3) u(2) u(1) 0;
    y(7:N-1) y(6:N-2) y(5:N-3) y(4:N-4) y(3:N-5) y(2:N-6) y(1:N-7) u(7:N-1) u(6:N-2) u(5:N-3) u(4:N-4) u(3:N-5) u(2:N-6) u(1:N-7)];
T_LS = inv(Phi'*Phi)*Phi'*y;
y_LS = zeros(N,1);

y_LS(1) = e(1);
y_LS(2) = T_LS(1)*y(1)+T_LS(8)*u(1)+e(2);
y_LS(3) = T_LS(1)*y(2)+T_LS(2)*y(1)+T_LS(8)*u(2)+T_LS(9)*u(1)+e(3);
y_LS(4) = T_LS(1)*y(3)+T_LS(2)*y(2)+T_LS(3)*y(1)+T_LS(8)*u(3)+T_LS(9)*u(2)+T_LS(10)*u(1)+e(4);
y_LS(5) = T_LS(1)*y(4)+T_LS(2)*y(3)+T_LS(3)*y(2)+T_LS(4)*y(1)+T_LS(8)*u(4)+T_LS(9)*u(3)+T_LS(10)*u(2)+T_LS(11)*u(1)+e(5);
y_LS(6) = T_LS(1)*y(5)+T_LS(2)*y(4)+T_LS(3)*y(3)+T_LS(4)*y(2)+T_LS(5)*y(1)+T_LS(8)*u(5)+T_LS(9)*u(4)+T_LS(10)*u(3)+T_LS(11)*u(2)+T_LS(12)*u(1)+e(6);
y_LS(7) = T_LS(1)*y(6)+T_LS(2)*y(5)+T_LS(3)*y(4)+T_LS(4)*y(3)+T_LS(5)*y(2)+T_LS(6)*y(1)+T_LS(8)*u(6)+T_LS(9)*u(5)+T_LS(10)*u(4)+T_LS(11)*u(3)+T_LS(12)*u(2)+T_LS(13)*u(1)+e(7);

for i = 8:N
    y_LS(i) = T_LS(1)*y(i-1)+T_LS(2)*y(i-2)+T_LS(3)*y(i-3)+T_LS(4)*y(i-4)+T_LS(5)*y(i-5)+T_LS(6)*y(i-6)+T_LS(7)*y(i-7)+...
        +T_LS(8)*u(i-1)+T_LS(9)*u(i-2)+T_LS(10)*u(i-3)+T_LS(11)*u(i-4)+T_LS(12)*u(i-5)+T_LS(13)*u(i-6)+T_LS(14)*u(i-7)+e(i);

end

figure(4)
plot(t, y, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]) % Actual system response in blue
hold on;
plot(t, y_LS, 'LineWidth', 2, 'Color', [0.8500, 0.3250, 0.0980]) % Estimated system response in orange
hold off;

title(['Actual System Response and Over Parameter Estimated System Response by ' s], 'Interpreter', 'none')
xlabel('Time (s)')
ylabel('System Output')
legend('Actual System', 'Estimated System', 'Location', 'best')
grid on
set(gca, 'FontSize', 12)
set(gcf, 'Color', 'w') % Sets the background to white
set(gca, 'Box', 'on') % Draws a box around the plot

clear method

method.Name = {'Real System'; 'Least Square'};
method.a1 = [a1; T_LS(1)];
method.a2 = [a2; T_LS(2)];
method.a3 = [a3; T_LS(3)];
method.a4 = [a4; T_LS(4)];
method.a5 = [a1; T_LS(5)];
method.a6 = [a2; T_LS(6)];
method.a7 = [a3; T_LS(7)];

method.b1 = [b1; T_LS(7)];
method.b2 = [b2; T_LS(8)];
method.b3 = [b3; T_LS(9)];
method.b4 = [b4; T_LS(10)];
method.b5 = [b2; T_LS(11)];
method.b6 = [b3; T_LS(12)];
method.b7 = [b4; T_LS(13)];

T = struct2table(method)
disp("Under param")

