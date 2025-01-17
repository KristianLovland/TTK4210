%% Initialize data
clear all; clc;
Import_Logging_lst;

close all;

%% Identify system
% State vector is y = [T_D; T_B], with corresponding input u = [L; V];
Y = [TC1015(:, 1), TC1088(:, 1)];
U = [FC1015(:, 1), LC1028(:, 1)];

% Scaling is done to keep all signals in the range [0, 1]
Y = (Y - 25*ones(size(Y))) * [1/25 0; 0 1/25];
U = U * [1/120 0; 0 1/100];

% Dimensional limit parameter
G = 3;

[A,B,C,D,CF,F,x0]=dsr(Y,U,G);

%% Open-loop stuff
sample_time = min(diff(Time));
disc_system = ss(A, B, C, D, sample_time);
cont_system = d2c(disc_system);

A = cont_system.A;
B = cont_system.B;
C = cont_system.C;
D = cont_system.D;

[G_i1_num, G_i1_den] = ss2tf(A, B, C, D, 1);
[G_i2_num, G_i2_den] = ss2tf(A, B, C, D, 2);

g11 = tf(G_i1_num(1, :), G_i1_den);
g12 = tf(G_i2_num(1, :), G_i2_den);
g21 = tf(G_i1_num(2, :), G_i1_den);
g22 = tf(G_i2_num(2, :), G_i2_den);

G = [g11 g12; g21 g22];
RGA = G .* inv(G).';

%% Controller synthesis
% Kp1 = 1;
% Ti1 = 5000;
% Kp2 = 1.7;
% Ti2 = 5000;
% 
% k1 = tf([Kp1*Ti1 1], [Ti1 0]);
% k2 = tf([Kp2*Ti2 1], [Ti2 0]);
% 
% %% Closed-loop stuff
% S1 = inv(1 + l11);
% S2 = inv(1 + l22);
% T1 = 1 - S1;
% T2 = 1 - S2;

% figure;
% bodemag(S1);
% hold on;
% bodemag(T1);
% legend('S1', 'T1');
% title('Magnitude');

%% Plot relevant variables
% figure
% subplot(2, 1, 1)
% plot(Time, Y(:, 1))
% title('T_D');
% subplot(2, 1, 2)
% plot(Time, U(:, 1))
% title('L');
% 
% figure
% subplot(2, 1, 1)
% plot(Time, Y(:, 2))
% title('T_B');
% subplot(2, 1, 2)
% plot(Time, U(:, 2))
% title('V');

figure;
bodemag(RGA);
set(findall(gcf,'type','line'),'linewidth', 2);
title('\textbf{RGA for LV system}', 'Interpreter', 'Latex');
grid on;

% figure;
% margin(g11);
% % title('\textrm{Transfer function from }$M\_{D, ref}$ to $M\_D$', 'Interpreter', 'Latex');
% grid on;
% 
% figure;
% margin(g22);
% % title('\textrm{Transfer function from }$M\_{B, ref}$ to $M\_B$', 'Interpreter', 'Latex');
% grid on;
% 
% figure;
% margin(k1*l11);
% % title('\textrm{Transfer function from }$M\_{D, ref}$ to $M\_D$ \textrm{with PI control}', 'Interpreter', 'Latex');
% grid on;
% 
% figure;
% margin(k2*l22);
% % title('\textrm{Transfer function from }$M\_{B, ref}$ to $M\_B$ \textrm{with PI control}', 'Interpreter', 'Latex');
% grid on;