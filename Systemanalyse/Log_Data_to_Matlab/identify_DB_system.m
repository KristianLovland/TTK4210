%% Initialize data
clear all; clc;
Import_Logging_lst;

close all;

%% Identify system
% State vector is y = [M_D; M_B], with corresponding input u = [D; B];
Y = [LC1016(:, 1) LC1015(:, 1)];
U = [FC1005(:, 1) FC1019(:, 1)];
% Dimensional limit parameter
L = 3;

[A,B,C,D,CF,F,x0]=dsr(Y,U,L);

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
Kp1 = 1;
Ti1 = 100000;
Kp2 = 1;
Ti2 = 100000;

k1 = tf([Kp1*Ti1 1], [Ti1 0]);
k2 = tf([Kp2*Ti2 1], [Ti2 0]);

%% Closed-loop stuff
l1 = g11*k1;
l2 = g22*k2;

S1 = inv(1 + l1);
S2 = inv(1 + l2);
T1 = 1 - S1;
T2 = 1 - S2;

figure;
bodemag(S1);
hold on;
bodemag(T1);
legend('S1', 'T1');
title('Magnitude');

%% Plot relevant variables
% figure;
% subplot(2, 2, 1);
% bode(RGA(1, 1));
% subplot(2, 2, 2);
% bode(RGA(1, 2));
% subplot(2, 2, 3);
% bode(RGA(2, 1));
% subplot(2, 2, 4);
% bode(RGA(2, 2));

% figure;
% margin(G(1, 1));
% 
% figure;
% margin(G(2, 2));