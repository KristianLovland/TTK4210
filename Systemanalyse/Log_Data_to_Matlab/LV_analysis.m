%% Analyse LV system
% Assume all variables are in workspace

%% Define controller
% New controller parameters
Kp1 = 12;
Ti1 = 1000;
% Bandwith limitation (?) gives 35 and not 135, maybe
Kp2 = 35;
Ti2 = 1000;

% Controller parameters used in previous experiments
% Kp1 = 12;
% Ti1 = 100;
% Kp2 = 25;
% Ti2 = 1000;

k1 = tf([Kp1*Ti1 1], [Ti1 0]);
k2 = tf([Kp2*Ti2 1], [Ti2 0]);

K = [k1 0; 0 k2];

%% Define other variables of interest
G;
G_d = ones(2, 1);
G_tilde = [g11 0; 0 g22];
S_tilde = inv(eye(2) + G_tilde*K);
PRGA = G_tilde/G;
CLDG = PRGA*G_d;

%% Plot stuff
close all;

% figure;
% bodemag(PRGA);
% title('\textbf{PRGA for LV system}', 'Interpreter', 'Latex');
% set(findall(gcf,'type','line'),'linewidth', 2);
% grid on;

figure;
bodemag(1+g11*k1);
hold on;
bodemag(PRGA(1, 1));
bodemag(PRGA(1, 2));
% bodemag(CLDG(1, 1));
% bodemag(CLDG(1, 2));
legend('1 + L_2', 'Gamma_{11}', 'Gamma_{12}'); %, 'G~_{d11}', 'G~_{d12}');
title('Performance requirements for L_1');
set(findall(gcf,'type','line'),'linewidth', 2);
grid on;

figure;
bodemag(1+g22*k2);
hold on;
bodemag(PRGA(2, 1));
bodemag(PRGA(2, 2));
% bodemag(CLDG(2, 1));
% bodemag(CLDG(2, 2));
legend('1 + L_1', 'Gamma_{21}', 'Gamma_{22}'); %, 'G~_{d21}', 'G~_{d22}');
title('Performance requirements for L_2');
set(findall(gcf,'type','line'),'linewidth', 2);
grid on;
