%% Analyse LV system
% Assume all variables are in workspace

%% Define controller
Kp1 = 12;
Ti1 = 100;
Kp2 = 25;
Ti2 = 1000;

k1 = tf([Kp1*Ti1 1], [Ti1 0]);
k2 = tf([Kp2*Ti2 1], [Ti2 0]);

K = [k1 0; 0 k2];

%% Define other variables of interest
G;
G_d = eye(2);
G_tilde = [g11 0; 0 g22];
S_tilde = inv(eye(2) + G_tilde*K);
PRGA = G_tilde/G;
CLDG = PRGA*G_d;

%% Plot stuff
close all;

figure;
bodemag(PRGA);
title('\textbf{PRGA for LV system}', 'Interpreter', 'Latex');
set(findall(gcf,'type','line'),'linewidth', 2);
grid on;

figure;
bodemag(1+g11*k1);
hold on;
bodemag(PRGA(1, 1));
bodemag(PRGA(1, 2));
% bodemag(CLDG(1, 1));
% bodemag(CLDG(1, 2));
legend('1 + L_i', 'Gamma_{11}', 'Gamma_{12}'); %, 'G~_{d11}', 'G~_{d12}');
title('PRGA and CLDG requirements for G_{11}');
set(findall(gcf,'type','line'),'linewidth', 2);
grid on;

figure;
bodemag(1+g22*k2);
hold on;
bodemag(PRGA(2, 1));
bodemag(PRGA(2, 2));
% bodemag(CLDG(2, 1));
% bodemag(CLDG(2, 2));
legend('1 + L_i', 'Gamma_{21}', 'Gamma_{22}'); %, 'G~_{d21}', 'G~_{d22}');
title('PRGA and CLDG requirements for G_{22}');
set(findall(gcf,'type','line'),'linewidth', 2);
grid on;

figure;
bodemag(g11);
hold on;
bodemag(CLDG(1, 1));
bodemag(CLDG(1, 2));
legend('G_{11}', 'G~_{d21}', 'G~_{d22}');
title('How to avoid input constraints');
set(findall(gcf,'type','line'),'linewidth', 2);
grid on;
