%% Analyse LV system
% Assume all variables are in workspace

%% Define controller
% % New controller parameters [NOT WORKING WELL]
% Kp1 = 12.5;
% Ti1 = 1000;
% % Bandwith limitation (?) gives 35 and not 135, maybe
% Kp2 = 35;
% Ti2 = 1000;

% Controller parameters used in previous experiments
Kp1 = 12;
Ti1 = 100;
Kp2 = 25;
Ti2 = 1000;

% New analysis, try to keep gain low, only require efficient control in
% frequency range below omega = 1e-2
% Kp1 = 7;
% Ti1 = 500;
% Kp2 = 15;
% Ti2 = 1000;

% Try and fail
% Kp1 = 1;
% Ti1 = 1000000000;
% Kp2 = 1;
% Ti2 = 1000000000;


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

L = G*K;
l11 = L(1, 1);
l12 = L(1, 2);
l21 = L(2, 1);
l22 = L(2, 2);

S1 = 1 / (1 + g11*k1);
S2 = 1 / (1 + g22*k2);

%% Singular value bandwith stuff
G_func = @(s) [( cell2mat(g11.num) * [s^2; s; 1] ) / ( cell2mat(g11.den) * [s^2; s; 1] ), ...
               ( cell2mat(g12.num) * [s^2; s; 1] ) / ( cell2mat(g12.den) * [s^2; s; 1] ); ...
               ( cell2mat(g21.num) * [s^2; s; 1] ) / ( cell2mat(g21.den) * [s^2; s; 1] ), ...
               ( cell2mat(g22.num) * [s^2; s; 1] ) / ( cell2mat(g22.den) * [s^2; s; 1] )];

K_func = @(s) [ Kp1*(1 + Ti1*s)/(Ti1*s)     0;
                0                           Kp2*(1 + Ti2*s)/(Ti2*s)];

% n = 200;
% singular_values = zeros(200, 2);
% omega = logspace(-5, 5, n);
% for i = 1:n
%     % S = inv(eye(2) + G_func(1j * omega(i)) * K_func (1j * omega(i)) );
%     % [U, singular_value, V] = svd(S);
%     [U, singular_value, V] = svd(G_func(1j * omega(i)));
%     singular_values(i, 1) = singular_value(1, 1);
%     singular_values(i, 2) = singular_value(2, 2);
% end

%% Plot other stuff
close all;

% Check stability
% 'Poles of loop transfer function L(s)'
% pole(L)
% nyquist((1+l11)*(1+l22) - l12*l21);
% xlim([-150 50]);
% ylim([-20 20]);
% title('Nyquist plot of L(s)');
% addPanZoom;

% figure;
% loglog(omega, singular_values(:, 1))
% hold on;
% loglog(omega, singular_values(:, 2))
% yline(0.7, '--');
% xlabel('omega');
% ylabel('Singular value');
% legend('Max singular value', 'Min singular value', 'y = 0,7');
% title('Singular values of some matrix');
% set(findall(gcf,'type','line'),'linewidth', 2);
% grid on;

% 
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
title('\textbf{Performance requirements for L\_1}', 'Interpreter', 'Latex');
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
title('\textbf{Performance requirements for L\_2}', 'Interpreter', 'Latex');
set(findall(gcf,'type','line'),'linewidth', 2);
grid on;
