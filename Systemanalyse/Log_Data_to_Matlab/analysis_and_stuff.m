%% Hey and welcome to this computer program
% I assume that a system G or g11, g22 or something is identified

close all;

Kp1 = 830;
Ti1 = 1000;
k1 = tf([Kp1*Ti1 1], [Ti1 0]);
l1 = g11*k1;
S1 = inv(1 + l1);
T1 = 1 - S1;
figure;
bodemag(S1);
title('S\_1','Interpreter', 'Latex');
grid on;

Kp2 = 1000;
Ti2 = 1000;
k2 = tf([Kp2*Ti2 1], [Ti2 0]);
l2 = g22*k2;
S2 = inv(1 + l2);
T2 = 1 - S2;
figure;
bodemag(S2);
title('S\_2','Interpreter', 'Latex');
grid on;

% figure;
% hold on;
% for k = [1, 10, 100, 500, 1000]
%     Kp1 = k;
%     Ti1 = 10000;
%     k1 = tf([Kp1*Ti1 1], [Ti1 0]);
%     
%     l1 = g11*k1;
%     S1 = inv(1 + l1);
%     T1 = 1 - S1;
%     bodemag(S1);
% end
% legend('Kp = 1', 'Kp = 10', 'Kp = 100', 'Kp = 500', 'Kp = 1000', 'Interpreter', 'Latex');
% title('S for different gains, T\_i = 100000', 'Interpreter', 'Latex');
% 
% figure;
% hold on;
% for k = [1, 10, 100, 1000]
%     Kp1 = 1;
%     Ti1 = k;
%     k1 = tf([Kp1*Ti1 1], [Ti1 0]);
%     
%     l1 = g11*k1;
%     S1 = inv(1 + l1);
%     T1 = 1 - S1;
%     bodemag(S1);
% end
% legend('Ti = 1', 'Ti = 10', 'Ti = 100', 'Ti = 1000');
% title('S for different integral times, Kp = 1');