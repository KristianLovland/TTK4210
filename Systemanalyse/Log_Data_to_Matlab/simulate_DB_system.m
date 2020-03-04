%% Initialize model
% Assume identify_DB_system has been run

%% Initialize parameters
Kp1 = -2000;
Ti1 = 100000;
Kp2 = -1000;
Ti1 = 100000;

u1_min = 0;
u1_max = 120;
y1_min = 0;
y1_max = 120;

u2_min = 0;
u2_max = 120;
y2_min = 0;
y2_max = 120;

T_r1 = 600;
T_r2 = 800;
r1_min = 50;
r1_max = 55;
r2_min = 40;
r2_max = 45;

y0 = [50; 40];
x0 = inv(C) * y0;

%% Simulation
% Simulate for half an hour
T = 60*30;
experiment = sim('DB_identified_system_sim', T);
t = experiment.get('tout');
y = experiment.get('y');
y = y.data;
r = experiment.get('r');
r = r.data;
u = experiment.get('u');
u = u.data;

%% Plot results

figure;
subplot(2, 1, 1);
plot(y(:, 1));
hold on;
plot(r(:, 1));
xlabel('t [s]');
ylabel('Level [%]');
legend('M_D', 'M_{D, ref}');
title('Reflux drum level control');
subplot(2, 1, 2);
plot(u(:, 1));
xlabel('t [s]');
ylabel('D flow [t/h]');
legend('D');
title('Reflux drum input');

figure;
subplot(2, 1, 1);
plot(y(:, 2));
hold on;
plot(r(:, 2));
xlabel('t [s]');
ylabel('Level [%]');
legend('M_B', 'M_{B, ref}');
title('Distillation column level control');
subplot(2, 1, 2);
plot(u(:, 2));
xlabel('t [s]');
ylabel('B flow [t/h]');
legend('B');
title('Distillation column input');
