%% Initialize data
clear all; clc;
Import_Logging_lst;

close all;

%% Decide which variables to plot
T_D = true;
T_B = true;
tile_stuff = false;
savefigures = false;


%% Plot stuff
% Use colors 3-7 for process values, and 8-10 for inputs
colors = parula(10);

if T_D
    % Plot temperature in top of distillation column, T_D
    figure;
    subplot(2, 1, 1);
    plot(Time, TC1015(:, 1), 'color', colors(3, :));
    hold on;
    plot(Time, TC1015(:, 2), 'color', colors(5, :));
    xlabel('t', 'Interpreter', 'Latex');
    ylabel('T\_D', 'Interpreter', 'Latex');
    legend('T\_D','T\_{D, ref}', 'Interpreter', 'Latex');
    title('\textbf{Distillation column temperature at top}', 'Interpreter', 'Latex');

    subplot(2, 1, 2);
    plot(Time, TC1015(:, 3), 'color', colors(8, :));
    xlabel('t', 'Interpreter', 'Latex');
    ylabel('\textrm{Process input}', 'Interpreter', 'Latex');
    title('\textbf{Distillation column temperature at top, input}', 'Interpreter', 'Latex');
    
    if savefigures
        saveas(gcf, 'TC1015', 'eps');
    end
end

if T_B
    % Plot temperature in top of distillation column, T_D
    figure;
    subplot(2, 1, 1);
    plot(Time, TC1088(:, 1), 'color', colors(3, :));
    hold on;
    plot(Time, TC1088(:, 2), 'color', colors(5, :));
    xlabel('t', 'Interpreter', 'Latex');
    ylabel('T\_B', 'Interpreter', 'Latex');
    legend('T\_B','T\_{B, ref}', 'Interpreter', 'Latex');
    title('\textbf{Distillation column temperature at top}', 'Interpreter', 'Latex');

    subplot(2, 1, 2);
    plot(Time, TC1088(:, 3), 'color', colors(8, :));
    xlabel('t', 'Interpreter', 'Latex');
    ylabel('\textrm{Process input}', 'Interpreter', 'Latex');
    title('\textbf{Distillation column temperature at top, input}', 'Interpreter', 'Latex');
    
    if savefigures
        saveas(gcf, 'TC1088', 'eps');
    end
end


if tile_stuff
    tilefigs;
end