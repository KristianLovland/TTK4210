%% Initialize data
clear all; clc;
Import_Logging_lst;

close all;

%% Decide which variables to plot
M_D = true;
M_B = true;
tile_stuff = false;
savefigures = false;


%% Plot stuff
% Use colors 3-7 for process values, and 8-10 for inputs
colors = parula(10);

if M_D
    % Plot reflux drum liquid level M_D
    figure;
    subplot(2, 1, 1);
    plot(Time, LC1016(:, 1), 'color', colors(3, :));
    hold on;
    plot(Time, LC1016(:, 2), 'color', colors(5, :));
    xlabel('t', 'Interpreter', 'Latex');
    ylabel('M\_D', 'Interpreter', 'Latex');
    legend('M\_D','M\_{D, ref}', 'Interpreter', 'Latex');
    title('\textbf{Reflux drum liquid level}', 'Interpreter', 'Latex');

    subplot(2, 1, 2);
    plot(Time, LC1016(:, 3), 'color', colors(8, :));
    xlabel('t', 'Interpreter', 'Latex');
    ylabel('\textrm{Process input}', 'Interpreter', 'Latex');
    title('\textbf{Reflux drum liquid level, input}', 'Interpreter', 'Latex');
    
    if savefigures
        saveas(gcf, 'LC1016', 'eps');
    end
end

if M_B
    % Plot distillation column liquid level M_B
    figure;
    subplot(2, 1, 1);
    plot(Time, LC1015(:, 1), 'color', colors(3, :));
    hold on;
    plot(Time, LC1015(:, 2), 'color', colors(5, :));
    xlabel('t', 'Interpreter', 'Latex');
    ylabel('M\_B', 'Interpreter', 'Latex');
    legend('M\_B','M\_{B, ref}', 'Interpreter', 'Latex');
    title('\textbf{Distillation column liquid level}', 'Interpreter', 'Latex');

    subplot(2, 1, 2);
    plot(Time, LC1015(:, 3), 'color', colors(8, :));
    xlabel('t', 'Interpreter', 'Latex');
    ylabel('\textrm{Process input}', 'Interpreter', 'Latex');
    title('\textbf{Distillation column liquid level, input}', 'Interpreter', 'Latex');
    
    if savefigures
        saveas(gcf, 'LC1015', 'eps');
    end
end


if tile_stuff
    tilefigs;
end