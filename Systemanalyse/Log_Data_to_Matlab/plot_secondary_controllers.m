%% Initialize data
clear all; clc;
Import_Logging_lst;

close all;

%% Decide which variables to plot
D = false;
L = false;
B = false;
V = false;
p = false;
tile_stuff = false;
addZoom = false;
savefigures = false;


%% Plot stuff
% Use colors 3-7 for process values, and 8-10 for inputs
colors = parula(10);

if D
    % Plot distillate flow D
    figure;
    subplot(2, 1, 1);
    plot(Time, FC1005(:, 1), 'color', colors(3, :));
    hold on;
    plot(Time, FC1005(:, 2), 'color', colors(5, :));
    xlim([28 32]);
    xlabel('t [s]', 'Interpreter', 'Latex');
    ylabel('D [t/h]', 'Interpreter', 'Latex');
    legend('D','D\_r', 'Interpreter', 'Latex');
    title('\textbf{Distillate flow}', 'Interpreter', 'Latex');

    subplot(2, 1, 2);
    plot(Time, FC1005(:, 3), 'color', colors(8, :));
    xlim([28 32]);
    ylim([0 100]);
    xlabel('t [s]', 'Interpreter', 'Latex');
    ylabel('\textrm{Process input} [\%]', 'Interpreter', 'Latex');
    title('\textbf{Distillate flow, input}', 'Interpreter', 'Latex');
    
    if savefigures
        saveas(gcf, 'FC1005', 'epsc');
    end
end

if L
    % Plot reflux flow L
    figure;
    subplot(2, 1, 1);
    plot(Time, FC1015(:, 1), 'color', colors(3, :));
    hold on;
    plot(Time, FC1015(:, 2), 'color', colors(5, :));
    xlim([31 35]);
    xlabel('t [s]', 'Interpreter', 'Latex');
    ylabel('L [t/h]', 'Interpreter', 'Latex');
    legend('L','L\_{ref}', 'Interpreter', 'Latex');
    title('\textbf{Reflux flow}', 'Interpreter', 'Latex');

    subplot(2, 1, 2);
    plot(Time, FC1015(:, 3), 'color', colors(8, :));
    xlim([31 35]);
    ylim([0 100]);
    xlabel('t [s]', 'Interpreter', 'Latex');
    ylabel('\textrm{Process input}', 'Interpreter', 'Latex');
    title('\textbf{Reflux flow, input}', 'Interpreter', 'Latex');
    
    if savefigures
        saveas(gcf, 'FC1015', 'epsc');
    end
end

if B
    % Plot bottom product flow B
    figure;
    subplot(2, 1, 1);
    plot(Time, FC1019(:, 1), 'color', colors(3, :));
    hold on;
    plot(Time, FC1019(:, 2), 'color', colors(5, :));
    xlim([30 34]);
    xlabel('t [s]', 'Interpreter', 'Latex');
    ylabel('B [t/h]', 'Interpreter', 'Latex');
    legend('B','B\_{ref}', 'Interpreter', 'Latex');
    title('\textbf{Distillate flow}', 'Interpreter', 'Latex');

    subplot(2, 1, 2);
    plot(Time, FC1019(:, 3), 'color', colors(8, :));
    xlim([30 34]);
    ylim([0 100]);
    xlabel('t [s]', 'Interpreter', 'Latex');
    ylabel('\textrm{Process input} [\%]', 'Interpreter', 'Latex');
    title('\textbf{Bottom product flow, input}', 'Interpreter', 'Latex');
    
    if savefigures
        saveas(gcf, 'FC1019', 'epsc');
    end
end

if V
    % Plot heat exchange area, indirectly affecting bottom vapour flowrate V
    figure;
    subplot(2, 1, 1);
    plot(Time, LC1028(:, 1), 'color', colors(3, :));
    hold on;
    %plot(Time, LC1028(:, 2), 'color', colors(5, :));
    xlim([600 1500]);
    xlabel('t [s]', 'Interpreter', 'Latex');
    ylabel('V\_{ish} [\%]', 'Interpreter', 'Latex');
    %legend('V\_{ish}','V\_{ish, ref}', 'Interpreter', 'Latex');
    title('\textbf{Heat exchange area flow}', 'Interpreter', 'Latex');

    subplot(2, 1, 2);
    plot(Time, LC1028(:, 3), 'color', colors(8, :));
    xlim([600 1500]);
    ylim([0 100]);
    xlabel('t [s]', 'Interpreter', 'Latex');
    ylabel('\textrm{Process input} [\%]', 'Interpreter', 'Latex');
    title('\textbf{Heat exchange area, input}', 'Interpreter', 'Latex');
    
    if savefigures
        saveas(gcf, 'LC1028', 'epsc');
    end
end

if p
    % Plot pressure p
    figure;
    subplot(2, 1, 1);
    plot(Time, PC1024(:, 1), 'color', colors(3, :));
    hold on;
    plot(Time, PC1024(:, 2), 'color', colors(5, :));
    xlim([250 850]);
    xlabel('t [s]', 'Interpreter', 'Latex');
    ylabel('p [\textrm{barg}]', 'Interpreter', 'Latex');
    legend('p','p\_{ref}', 'Interpreter', 'Latex');
    title('\textbf{Pressure}', 'Interpreter', 'Latex');

    subplot(2, 1, 2);
    plot(Time, PC1024(:, 3), 'color', colors(8, :));
    xlim([250 850]);
    ylim([0 100]);
    xlabel('t [s]', 'Interpreter', 'Latex');
    ylabel('\textrm{Process input} [\%]', 'Interpreter', 'Latex');
    title('\textbf{Pressure, input}', 'Interpreter', 'Latex');
    
    if savefigures
        saveas(gcf, 'PC1024', 'epsc');
    end
end


if tile_stuff
    tilefigs;
end
if addZoom
    addPanZoom;
end