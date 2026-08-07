

clear
clc

file_path = '/Users/janaa/Desktop/MS2 transcription/Sec_tert_for CDF.xlsx';

CTRL_ON_dist_sec = readmatrix(file_path,'Sheet','COMPILED','Range','I2:I700');
CTRL_ON_dist_sec = CTRL_ON_dist_sec(~isnan(CTRL_ON_dist_sec));

CTRL_ON_dist_tert = readmatrix(file_path,'Sheet','COMPILED','Range','J2:J700');
CTRL_ON_dist_tert = CTRL_ON_dist_tert(~isnan(CTRL_ON_dist_tert));

CTRL_OFF_dist_sec = readmatrix(file_path,'Sheet','COMPILED','Range','D2:D700');
CTRL_OFF_dist_sec = CTRL_OFF_dist_sec(~isnan(CTRL_OFF_dist_sec));

CTRL_OFF_dist_tert = readmatrix(file_path,'Sheet','COMPILED','Range','E2:E700');
CTRL_OFF_dist_tert = CTRL_OFF_dist_tert(~isnan(CTRL_OFF_dist_tert));


[s1 p1]=kstest2(CTRL_ON_dist_sec,CTRL_OFF_dist_sec,'Alpha',0.05);
[s2 p2]=ttest2(CTRL_ON_dist_tert,CTRL_OFF_dist_tert,'Alpha',0.05);

% figure;
fig = figure; % Create a new figure window
set(fig, 'Color', 'w');   % Set figure background black OR WHITE
set(fig, 'Position', [100, 100, 800, 450]); % [left, bottom, width, height]

% % % %%%% ======================ON DURATION CDF================
bw=0.1;
[f, xi] = ksdensity(CTRL_ON_dist_sec, 'function', 'cdf', 'Bandwidth', bw);  %
% plot(xi, f, 'LineWidth', 8, 'Color', 'k');
% hold on

[f2, xi2] = ksdensity(CTRL_OFF_dist_sec, 'function', 'cdf', 'Bandwidth', bw);  % 
% plot(xi2, f2, 'LineWidth', 8, 'Color', '[0.6 0.6 0.6]');
% hold on

[f3, xi3] = ksdensity(CTRL_ON_dist_tert, 'function', 'cdf', 'Bandwidth', bw);  %
% plot(xi, f, 'LineWidth', 8, 'Color', 'k');
% hold on

[f4, xi4] = ksdensity(CTRL_OFF_dist_tert, 'function', 'cdf', 'Bandwidth', bw);  % 
% plot(xi2, f2, 'LineWidth', 8, 'Color', '[0.6 0.6 0.6]');
% hold on

% 
% Define x limits
x_min1 = 0;
x_max1 = 4000;
% 
% % Extend CDFs to xlim using interpolation and extrapolation
x_query = linspace(x_min1, x_max1, 2000);  % High-res x values for smooth curves
f_interp = interp1(xi, f, x_query, 'linear', 'extrap');
f2_interp = interp1(xi2, f2, x_query, 'linear', 'extrap');
f3_interp = interp1(xi3, f3, x_query, 'linear', 'extrap');
f4_interp = interp1(xi4, f4, x_query, 'linear', 'extrap');
% 
% % Ensure the CDF is clipped between 0 and 1 (since extrapolation may go beyond)
f_interp = min(max(f_interp, 0), 1);
f2_interp = min(max(f2_interp, 0), 1);
f3_interp = min(max(f3_interp, 0), 1);
f4_interp = min(max(f4_interp, 0), 1);

plot(x_query, f3_interp, 'LineWidth', 7, 'Color', '[0.9 0 0.9]');        % CTRL ON
hold on
plot(x_query, f4_interp, 'LineWidth', 7, 'Color', '[0.5 0.5 0.5]');   % ctrl off


% Labels and formatting with BLACK text
ax = gca;
ax.Color = 'w';          % Axes background WHITE
ax.XColor = 'k';         % X-axis ticks BLACK
ax.YColor = 'k';         % Y-axis ticks BLACK
ax.LineWidth = 1.5;      % Make axes lines visible
box on;                  % Add rectangular box around plot

xlabel('Enhancer- BRD4 condensate 3D distance (nm)', 'FontSize', 24, 'FontWeight', 'bold', 'Color', 'k');
ylabel('Cumulative frequency', 'FontSize', 24, 'FontWeight', 'bold', 'Color', 'k');
set(gca,'FontSize', 30, 'FontWeight', 'bold');
xlim([x_min1 x_max1]);

h_legend = legend('gene ON', 'gene OFF','Location', 'southeast');
set(h_legend, ...
    'Color',      'w', ... % White background
    'TextColor',  'k', ... % Black text
    'EdgeColor',  'k');    % Black border

xticks(0:500:4000);