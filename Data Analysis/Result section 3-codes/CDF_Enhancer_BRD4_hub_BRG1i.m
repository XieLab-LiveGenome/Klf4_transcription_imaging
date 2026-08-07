clear
clc


file_path = '/Users/janaa/Desktop/MS2 transcription/Condensate datasheets/Brg1 10um condensate.xlsx';

CTRL_dist = readmatrix(file_path,'Sheet','ms intensity vs distance','Range','C2:C700');
CTRL_dist = CTRL_dist(~isnan(CTRL_dist));

CTRL_MS2 = readmatrix(file_path,'Sheet','ms intensity vs distance','Range','B2:B700');
CTRL_MS2 = CTRL_MS2(~isnan(CTRL_MS2));

thr = 4000;  % MS2 intensity
thr_d = 1000; % distance threshold

CTRL_ON_dist = CTRL_dist(CTRL_MS2>=thr);
CTRL_OFF_dist = CTRL_dist(CTRL_MS2<thr);

CTRL_ON_pos = CTRL_ON_dist(CTRL_ON_dist <= thr_d);
CTRL_ON_neg = CTRL_ON_dist(CTRL_ON_dist > thr_d);

CTRL_OFF_pos = CTRL_OFF_dist(CTRL_OFF_dist >= thr_d);
CTRL_OFF_neg = CTRL_OFF_dist(CTRL_OFF_dist < thr_d);

percent_ON_pos = (size(CTRL_ON_pos,1)/size(CTRL_dist,1))*100;
percent_ON_neg = (size(CTRL_ON_neg,1)/size(CTRL_dist,1))*100;

percent_ON_pos_norm = (size(CTRL_ON_pos,1)/size(CTRL_ON_dist,1))*100;
percent_ON_neg_norm = (size(CTRL_ON_neg,1)/size(CTRL_ON_dist,1))*100;

percent_OFF_pos = (size(CTRL_OFF_pos,1)/size(CTRL_dist,1))*100;
percent_OFF_neg = (size(CTRL_OFF_neg,1)/size(CTRL_dist,1))*100;

frac = [percent_ON_pos, percent_ON_neg, percent_OFF_pos, percent_OFF_neg] ;
frac_pos = percent_ON_pos + percent_OFF_pos;

[s1 p1]=kstest2(CTRL_ON_dist,CTRL_OFF_dist,'Alpha',0.05);
[s2 p2]=ttest2(CTRL_ON_dist,CTRL_OFF_dist,'Alpha',0.05);

% figure;
fig = figure; % Create a new figure window
set(fig, 'Color', 'w');   % Set figure background black OR WHITE
set(fig, 'Position', [100, 100, 800, 450]); % [left, bottom, width, height]

% % % %%%% ======================ON DURATION CDF================
bw=0.1;
[f, xi] = ksdensity(CTRL_ON_dist, 'function', 'cdf', 'Bandwidth', bw);  %
% plot(xi, f, 'LineWidth', 8, 'Color', 'k');
% hold on

[f2, xi2] = ksdensity(CTRL_OFF_dist, 'function', 'cdf', 'Bandwidth', bw);  % 
% plot(xi2, f2, 'LineWidth', 8, 'Color', '[0.6 0.6 0.6]');
% hold on

% 
% Define x limits
x_min1 = 0;
x_max1 = 6000;
% 
% % Extend CDFs to xlim using interpolation and extrapolation
x_query = linspace(x_min1, x_max1, 2000);  % High-res x values for smooth curves
f_interp = interp1(xi, f, x_query, 'linear', 'extrap');
f2_interp = interp1(xi2, f2, x_query, 'linear', 'extrap');
% 
% % Ensure the CDF is clipped between 0 and 1 (since extrapolation may go beyond)
f_interp = min(max(f_interp, 0), 1);
f2_interp = min(max(f2_interp, 0), 1);

plot(x_query, f_interp, 'LineWidth', 7, 'Color', '[0.9 0 0.9]');        % CTRL ON
hold on
plot(x_query, f2_interp, 'LineWidth', 7, 'Color', '[0.5 0.5 0.5]');   % ctrl off


% Labels and formatting with BLACK text
ax = gca;
ax.Color = 'w';          % Axes background WHITE
ax.XColor = 'k';         % X-axis ticks BLACK
ax.YColor = 'k';         % Y-axis ticks BLACK
ax.LineWidth = 1.5;      % Make axes lines visible
box on;                  % Add rectangular box around plot

xlabel('Enhancer- MED14 condensate 3D distance (nm)', 'FontSize', 24, 'FontWeight', 'bold', 'Color', 'k');
ylabel('Cumulative frequency', 'FontSize', 24, 'FontWeight', 'bold', 'Color', 'k');
set(gca,'FontSize', 30, 'FontWeight', 'bold');
xlim([x_min1 x_max1]);

h_legend = legend('gene ON', 'gene OFF','Location', 'southeast');
set(h_legend, ...
    'Color',      'w', ... % White background
    'TextColor',  'k', ... % Black text
    'EdgeColor',  'k');    % Black border

xticks(0:500:6000);