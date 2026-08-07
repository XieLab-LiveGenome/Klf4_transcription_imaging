clear
clc

% output_file = 'Rad21 compiled';
file_path = '/Users/janaa/Desktop/MS2 transcription/4 color compiled/4 COLOR_6-29-2026.xlsx';
file_path2 = '/Users/janaa/Desktop/MS2 transcription/4 color compiled/4 COLOR R21 depl.xlsx';


CTRL_dist_e = readmatrix(file_path,'Sheet','ms intensity vs distance','Range','e2:e178');
CTRL_dist_e = CTRL_dist_e(~isnan(CTRL_dist_e));

CTRL_dist_p = readmatrix(file_path,'Sheet','ms intensity vs distance','Range','k2:k300');
CTRL_dist_p = CTRL_dist_p(~isnan(CTRL_dist_p));

CTRL_MS2 = readmatrix(file_path,'Sheet','ms intensity vs distance','Range','d2:d178');
CTRL_MS2 = CTRL_MS2(~isnan(CTRL_MS2));


R21_dist_e = readmatrix(file_path2,'Sheet','ms intensity vs distance','Range','e2:e300');
R21_dist_e = R21_dist_e(~isnan(R21_dist_e));

R21_dist_p = readmatrix(file_path2,'Sheet','ms intensity vs distance','Range','k2:k300');
R21_dist_p = R21_dist_p(~isnan(R21_dist_p));

R21_MS2 = readmatrix(file_path2,'Sheet','ms intensity vs distance','Range','d2:d300');
R21_MS2 = R21_MS2(~isnan(R21_MS2));


thr = 0.2*65536;  % MS2 intensity
% thr = 10000;
thr_d = 1000; % distance threshold

CTRL_ON_dist_e = CTRL_dist_e(CTRL_MS2>=thr);
CTRL_OFF_dist_e = CTRL_dist_e(CTRL_MS2<thr);

% CTRL_ON_dist_p = CTRL_dist_p(CTRL_MS2>=thr);
% CTRL_OFF_dist_p = CTRL_dist_p(CTRL_MS2<thr);


R21_ON_dist_e = R21_dist_e(R21_MS2>=thr);
R21_OFF_dist_e = R21_dist_e(R21_MS2<thr);

R21_ON_dist_p = R21_dist_p(R21_MS2>=thr);
R21_OFF_dist_p = R21_dist_p(R21_MS2<thr);



[s1 p1]=kstest2(CTRL_dist_e,R21_dist_e,'Alpha',0.05);
[s2 p2]=kstest2(CTRL_dist_e,R21_dist_e,'Alpha',0.05);

[s3 p3]=kstest2(CTRL_ON_dist_e,R21_ON_dist_e,'Alpha',0.05);
[s4 p4]=kstest2(CTRL_OFF_dist_e,R21_OFF_dist_e,'Alpha',0.05);

% [s5 p5]=kstest2(CTRL_ON_dist_p,R21_ON_dist_p,'Alpha',0.05);
% [s6 p6]=kstest2(CTRL_OFF_dist_p,R21_OFF_dist_p,'Alpha',0.05);
% 
% [s7 p7]=kstest2(CTRL_ON_dist_e,CTRL_OFF_dist_e,'Alpha',0.05);
% [s8 p8]=kstest2(CTRL_ON_dist_p,CTRL_OFF_dist_p,'Alpha',0.05);
% 
% figure;
fig = figure; % Create a new figure window
set(fig, 'Color', 'w');   % Set figure background black OR WHITE
set(fig, 'Position', [100, 100, 800, 450]); % [left, bottom, width, height]

% % % %%%% ======================ON DURATION CDF================

bw=0.1;
[f, xi] = ksdensity(CTRL_dist_e, 'function', 'cdf', 'Bandwidth', bw); 
[f2, xi2] = ksdensity(R21_dist_e, 'function', 'cdf', 'Bandwidth', bw);  % 

[f3, xi3] = ksdensity(CTRL_ON_dist_e, 'function', 'cdf', 'Bandwidth', bw); 
[f4, xi4] = ksdensity(R21_ON_dist_e, 'function', 'cdf', 'Bandwidth', bw);  % 

[f5, xi5] = ksdensity(CTRL_OFF_dist_e, 'function', 'cdf', 'Bandwidth', bw); 
[f6, xi6] = ksdensity(R21_OFF_dist_e, 'function', 'cdf', 'Bandwidth', bw);  % 

[f7, xi7] = ksdensity(CTRL_ON_dist_e, 'function', 'cdf', 'Bandwidth', bw); 
[f8, xi8] = ksdensity(CTRL_OFF_dist_e, 'function', 'cdf', 'Bandwidth', bw);  % 

[f9, xi9] = ksdensity(CTRL_ON_dist_p, 'function', 'cdf', 'Bandwidth', bw); 
[f10, xi10] = ksdensity(CTRL_OFF_dist_p, 'function', 'cdf', 'Bandwidth', bw);  % 
% 
% Define x limits
x_min1 = 0;
x_max1 = 3000;
% 
% % Extend CDFs to xlim using interpolation and extrapolation
x_query = linspace(x_min1, x_max1, 2000);  % High-res x values for smooth curves
f_interp = interp1(xi, f, x_query, 'linear', 'extrap');
f2_interp = interp1(xi2, f2, x_query, 'linear', 'extrap');
f3_interp = interp1(xi3, f3, x_query, 'linear', 'extrap');
f4_interp = interp1(xi4, f4, x_query, 'linear', 'extrap');
f5_interp = interp1(xi5, f5, x_query, 'linear', 'extrap');
f6_interp = interp1(xi6, f6, x_query, 'linear', 'extrap');
f7_interp = interp1(xi7, f7, x_query, 'linear', 'extrap');
f8_interp = interp1(xi8, f8, x_query, 'linear', 'extrap');
f9_interp = interp1(xi9, f9, x_query, 'linear', 'extrap');
f10_interp = interp1(xi10, f10, x_query, 'linear', 'extrap');


% 
% % Ensure the CDF is clipped between 0 and 1 (since extrapolation may go beyond)
f_interp = min(max(f_interp, 0), 1);
f2_interp = min(max(f2_interp, 0), 1);
f3_interp = min(max(f3_interp, 0), 1);
f4_interp = min(max(f4_interp, 0), 1);
f5_interp = min(max(f5_interp, 0), 1);
f6_interp = min(max(f6_interp, 0), 1);
f7_interp = min(max(f7_interp, 0), 1);
f8_interp = min(max(f8_interp, 0), 1);
f9_interp = min(max(f9_interp, 0), 1);
f10_interp = min(max(f10_interp, 0), 1);
 
plot(x_query, f_interp, 'LineWidth', 7, 'Color', '[0.5 0.5 0.5]');        % CTRL ON
hold on
plot(x_query, f2_interp, 'LineWidth', 7, 'Color', '[0.9 0.1 0.1]');   % ctrl off


% Labels and formatting with BLACK text
ax = gca;
ax.Color = 'w';          % Axes background WHITE
ax.XColor = 'k';         % X-axis ticks BLACK
ax.YColor = 'k';         % Y-axis ticks BLACK
ax.LineWidth = 1.5;      % Make axes lines visible
box on;                  % Add rectangular box around plot

xlabel('Enhancer- BRD4 hub 3D distance (nm)', 'FontSize', 24, 'FontWeight', 'bold', 'Color', 'k');
ylabel('Cumulative frequency', 'FontSize', 24, 'FontWeight', 'bold', 'Color', 'k');
set(gca,'FontSize', 30, 'FontWeight', 'bold');
xlim([x_min1 x_max1]);

h_legend = legend('Ctrl', 'RAD21 depletion','Location', 'southeast');
set(h_legend, ...
    'Color',      'w', ... % White background
    'TextColor',  'k', ... % Black text
    'EdgeColor',  'k');    % Black border

xticks(0:500:3000);