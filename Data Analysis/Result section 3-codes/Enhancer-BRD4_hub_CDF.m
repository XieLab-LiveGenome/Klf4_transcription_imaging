%% GENERATE CDF PLOTS OF ENHANCER-BRD4 HUB DISTANCES FOR DIFFERENT CONDITIONS

clear
clc

% output_file = 'Rad21 compiled';
file_path = '/Users/janaa/Desktop/MS2 transcription/Condensate datasheets/Condensate_all perturbations.xlsx';


CTRL_dist = readmatrix(file_path,'Sheet','CTRL ALL','Range','C2:C500');
CTRL_dist = CTRL_dist(~isnan(CTRL_dist));

P300_dist = readmatrix(file_path,'Sheet','P300 ALL','Range','C2:C500');
P300_dist = P300_dist(~isnan(P300_dist));

JQ1_dist = readmatrix(file_path,'Sheet','JQ1 ALL','Range','C2:C500');
JQ1_dist = JQ1_dist(~isnan(JQ1_dist));

BRG_dist = readmatrix(file_path,'Sheet','BRG ALL','Range','C2:C500');
BRG_dist = BRG_dist(~isnan(BRG_dist));

CTRL_dist_4C = readmatrix(file_path,'Sheet','CTRL 4C ALL','Range','C2:C500');
CTRL_dist_4C = CTRL_dist_4C(~isnan(CTRL_dist_4C));

G67_dist = readmatrix(file_path,'Sheet','G67 ALL','Range','C2:C500');
G67_dist = G67_dist(~isnan(G67_dist));

R21_dist = readmatrix(file_path,'Sheet','R21 ALL','Range','C2:C500');
R21_dist = R21_dist(~isnan(R21_dist));

R21_G67_dist = readmatrix(file_path,'Sheet','R21_G67 ALL','Range','C2:C500');
R21_G67_dist = R21_G67_dist(~isnan(R21_G67_dist));

MED14_dist = readmatrix(file_path,'Sheet','MED14 ALL','Range','C2:C500');
MED14_dist = MED14_dist(~isnan(MED14_dist));

% 
[s1 p1]=kstest2(CTRL_dist,P300_dist,'Alpha',0.05);
% 
[s2 p2]=kstest2(CTRL_dist,JQ1_dist,'Alpha',0.05);
% 
[s3 p3]=kstest2(CTRL_dist,BRG_dist,'Alpha',0.05);
% 
[s4 p4]=kstest2(CTRL_dist_4C,G67_dist,'Alpha',0.05);
% 
[s5 p5]=kstest2(CTRL_dist_4C,R21_dist,'Alpha',0.05);

[s6 p6]=kstest2(CTRL_dist_4C,R21_G67_dist,'Alpha',0.05);

[s7 p7]=kstest2(CTRL_dist_4C,MED14_dist,'Alpha',0.05);


% % 
% % Stats_all_dist = [p1 p2 p3 p4 p5];

% 
% % figure;
% % % % %%%% ======================ON DURATION CDF================
% 
bw=0.1;
[f, xi] = ksdensity(CTRL_dist, 'function', 'cdf', 'Bandwidth', bw);  %
% plot(xi, f, 'LineWidth', 8, 'Color', 'k');
% hold on

[f2, xi2] = ksdensity(P300_dist, 'function', 'cdf', 'Bandwidth', bw);  % 
% plot(xi2, f2, 'LineWidth', 8, 'Color', '[0.6 0.6 0.6]');
% hold on

[f3, xi3] = ksdensity(JQ1_dist, 'function', 'cdf', 'Bandwidth', bw);  % Try different values
% plot(xi3, f3, 'LineWidth', 8, 'Color', '[1 0.6 0.6]');
% hold on

[f4, xi4] = ksdensity(BRG_dist, 'function', 'cdf', 'Bandwidth', bw);  % Try different values
% plot(xi4, f4, 'LineWidth', 8, 'Color', '[0 0.7 0.1]');
% hold on
% 
[f5, xi5] = ksdensity(CTRL_dist_4C, 'function', 'cdf', 'Bandwidth', bw);  % Try different values
plot(xi5, f5, 'LineWidth', 8, 'Color', '[0.7 0.1 0.8]');
hold on

[f6, xi6] = ksdensity(G67_dist, 'function', 'cdf', 'Bandwidth', bw);  % Try different values
plot(xi6, f6, 'LineWidth', 8, 'Color', 'b');
hold on

[f7, xi7] = ksdensity(R21_dist, 'function', 'cdf', 'Bandwidth', bw);  % Try different values

[f8, xi8] = ksdensity(R21_G67_dist, 'function', 'cdf', 'Bandwidth', bw);  % Try different values

[f9, xi9] = ksdensity(MED14_dist, 'function', 'cdf', 'Bandwidth', bw);  % Try different values

% plot(xi6, f6, 'LineWidth', 8, 'Color', 'b');
% hold on
% % 
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
f5_interp = interp1(xi5, f5, x_query, 'linear', 'extrap');
f6_interp = interp1(xi6, f6, x_query, 'linear', 'extrap');
f7_interp = interp1(xi7, f7, x_query, 'linear', 'extrap');
f8_interp = interp1(xi8, f8, x_query, 'linear', 'extrap');
f9_interp = interp1(xi9, f9, x_query, 'linear', 'extrap');
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

% 

figure
plot(x_query, f5_interp, 'LineWidth', 8, 'Color', '[0.6 0.6 0.6]');    %CTRL
hold on
plot(x_query, f9_interp, 'LineWidth', 8, 'Color', '[0.9 0.2 0.1]');    %R21



% Labels and formatting
xlabel('Enhancer- BRD4 condensate distance (nm)', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('Cumulative frequency', 'FontSize', 24, 'FontWeight', 'bold');
set(gca,'FontSize', 20, 'FontWeight', 'bold');
xlim([x_min1 x_max1]);

legend('Ctrl', 'MED14 depl.', 'Location', 'southeast','FontSize', 24, 'FontWeight', 'bold');


xticks(0:1000:6000);

