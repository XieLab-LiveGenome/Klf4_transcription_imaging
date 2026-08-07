clear
clc

file_path = '/Volumes/xiel2lab/Aniket/Condensate datasheets/Condensate size mat files /concatenated_cond_size_filtered.mat';
% file_path2 = '/Users/janaa/Desktop/MS2 transcription/Condensate datasheets/brd4 ctrl condensate sizes.xlsx';

data  = load(file_path, 'condensate_feret', 'condensate_AR','condensate_FWHM','condensate_FWHM_r');

ctrl_BRD4 = data.condensate_feret;
ctrl_BRD4_AR = data.condensate_AR;

% ctrl_M14 = readmatrix(file_path,'Sheet','CTRL_ALL', 'Range','b2:e250');
% ctrl_M14 = ctrl_M14(:);
% ctrl_M14 = ctrl_M14(~isnan(ctrl_M14));
% 
% ctrl_BRD4 = readmatrix(file_path2,'Sheet','CTRL_ALL', 'Range','b2:b500');
% ctrl_BRD4 = ctrl_BRD4(:);
% ctrl_BRD4 = ctrl_BRD4(~isnan(ctrl_BRD4));

figure

x_min = 0;
x_max = 1000;
s=25;

% h_ctrl_M14 = histogram(ctrl_M14, 'LineWidth', 4.0, 'BinWidth', s, 'FaceColor', '[0.0 0.7 0.7]');
% xlim([x_min x_max]);

h_ctrl_BRD4 = histogram(ctrl_BRD4, 'LineWidth', 4.0, 'BinWidth', s, 'FaceColor', '[0.7 0.0 0.7]');
xlim([x_min x_max]);

xlabel('Condensate diameter (nm)', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Counts', 'FontSize', 20, 'FontWeight', 'bold');
set(gca,'FontSize', 20, 'FontWeight', 'bold');