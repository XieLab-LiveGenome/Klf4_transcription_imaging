clear
clc

% --- Data Loading (Control only) ---
file_path = '/Users/janaa/Desktop/MS2 transcription/4 color compiled/4 COLOR_6-29-2026.xlsx';

EC_dist_ctrl = readmatrix(file_path,'Sheet','ms intensity vs distance', 'Range','E3:E300');
EP_dist_ctrl = readmatrix(file_path,'Sheet','ms intensity vs distance', 'Range','G3:G300');
PC_dist_ctrl = readmatrix(file_path,'Sheet','ms intensity vs distance', 'Range','F3:F300');

% --- Plot Setup ---
ctrl_color = [0.400, 0.400, 0.400];   % Dark grey

figure('Name', 'Spatial Distances - Control', 'Position', [100, 100, 1200, 380]);

% Pearson correlations
r_EC_PC = corrcoef(EC_dist_ctrl, PC_dist_ctrl, 'Rows', 'complete');
r_EC_EP = corrcoef(EC_dist_ctrl, EP_dist_ctrl, 'Rows', 'complete');
r_PC_EP = corrcoef(PC_dist_ctrl, EP_dist_ctrl, 'Rows', 'complete');

% --- Panel 1: E-C vs P-C ---
subplot(1, 3, 1);
scatter(EC_dist_ctrl, PC_dist_ctrl, 30, ctrl_color);
hold on;
valid = ~isnan(EC_dist_ctrl) & ~isnan(PC_dist_ctrl);
if any(valid)
    p = polyfit(EC_dist_ctrl(valid), PC_dist_ctrl(valid), 1);
    x_fit = [min(EC_dist_ctrl(valid)), max(EC_dist_ctrl(valid))];
    plot(x_fit, polyval(p, x_fit), 'Color', ctrl_color, 'LineStyle', ':', 'LineWidth', 2);
end
hold off;
xlim([0 3000]); ylim([0 3000]);
xlabel('E-C 3D distance (nm)', 'FontWeight', 'bold');
ylabel('P-C 3D distance (nm)', 'FontWeight', 'bold');
text(250, 4600, sprintf('Ctrl (r = %.2f)', r_EC_PC(1,2)), ...
    'FontSize', 11, 'FontWeight', 'bold', 'Color', ctrl_color);

% --- Panel 2: E-C vs E-P ---
subplot(1, 3, 2);
scatter(EC_dist_ctrl, EP_dist_ctrl, 30, ctrl_color);
hold on;
valid = ~isnan(EC_dist_ctrl) & ~isnan(EP_dist_ctrl);
if any(valid)
    p = polyfit(EC_dist_ctrl(valid), EP_dist_ctrl(valid), 1);
    x_fit = [min(EC_dist_ctrl(valid)), max(EC_dist_ctrl(valid))];
    plot(x_fit, polyval(p, x_fit), 'Color', ctrl_color, 'LineStyle', ':', 'LineWidth', 2);
end
hold off;
xlim([0 3000]); ylim([0 1500]);
xlabel('E-C 3D distance (nm)', 'FontWeight', 'bold');
ylabel('E-P 3D distance (nm)', 'FontWeight', 'bold');
text(250, 2300, sprintf('Ctrl (r = %.2f)', r_EC_EP(1,2)), ...
    'FontSize', 11, 'FontWeight', 'bold', 'Color', ctrl_color);

% --- Panel 3: P-C vs E-P ---
subplot(1, 3, 3);
scatter(PC_dist_ctrl, EP_dist_ctrl, 30, ctrl_color);
hold on;
valid = ~isnan(PC_dist_ctrl) & ~isnan(EP_dist_ctrl);
if any(valid)
    p = polyfit(PC_dist_ctrl(valid), EP_dist_ctrl(valid), 1);
    x_fit = [min(PC_dist_ctrl(valid)), max(PC_dist_ctrl(valid))];
    plot(x_fit, polyval(p, x_fit), 'Color', ctrl_color, 'LineStyle', ':', 'LineWidth', 2);
end
hold off;
xlim([0 3000]); ylim([0 1500]);
xlabel('P-C 3D distance (nm)', 'FontWeight', 'bold');
ylabel('E-P 3D distance (nm)', 'FontWeight', 'bold');
text(250, 2300, sprintf('Ctrl (r = %.2f)', r_PC_EP(1,2)), ...
    'FontSize', 11, 'FontWeight', 'bold', 'Color', ctrl_color);

sgtitle('Enhancer, Promoter, and Condensate Spatial Distances (Ctrl)', ...
    'FontSize', 16, 'FontWeight', 'bold');