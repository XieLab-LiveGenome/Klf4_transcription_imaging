clear
clc
file_path = '/Users/janaa/Desktop/MS2 transcription/4 color compiled/4 COLOR_6-29-2026.xlsx';

% ---- Read all three columns together and filter NaNs jointly ----
% D = MS2 intensity, E = enhancer-hub distance, K = promoter-hub distance
raw = readmatrix(file_path, 'Sheet', 'ms intensity vs distance', 'Range', 'D2:K300');
CTRL_MS2    = raw(:, 1);   % column D
CTRL_dist_e = raw(:, 2);   % column E
CTRL_dist_p = raw(:, 8);   % column K

good = ~isnan(CTRL_MS2) & ~isnan(CTRL_dist_e) & ~isnan(CTRL_dist_p);
CTRL_MS2    = CTRL_MS2(good);
CTRL_dist_e = CTRL_dist_e(good);
CTRL_dist_p = CTRL_dist_p(good);

fprintf('Retained %d of %d rows after joint NaN removal.\n', nnz(good), numel(good));

% ---- Thresholding (rows now aligned) ----
thr = 0.2 * 65536;  % MS2 intensity threshold for BFP

CTRL_ON_dist_e  = CTRL_dist_e(CTRL_MS2 >= thr);
CTRL_OFF_dist_e = CTRL_dist_e(CTRL_MS2 <  thr);
CTRL_ON_dist_p  = CTRL_dist_p(CTRL_MS2 >= thr);
CTRL_OFF_dist_p = CTRL_dist_p(CTRL_MS2 <  thr);

[s1, p1] = kstest2(CTRL_ON_dist_e, CTRL_OFF_dist_e, 'Alpha', 0.05);
[s2, p2] = kstest2(CTRL_ON_dist_p, CTRL_OFF_dist_p, 'Alpha', 0.05);

fprintf('Enhancer-hub  KS: p = %.3g (n_ON=%d, n_OFF=%d)\n', p1, numel(CTRL_ON_dist_e), numel(CTRL_OFF_dist_e));
fprintf('Promoter-hub  KS: p = %.3g (n_ON=%d, n_OFF=%d)\n', p2, numel(CTRL_ON_dist_p), numel(CTRL_OFF_dist_p));

% ---- Common CDF computation ----
bw     = 0.1;
x_min1 = 0;
x_max1 = 3000;
x_query = linspace(x_min1, x_max1, 2000);

    function y = cdf_on_grid(data, x, bw)
        [f, xi] = ksdensity(data, 'function', 'cdf', 'Bandwidth', bw);
        y = interp1(xi, f, x, 'linear', 'extrap');
        y = min(max(y, 0), 1);
    end

f_e_on  = cdf_on_grid(CTRL_ON_dist_e,  x_query, bw);
f_e_off = cdf_on_grid(CTRL_OFF_dist_e, x_query, bw);
f_p_on  = cdf_on_grid(CTRL_ON_dist_p,  x_query, bw);
f_p_off = cdf_on_grid(CTRL_OFF_dist_p, x_query, bw);

% ---- Figure with two panels ----
fig = figure;
set(fig, 'Color', 'w');
set(fig, 'Position', [100, 100, 1600, 450]);   % wider for 2 panels

col_on  = [0.5 0.5 0.5];
col_off = [0.1 0.7 0.1];

% --- Panel 1: enhancer-hub ---
subplot(1, 2, 1);
plot(x_query, f_e_on,  'LineWidth', 7, 'Color', col_on);  hold on
plot(x_query, f_e_off, 'LineWidth', 7, 'Color', col_off);

ax = gca;
ax.Color = 'w'; ax.XColor = 'k'; ax.YColor = 'k'; ax.LineWidth = 1.5;
box on;
xlabel('Enhancer- BRD4 hub 3D distance (nm)', 'FontSize', 22, 'FontWeight', 'bold', 'Color', 'k');
ylabel('Cumulative frequency', 'FontSize', 22, 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'FontSize', 24, 'FontWeight', 'bold');
xlim([x_min1 x_max1]); ylim([0 1]);
xticks(0:500:3000);

h1 = legend('gene ON', 'gene OFF', 'Location', 'southeast');
set(h1, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', 'k');

% --- Panel 2: promoter-hub ---
subplot(1, 2, 2);
plot(x_query, f_p_on,  'LineWidth', 7, 'Color', col_on);  hold on
plot(x_query, f_p_off, 'LineWidth', 7, 'Color', col_off);

ax = gca;
ax.Color = 'w'; ax.XColor = 'k'; ax.YColor = 'k'; ax.LineWidth = 1.5;
box on;
xlabel('Promoter- BRD4 hub 3D distance (nm)', 'FontSize', 22, 'FontWeight', 'bold', 'Color', 'k');
ylabel('Cumulative frequency', 'FontSize', 22, 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'FontSize', 24, 'FontWeight', 'bold');
xlim([x_min1 x_max1]); ylim([0 1]);
xticks(0:500:3000);

h2 = legend('gene ON', 'gene OFF', 'Location', 'southeast');
set(h2, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', 'k');