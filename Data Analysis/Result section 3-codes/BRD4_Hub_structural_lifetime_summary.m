%% CONDENSATE/HUB STRUCTURAL LIFETIME ANALYSIS

clear; clc; close all;

% ========================= USER CONFIG ===================================
data_folder   = '/Users/janaa/Desktop/MS2 transcription/Compiled data 7-20-2026/BRD4 hub tracking- 1 color/Final';   
output_folder = '/Users/janaa/Desktop/MS2 transcription/Compiled data 7-20-2026/BRD4 hub tracking- 1 color';
% =========================================================================

if ~exist(output_folder, 'dir'), mkdir(output_folder); end

%% 1. LOAD mat files
mat_files = dir(fullfile(data_folder, '*.mat'));
if isempty(mat_files), error('No .mat files in %s', data_folder); end

cell_data = struct('id',{},'name',{},'lifetimes',{},'censored',{},'n',{});
for i = 1:numel(mat_files)
    S = load(fullfile(mat_files(i).folder, mat_files(i).name));

    if ~isfield(S, 'structural_lifetime_sec')
        warning('Skipping %s — no structural_lifetime_sec.', mat_files(i).name);
        continue;
    end

    % 1) CONVERT TO MINUTES
    lt   = S.structural_lifetime_sec(:) / 60; 

    keep = isfinite(lt) & lt > 0;
    lt   = lt(keep);

    if      isfield(S, 'censored'),    c = logical(S.censored(:));    c = c(keep);
    elseif  isfield(S, 'is_censored'), c = logical(S.is_censored(:)); c = c(keep);
    else,                              c = false(size(lt));
    end

    k = numel(cell_data) + 1;
    cell_data(k).id        = k;
    cell_data(k).name      = mat_files(i).name;
    cell_data(k).lifetimes = lt;
    cell_data(k).censored  = c;
    cell_data(k).n         = numel(lt);

    fprintf('Cell %d (%s): %d tracks, %d censored\n', k, mat_files(i).name, numel(lt), sum(c));
end

n_cells  = numel(cell_data);
all_lt   = vertcat(cell_data.lifetimes);
all_c    = vertcat(cell_data.censored);
n_tracks = numel(all_lt);

fprintf('\nTotal: %d tracks across %d cells.\n\n', n_tracks, n_cells);

%% 2. SUMMARY STATS

% Pooled distribution
mean_all   = mean(all_lt);
sd_all     = std(all_lt);
median_all = median(all_lt);

% Bootstrap CI on pooled median
rng(42);
boot_med  = bootstrp(5000, @median, all_lt);
median_ci = prctile(boot_med, [2.5 97.5]);

% Per-cell means and medians
cell_means   = arrayfun(@(c) mean(c.lifetimes),   cell_data);
cell_medians = arrayfun(@(c) median(c.lifetimes), cell_data);

% Mean of per-cell means +/- SD
mean_of_means = mean(cell_means);
sd_of_means   = std(cell_means);
sem_of_means  = std(cell_means) / sqrt(n_cells);

% Mean of per-cell medians +/- SD
mean_of_medians = mean(cell_medians);
sd_of_medians   = std(cell_medians);
sem_of_medians  = std(cell_medians) / sqrt(n_cells);

fprintf('--- Pooled (all condensates) ---\n');
fprintf('  Mean   = %.2f +/- %.2f min (SD, n=%d condensates)\n', mean_all, sd_all, numel(all_lt));
fprintf('  Median = %.2f min  [95%% CI: %.2f, %.2f]\n', median_all, median_ci(1), median_ci(2));
fprintf('\n--- Per-cell means (N=%d cells) ---\n', n_cells);
fprintf('  Individual: %s\n', sprintf('%.2f  ', cell_means));
fprintf('  Mean of means   = %.2f +/- %.2f min (SD)   SEM = %.2f\n', ...
        mean_of_means, sd_of_means, sem_of_means);
fprintf('\n--- Per-cell medians (N=%d cells) ---\n', n_cells);
fprintf('  Individual: %s\n', sprintf('%.2f  ', cell_medians));
fprintf('  Mean of medians = %.2f +/- %.2f min (SD)   SEM = %.2f\n\n', ...
        mean_of_medians, sd_of_medians, sem_of_medians);

% Colors
col_blue    = [0.25 0.40 0.70];
col_dark    = [0.15 0.25 0.55];
col_red     = [0.75 0.20 0.20];
col_gray    = [0.70 0.70 0.70];
col_magenta = [0.90 0.10 0.90];

%% 3. FIGURE — HISTOGRAM + SMOOTHED KM  (1 x 2 panel)

fig = figure('Name','Lifetime_summary','Color','w','Position',[100 100 1100 480]);

% --- Panel 1: Histogram -----------------------------------------------
subplot(1,2,1);
histogram(all_lt, 'BinWidth', 5, 'FaceColor', col_magenta, ...
          'EdgeColor', 'none', 'FaceAlpha', 0.85);
hold on; yl = ylim;
plot([median_all median_all], yl, '--', 'Color', col_red, 'LineWidth', 2.2);
plot([mean_all   mean_all],   yl, '--', 'Color', 'k',     'LineWidth', 2.2);
xlabel('Structural lifetime (min)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Count',                      'FontSize', 14, 'FontWeight', 'bold');
title('Lifetime distribution',       'FontSize', 15, 'FontWeight', 'bold');
legend({sprintf('n = %d tracks', n_tracks), ...
        sprintf('Median = %.1f min', median_all), ...
        sprintf('Mean \\pm SD = %.1f \\pm %.1f min', mean_all, sd_all)}, ...
       'Location', 'best', 'Box', 'off', 'FontSize', 11, 'FontWeight', 'bold');
xlim([0 60]);
set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'TickDir', 'out', 'Box', 'off');

% --- Panel 2: Smoothed KM survival ------------------------------------
subplot(1,2,2);

bw_surv = 0.4;                                    % bandwidth (min); tune as needed
x_query = linspace(0, 60, 2000);                   % high-res x grid

% Per-cell smoothed survival curves (gray)
S_smooth_cells = nan(numel(x_query), n_cells);
for i = 1:n_cells
    lt_i = cell_data(i).lifetimes;
    if numel(lt_i) < 5; continue; end              % skip tiny cells
    [f_i, xi_i] = ksdensity(lt_i, 'Function', 'cdf', 'Bandwidth', bw_surv);
    s_interp = interp1(xi_i, 1 - f_i, x_query, 'linear', 'extrap');
    s_interp = min(max(s_interp, 0), 1);           % clip to [0  1]
    S_smooth_cells(:, i) = s_interp(:);
    plot(x_query, s_interp, 'Color', col_gray, 'LineWidth', 1.0);
    hold on;
end

% Pooled smoothed survival (magenta bold)
[f_all, xi_all] = ksdensity(all_lt, 'Function', 'cdf', 'Bandwidth', bw_surv);
S_pooled = interp1(xi_all, 1 - f_all, x_query, 'linear', 'extrap');
S_pooled = min(max(S_pooled, 0), 1);

% Mean +/- 95% CI across per-cell curves (shaded band)
n_at_x           = sum(~isnan(S_smooth_cells), 2);
mean_S           = mean(S_smooth_cells, 2, 'omitnan');
sem_S            = std(S_smooth_cells, 0, 2, 'omitnan') ./ sqrt(max(n_at_x, 1));
t_crit           = tinv(0.975, max(n_at_x - 1, 1));
ci_lo            = max(mean_S - t_crit .* sem_S, 0);
ci_up            = min(mean_S + t_crit .* sem_S, 1);
valid            = ~isnan(mean_S) & n_at_x >= 2;

fill([x_query(valid)'; flipud(x_query(valid)')], ...
     [ci_lo(valid);     flipud(ci_up(valid))], ...
     col_magenta, 'FaceAlpha', 0.20, 'EdgeColor', 'none');
plot(x_query, S_pooled, 'Color', col_magenta, 'LineWidth', 2.8);

% Median survival annotation
t50 = interp1(S_pooled, x_query, 0.5, 'linear');
if ~isnan(t50)
    plot([t50 t50], [0 0.5], ':k', 'LineWidth', 1.5);
    plot([0  t50],  [0.5 0.5], ':k', 'LineWidth', 1.5);
    text(t50 + 0.5, 0.53, sprintf('t_{50} = %.1f min', t50), ...
         'FontSize', 12, 'FontWeight', 'bold');
end

xlabel('Time (min)',              'FontSize', 14, 'FontWeight', 'bold');
ylabel('Survival probability',    'FontSize', 14, 'FontWeight', 'bold');
title('Kaplan-Meier (smoothed)',  'FontSize', 15, 'FontWeight', 'bold');
xlim([0 60]); ylim([0 1.02]);
set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'TickDir', 'out', 'Box', 'off');
grid off;

%% 4. SAVE
saveas(fig, fullfile(output_folder, 'Lifetime_summary.png'));
savefig(fig, fullfile(output_folder, 'Lifetime_summary.fig'));
fprintf('Saved to %s\n', output_folder);