clear; clc; 
set(0, 'DefaultAxesFontName', 'Arial', ...
       'DefaultTextFontName', 'Arial');

%%============================== PROJECT DIRECTORY SETUP ========================
% Required for EP_contact_param.m and its helpers (traj_clean, binaryv2)
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'functions'));
%%==============================================================================

%% ==================== SHARED SETUP ============================
ctrl_file   = '/Users/janaa/Desktop/MS2 transcription/EP datasheets/Fast EP /EP/Rad21 EP ctrl vs depletion 1-7-26.xlsx';
treat_file  = '/Users/janaa/Desktop/MS2 transcription/EP datasheets/Fast EP /EP/Rad21 EP ctrl vs depletion 1-7-26.xlsx';

sheet_ctrl  = 'CTRL R21 compiled';    % sheet name (e.g. 'Sheet1') or index
sheet_treat = 'Distance compiled_R21_traj';

time_int    = 5;      % frame interval (seconds)
thr         = 120;    % contact threshold (nm)

% Labels
ctrl_label  = 'Ctrl';
treat_label = 'RAD21(-)';
fig_title   = 'RAD21(-)';
panel_label = 'F';

% Colors
color_ctrl  = [0.65 0.65 0.65];
color_treat = [0.8  0.2  0.2];

% ---- PART A parameters (average E-P distance) -----------------
x_min_hist = 0;   x_max_hist = 800;   bin_w = 20;
y_max_hist = [];                      % [] = auto

% ---- PART B plot limits --------------------------------------
x_min_OFF = 0;    x_max_OFF = 600;
x_min_D   = 0;    x_max_D   = 30;

% ---- Load once, used by both parts ---------------------------
fprintf('===== LOADING =====\n');
ctrl_block  = load_EP_block_xlsx(ctrl_file,  sheet_ctrl);
treat_block = load_EP_block_xlsx(treat_file, sheet_treat);


%% ##############################################################
%  PART A — AVERAGE E-P DISTANCE
%  Every valid frame is used. Gaps are irrelevant here.
%  ##############################################################
fprintf('\n===== PART A: AVERAGE E-P DISTANCE =====\n');

% --- A1) Pooled distances (frame-level) ------------------------
ctrl_all  = ctrl_block(:);   ctrl_all(isnan(ctrl_all))   = [];
treat_all = treat_block(:);  treat_all(isnan(treat_all)) = [];

fprintf('Pooled across all frames:\n');
fprintf('  %-6s : n = %d frames, %.1f +/- %.1f nm (mean +/- SD)\n', ...
    ctrl_label,  numel(ctrl_all),  mean(ctrl_all),  std(ctrl_all));
fprintf('  %-6s : n = %d frames, %.1f +/- %.1f nm (mean +/- SD)\n', ...
    treat_label, numel(treat_all), mean(treat_all), std(treat_all));

% --- A2) Per-trajectory means + t-test -------------------------
ctrl_traj_means  = mean(ctrl_block,  1, 'omitnan');
treat_traj_means = mean(treat_block, 1, 'omitnan');

ctrl_n_frames    = sum(~isnan(ctrl_block),  1);
treat_n_frames   = sum(~isnan(treat_block), 1);

ctrl_traj_means  = ctrl_traj_means(~isnan(ctrl_traj_means));
treat_traj_means = treat_traj_means(~isnan(treat_traj_means));

n_ctrl    = numel(ctrl_traj_means);
n_treat   = numel(treat_traj_means);
sem_ctrl  = std(ctrl_traj_means)  / sqrt(n_ctrl);
sem_treat = std(treat_traj_means) / sqrt(n_treat);

[~, p_ttest,    ci_ttest, stats_ttest]    = ttest2(ctrl_traj_means, treat_traj_means, ...
    'Alpha', 0.05, 'Vartype', 'unequal');
[~, p_ttest_eq, ~,        stats_ttest_eq] = ttest2(ctrl_traj_means, treat_traj_means, ...
    'Alpha', 0.05, 'Vartype', 'equal');
p_rank = ranksum(ctrl_traj_means, treat_traj_means);

fprintf('\nPer-trajectory averages (n = trajectories):\n');
fprintf('  %-6s : n = %2d traj, mean = %.1f nm, SD = %.1f, SEM = %.1f  [%d-%d valid frames/traj]\n', ...
    ctrl_label,  n_ctrl,  mean(ctrl_traj_means),  std(ctrl_traj_means),  sem_ctrl, ...
    min(ctrl_n_frames(ctrl_n_frames>0)),  max(ctrl_n_frames));
fprintf('  %-6s : n = %2d traj, mean = %.1f nm, SD = %.1f, SEM = %.1f  [%d-%d valid frames/traj]\n', ...
    treat_label, n_treat, mean(treat_traj_means), std(treat_traj_means), sem_treat, ...
    min(treat_n_frames(treat_n_frames>0)), max(treat_n_frames));
fprintf('  Difference (%s - %s) = %.2f nm\n', ...
    ctrl_label, treat_label, mean(ctrl_traj_means) - mean(treat_traj_means));
fprintf('  Welch''s   t-test : t(%.2f) = %.3f, p = %.4e\n', stats_ttest.df,    stats_ttest.tstat,    p_ttest);
fprintf('  Student''s t-test : t(%d)    = %.3f, p = %.4e\n', stats_ttest_eq.df, stats_ttest_eq.tstat, p_ttest_eq);
fprintf('  Mann-Whitney U   : p = %.4e\n', p_rank);
fprintf('  95%% CI of difference (Welch) : [%.2f, %.2f] nm\n', ci_ttest(1), ci_ttest(2));

% --- A3) KS test on pooled distance distributions --------------
[~, p_dist] = kstest2(ctrl_all, treat_all, 'Alpha', 0.05);
fprintf('  KS test, pooled distance distributions : p = %.4e\n', p_dist);


%% ##############################################################
%  PART B — E-P CONTACT PARAMETERS
%  Uses the original EP_contact_param.m (functions folder).
%  Gap handling and trajectory trimming are done inside traj_clean;
%  no additional gap logic is applied here.
%    out{1} = pooled ON  (contact)    durations, all cells
%    out{2} = pooled OFF (separation) durations, all cells
%    out{3} = per-cell mean ON  duration
%    out{4} = per-cell mean OFF duration
%  ##############################################################
fprintf('\n===== PART B: E-P CONTACT PARAMETERS =====\n');

out_ctrl  = EP_contact_param(ctrl_block,  thr, time_int);
out_treat = EP_contact_param(treat_block, thr, time_int);

ctrl_contact_dur  = out_ctrl{1};    ctrl_contact_off  = out_ctrl{2};
treat_contact_dur = out_treat{1};   treat_contact_off = out_treat{2};

ctrl_on_cellavg   = out_ctrl{3};    ctrl_off_cellavg  = out_ctrl{4};
treat_on_cellavg  = out_treat{3};   treat_off_cellavg = out_treat{4};

% Cells with no ON (or no OFF) events return NaN -> drop for stats
ctrl_on_cellavg   = ctrl_on_cellavg(~isnan(ctrl_on_cellavg));
treat_on_cellavg  = treat_on_cellavg(~isnan(treat_on_cellavg));
ctrl_off_cellavg  = ctrl_off_cellavg(~isnan(ctrl_off_cellavg));
treat_off_cellavg = treat_off_cellavg(~isnan(treat_off_cellavg));

% --- B2) Survival curves (pooled events) ------------------------
S_off_ctrl  = survival_prob(ctrl_contact_off);
S_off_treat = survival_prob(treat_contact_off);
S_dur_ctrl  = survival_prob(ctrl_contact_dur);
S_dur_treat = survival_prob(treat_contact_dur);

% --- B3) Statistics --------------------------------------------
% Pooled events (KS, n = events)
[~, p_off] = kstest2(ctrl_contact_off, treat_contact_off, 'Alpha', 0.05);
[~, p_dur] = kstest2(ctrl_contact_dur, treat_contact_dur, 'Alpha', 0.05);

% Per cell (Welch t-test, n = cells)
[~, p_dur_cell] = ttest2(ctrl_on_cellavg,  treat_on_cellavg,  'Alpha',0.05, 'Vartype','unequal');
[~, p_off_cell] = ttest2(ctrl_off_cellavg, treat_off_cellavg, 'Alpha',0.05, 'Vartype','unequal');

% --- B4) Duration summary statistics ---------------------------
st_dur_ctrl  = dur_stats(ctrl_contact_dur);
st_dur_treat = dur_stats(treat_contact_dur);
st_off_ctrl  = dur_stats(ctrl_contact_off);
st_off_treat = dur_stats(treat_contact_off);

sc_dur_ctrl  = dur_stats(ctrl_on_cellavg);
sc_dur_treat = dur_stats(treat_on_cellavg);
sc_off_ctrl  = dur_stats(ctrl_off_cellavg);
sc_off_treat = dur_stats(treat_off_cellavg);

fprintf('\nContact (ON) duration, pooled events:\n');
fprintf('  %-6s : n = %3d, mean = %.1f s, median = %.1f s, SD = %.1f, SEM = %.2f\n', ...
    ctrl_label,  st_dur_ctrl.n,  st_dur_ctrl.mean,  st_dur_ctrl.median,  st_dur_ctrl.sd,  st_dur_ctrl.sem);
fprintf('  %-6s : n = %3d, mean = %.1f s, median = %.1f s, SD = %.1f, SEM = %.2f\n', ...
    treat_label, st_dur_treat.n, st_dur_treat.mean, st_dur_treat.median, st_dur_treat.sd, st_dur_treat.sem);
fprintf('  KS test (pooled events) : p = %.4e\n', p_dur);

fprintf('\nContact (ON) duration, per cell:\n');
fprintf('  %-6s : n = %2d cells, mean = %.1f s, median = %.1f s, SEM = %.2f\n', ...
    ctrl_label,  sc_dur_ctrl.n,  sc_dur_ctrl.mean,  sc_dur_ctrl.median,  sc_dur_ctrl.sem);
fprintf('  %-6s : n = %2d cells, mean = %.1f s, median = %.1f s, SEM = %.2f\n', ...
    treat_label, sc_dur_treat.n, sc_dur_treat.mean, sc_dur_treat.median, sc_dur_treat.sem);
fprintf('  Welch t-test (per cell) : p = %.4e\n', p_dur_cell);

fprintf('\nSeparation (OFF) duration, pooled events:\n');
fprintf('  %-6s : n = %3d, mean = %.1f s, median = %.1f s, SD = %.1f, SEM = %.2f\n', ...
    ctrl_label,  st_off_ctrl.n,  st_off_ctrl.mean,  st_off_ctrl.median,  st_off_ctrl.sd,  st_off_ctrl.sem);
fprintf('  %-6s : n = %3d, mean = %.1f s, median = %.1f s, SD = %.1f, SEM = %.2f\n', ...
    treat_label, st_off_treat.n, st_off_treat.mean, st_off_treat.median, st_off_treat.sd, st_off_treat.sem);
fprintf('  KS test (pooled events) : p = %.4e\n', p_off);

fprintf('\nSeparation (OFF) duration, per cell:\n');
fprintf('  %-6s : n = %2d cells, mean = %.1f s, median = %.1f s, SEM = %.2f\n', ...
    ctrl_label,  sc_off_ctrl.n,  sc_off_ctrl.mean,  sc_off_ctrl.median,  sc_off_ctrl.sem);
fprintf('  %-6s : n = %2d cells, mean = %.1f s, median = %.1f s, SEM = %.2f\n', ...
    treat_label, sc_off_treat.n, sc_off_treat.mean, sc_off_treat.median, sc_off_treat.sem);
fprintf('  Welch t-test (per cell) : p = %.4e\n', p_off_cell);


%% ##############################################################
%  FIGURE — panel 1 from Part A, panels 2-3 from Part B
%  ##############################################################
fs_lab  = 15;   fs_tick = 13;   fs_inset = 16;   fs_stat = 11.5;
lw_axis = 1.5;
col_p   = [0.55 0.55 0.55];

figure('Color','w', 'Position',[100 40 620 1180]);

% --- Panel 1: distance histogram (PART A) ----------------------
ax1 = axes('Position',[0.17 0.715 0.76 0.215]); hold(ax1,'on');
histogram(ax1, ctrl_all, 'Normalization','probability', 'DisplayStyle','stairs', ...
    'LineWidth',3, 'BinWidth',bin_w, 'EdgeColor',color_ctrl, 'DisplayName',ctrl_label);
fitGaussian_histogram(ctrl_all,  x_min_hist, x_max_hist, bin_w, color_ctrl,  '-');
histogram(ax1, treat_all, 'Normalization','probability', 'DisplayStyle','stairs', ...
    'LineWidth',3, 'BinWidth',bin_w, 'EdgeColor',color_treat, 'DisplayName',treat_label);
fitGaussian_histogram(treat_all, x_min_hist, x_max_hist, bin_w, color_treat, '-');

xlim(ax1,[x_min_hist x_max_hist]);
if isempty(y_max_hist); y_max_hist = ceil(max(ax1.YLim)*100)/100; end
ylim(ax1,[0 y_max_hist]);
set(ax1, 'YTick',0:0.01:y_max_hist, 'XTick',x_min_hist:100:x_max_hist, ...
         'FontSize',fs_tick, 'FontWeight','bold', 'LineWidth',lw_axis, 'Box','on', 'TickDir','in');
xlabel(ax1,'E-P 3D distance (nm)', 'FontSize',fs_lab, 'FontWeight','bold');
ylabel(ax1,'Probability',          'FontSize',fs_lab, 'FontWeight','bold');
title(ax1, fig_title, 'FontSize',18, 'FontWeight','bold');
legend(ax1, 'Location','northeast', 'FontSize',13, 'Box','off');

text(ax1, 0.97, 0.60, { sprintf('p (pooled)   = %s', pnum(p_dist)), ...
                         sprintf('p (per cell) = %s', pnum(p_ttest)) }, ...
    'Units','normalized', 'HorizontalAlignment','right', 'VerticalAlignment','top', ...
    'FontSize',fs_stat, 'FontWeight','bold', 'Color',col_p);

% --- Panel 2: separation duration survival (PART B) ------------
ax2 = axes('Position',[0.17 0.395 0.76 0.225]); hold(ax2,'on');
plot(ax2, S_off_ctrl{2},  S_off_ctrl{1},  'Color',color_ctrl,  'LineWidth',3);
plot(ax2, S_off_treat{2}, S_off_treat{1}, 'Color',color_treat, 'LineWidth',3);
xlim(ax2,[x_min_OFF x_max_OFF]); ylim(ax2,[0 1]);
set(ax2, 'YTick',0:0.2:1, 'FontSize',fs_tick, 'FontWeight','bold', ...
         'LineWidth',lw_axis, 'Box','on', 'TickDir','in');
xlabel(ax2,'Time (s)',             'FontSize',fs_lab, 'FontWeight','bold');
ylabel(ax2,'Survival Probability', 'FontSize',fs_lab, 'FontWeight','bold');
text(ax2, 0.55, 0.93, {'E-P separate','duration'}, 'Units','normalized', ...
    'HorizontalAlignment','center', 'VerticalAlignment','top', ...
    'FontSize',fs_inset, 'FontWeight','bold');
text(ax2, 0.95, 0.20, pstr(p_off), 'Units','normalized', ...
    'HorizontalAlignment','right', 'VerticalAlignment','middle', ...
    'FontSize',fs_inset, 'FontWeight','bold', 'Color',col_p);

text(ax2, 0.95, 0.62, 'mean / median (s)', 'Units','normalized', ...
    'HorizontalAlignment','right', 'FontSize',fs_stat, 'FontWeight','bold', 'Color',col_p);
text(ax2, 0.95, 0.52, sprintf('%s: %.0f / %.0f', ctrl_label, st_off_ctrl.mean, st_off_ctrl.median), ...
    'Units','normalized', 'HorizontalAlignment','right', ...
    'FontSize',fs_stat, 'FontWeight','bold', 'Color',color_ctrl);
text(ax2, 0.95, 0.42, sprintf('%s: %.0f / %.0f', treat_label, st_off_treat.mean, st_off_treat.median), ...
    'Units','normalized', 'HorizontalAlignment','right', ...
    'FontSize',fs_stat, 'FontWeight','bold', 'Color',color_treat);

% --- Panel 3: contact duration survival (PART B) ---------------
ax3 = axes('Position',[0.17 0.075 0.76 0.225]); hold(ax3,'on');
plot(ax3, S_dur_ctrl{2},  S_dur_ctrl{1},  'Color',color_ctrl,  'LineWidth',3);
plot(ax3, S_dur_treat{2}, S_dur_treat{1}, 'Color',color_treat, 'LineWidth',3);
xlim(ax3,[x_min_D x_max_D]); ylim(ax3,[0 1]);
set(ax3, 'YTick',0:0.2:1, 'XTick',x_min_D:5:x_max_D, 'FontSize',fs_tick, ...
         'FontWeight','bold', 'LineWidth',lw_axis, 'Box','on', 'TickDir','in');
xlabel(ax3,'Time (s)',             'FontSize',fs_lab, 'FontWeight','bold');
ylabel(ax3,'Survival Probability', 'FontSize',fs_lab, 'FontWeight','bold');
text(ax3, 0.55, 0.93, {'E-P contact','duration'}, 'Units','normalized', ...
    'HorizontalAlignment','center', 'VerticalAlignment','top', ...
    'FontSize',fs_inset, 'FontWeight','bold');
text(ax3, 0.95, 0.20, pstr(p_dur), 'Units','normalized', ...
    'HorizontalAlignment','right', 'VerticalAlignment','middle', ...
    'FontSize',fs_inset, 'FontWeight','bold', 'Color',col_p);

text(ax3, 0.95, 0.62, 'mean / median (s)', 'Units','normalized', ...
    'HorizontalAlignment','right', 'FontSize',fs_stat, 'FontWeight','bold', 'Color',col_p);
text(ax3, 0.95, 0.52, sprintf('%s: %.1f / %.1f', ctrl_label, st_dur_ctrl.mean, st_dur_ctrl.median), ...
    'Units','normalized', 'HorizontalAlignment','right', ...
    'FontSize',fs_stat, 'FontWeight','bold', 'Color',color_ctrl);
text(ax3, 0.95, 0.42, sprintf('%s: %.1f / %.1f', treat_label, st_dur_treat.mean, st_dur_treat.median), ...
    'Units','normalized', 'HorizontalAlignment','right', ...
    'FontSize',fs_stat, 'FontWeight','bold', 'Color',color_treat);

if ~isempty(panel_label)
    annotation('textbox',[0.01 0.945 0.08 0.05], 'String',panel_label, ...
        'FontSize',30, 'FontWeight','bold', 'EdgeColor','none', ...
        'HorizontalAlignment','left', 'VerticalAlignment','middle');
end

% exportgraphics(gcf, 'EP_contact_panel.pdf', 'ContentType','vector');


%% ==================== SAVE COMPILED RESULTS ===================
out_dir   = fileparts(treat_file);
safe_lab  = regexprep(treat_label, '[^\w-]', '_');
save_name = fullfile(out_dir, sprintf('%s_compiled.mat', safe_lab));

C = struct();

% ---- metadata -------------------------------------------------
C.meta.date_run    = datetime('now');
C.meta.ctrl_file   = ctrl_file;
C.meta.treat_file  = treat_file;
C.meta.ctrl_label  = ctrl_label;
C.meta.treat_label = treat_label;
C.meta.time_int_s  = time_int;
C.meta.thr_nm      = thr;
C.meta.method      = 'EP_contact_param.m (traj_clean gap handling)';

% ---- raw blocks -----------------------------------------------
C.raw.ctrl_block  = ctrl_block;
C.raw.treat_block = treat_block;

% ---- PART A: distances ----------------------------------------
C.partA.ctrl.pooled       = ctrl_all;
C.partA.ctrl.traj_means   = ctrl_traj_means(:);
C.partA.ctrl.n_traj       = n_ctrl;
C.partA.ctrl.mean_pooled  = mean(ctrl_all);
C.partA.ctrl.sd_pooled    = std(ctrl_all);
C.partA.ctrl.mean_percell = mean(ctrl_traj_means);
C.partA.ctrl.sd_percell   = std(ctrl_traj_means);
C.partA.ctrl.sem_percell  = sem_ctrl;

C.partA.treat.pooled       = treat_all;
C.partA.treat.traj_means   = treat_traj_means(:);
C.partA.treat.n_traj       = n_treat;
C.partA.treat.mean_pooled  = mean(treat_all);
C.partA.treat.sd_pooled    = std(treat_all);
C.partA.treat.mean_percell = mean(treat_traj_means);
C.partA.treat.sd_percell   = std(treat_traj_means);
C.partA.treat.sem_percell  = sem_treat;

C.partA.stats.p_KS_pooled   = p_dist;
C.partA.stats.p_ttest_welch = p_ttest;
C.partA.stats.p_ttest_equal = p_ttest_eq;
C.partA.stats.p_ranksum     = p_rank;
C.partA.stats.ci_welch      = ci_ttest;
C.partA.stats.tstat_welch   = stats_ttest.tstat;
C.partA.stats.df_welch      = stats_ttest.df;

% ---- PART B: contact kinetics ---------------------------------
C.partB.ctrl.contact_dur      = ctrl_contact_dur;
C.partB.ctrl.off_dur          = ctrl_contact_off;
C.partB.ctrl.on_cellavg       = ctrl_on_cellavg;
C.partB.ctrl.off_cellavg      = ctrl_off_cellavg;
C.partB.ctrl.stats_dur_pooled = st_dur_ctrl;
C.partB.ctrl.stats_off_pooled = st_off_ctrl;
C.partB.ctrl.stats_dur_cell   = sc_dur_ctrl;
C.partB.ctrl.stats_off_cell   = sc_off_ctrl;
C.partB.ctrl.survival_dur     = S_dur_ctrl;
C.partB.ctrl.survival_off     = S_off_ctrl;

C.partB.treat.contact_dur      = treat_contact_dur;
C.partB.treat.off_dur          = treat_contact_off;
C.partB.treat.on_cellavg       = treat_on_cellavg;
C.partB.treat.off_cellavg      = treat_off_cellavg;
C.partB.treat.stats_dur_pooled = st_dur_treat;
C.partB.treat.stats_off_pooled = st_off_treat;
C.partB.treat.stats_dur_cell   = sc_dur_treat;
C.partB.treat.stats_off_cell   = sc_off_treat;
C.partB.treat.survival_dur     = S_dur_treat;
C.partB.treat.survival_off     = S_off_treat;

C.partB.stats.p_KS_contact_dur    = p_dur;
C.partB.stats.p_KS_off_dur        = p_off;
C.partB.stats.p_ttest_contact_cell = p_dur_cell;
C.partB.stats.p_ttest_off_cell     = p_off_cell;

save(save_name, '-struct', 'C', '-v7.3');
fprintf('\nSaved compiled results to:\n  %s\n', save_name);


%% ==================== LOCAL FUNCTIONS =========================
% NOTE: EP_contact_param, traj_clean and binaryv2 come from the
% 'functions' folder on the path — they are NOT redefined here.

function s = pstr(p)
    if p >= 1e-4
        s = sprintf('p=%.4f', p);
    else
        s = sprintf('p=%.2e', p);
    end
end


function s = pnum(p)
    if p >= 1e-4
        s = sprintf('%.4f', p);
    else
        s = sprintf('%.2e', p);
    end
end


function s = dur_stats(v)
    v = v(:); v = v(~isnan(v));
    s.n      = numel(v);
    s.mean   = mean(v);
    s.median = median(v);
    s.sd     = std(v);
    s.sem    = s.sd / sqrt(max(s.n,1));
    s.values = v;
end


function block = load_EP_block_xlsx(file_path, sheet)
% Each column = one trajectory, row 1 = header. Reads cell-by-cell so that
% text-formatted numbers, mixed types and blank cells are all handled.
% 0-valued and empty cells -> NaN. All-NaN columns are reported and removed.

    if ~isfile(file_path)
        error('Excel file not found: %s', file_path);
    end

    raw = readcell(file_path, 'Sheet', sheet);
    if size(raw,1) < 2
        error('Sheet has no data rows below the header: %s', file_path);
    end

    hdr = raw(1,:);
    dat = raw(2:end,:);
    [nr, nc] = size(dat);

    block = NaN(nr, nc);
    for c = 1:nc
        for r = 1:nr
            v = dat{r,c};
            if isnumeric(v) && isscalar(v) && ~isnan(v)
                block(r,c) = v;
            elseif ischar(v) || isstring(v)
                x = str2double(v);
                if ~isnan(x); block(r,c) = x; end
            end
        end
    end

    n_zeros = sum(block(:) == 0);
    block(block == 0) = NaN;

    col_valid = sum(~isnan(block), 1);
    dead = find(col_valid == 0);
    [~, fname] = fileparts(file_path);
    if ~isempty(dead)
        names = cell(1, numel(dead));
        for i = 1:numel(dead)
            h = hdr{dead(i)};
            if ischar(h) || isstring(h); names{i} = char(h);
            elseif isnumeric(h) && isscalar(h); names{i} = num2str(h);
            else; names{i} = '<blank>'; end
        end
        fprintf('  [%s] WARNING: %d column(s) had no usable numbers and were removed: %s\n', ...
            fname, numel(dead), strjoin(compose('col %d (%s)', dead(:), string(names(:))), ', '));
    end
    block(:, dead)  = [];
    col_valid(dead) = [];

    last_row = find(any(~isnan(block), 2), 1, 'last');
    if ~isempty(last_row); block = block(1:last_row, :); end

    fprintf('  [%s] %d trajectories x %d frames | valid frames/traj: min %d, median %d, max %d | %d zeros dropped\n', ...
        fname, size(block,2), size(block,1), min(col_valid), round(median(col_valid)), ...
        max(col_valid), n_zeros);
end


function S = survival_prob(durations)
    durations = durations(:);
    durations = durations(~isnan(durations));
    [f, x] = ecdf(durations);
    S = {1 - f, x};
end


function h_fit = fitGaussian_histogram(data, x_min, x_max, bin_w, color_rgb, ls)
    data = data(:);
    data = data(~isnan(data));

    edges   = x_min:bin_w:x_max;
    centres = edges(1:end-1) + bin_w/2;
    counts  = histcounts(data, edges, 'Normalization','probability');

    mu0 = mean(data); sig0 = std(data); amp0 = max(counts);

    gauss = @(p,x) p(1) .* exp(-((x - p(2)).^2) ./ (2*p(3)^2));
    opts  = optimset('Display','off');
    p_fit = lsqcurvefit(gauss, [amp0, mu0, sig0], centres, counts, ...
                        [0, x_min, 1], [Inf, x_max, Inf], opts);

    x_fit = linspace(x_min, x_max, 300);
    h_fit = plot(x_fit, gauss(p_fit, x_fit), 'Color', color_rgb, 'LineStyle', ls, ...
                 'LineWidth', 5, 'HandleVisibility','off');
end