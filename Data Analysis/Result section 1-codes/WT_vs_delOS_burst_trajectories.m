clear
clc
tic
% =========================================================================
% Plot 5 MS2 trajectories from WT (CTRL) vs del OCT4/SOX2 (G67) with HMM fits
% 2 columns x 5 rows layout
% =========================================================================

% ---- USER INPUTS --------------------------------------------------------
xlsx_file = '/Users/janaa/Desktop/MS2 transcription/MS2 data sheets/Master sheets/Latest/ctrl vs g67 bursting % 7-14-2025.xlsx';

sheet_WT  = 'CTRL TRAJECTORIES';
sheet_DEL = 'G67 TRAJECTORIES';

traj_WT  = [4 6 14 20 22];     
traj_DEL = [8 12 15 20 25];   

time_int_min = 3;             % min per frame
t_max_hr     = 6;             % clip at 6 hours
y_max        = 350;           % common y-axis limit
fs_label     = 16;            % axis-label font size
fs_tick      = 13;            % tick font size
fs_title     = 18;            % column-title font size

save_fig = 'WT_vs_delOS_traj_HMM.fig';

% ---- READ DATA ----------------------------------------------------------
WT_all  = readmatrix(xlsx_file, 'Sheet', sheet_WT);
DEL_all = readmatrix(xlsx_file, 'Sheet', sheet_DEL);

WT_all  = WT_all(any(~isnan(WT_all),2), any(~isnan(WT_all),1));
DEL_all = DEL_all(any(~isnan(DEL_all),2), any(~isnan(DEL_all),1));

% Frames corresponding to 6 hours
n_frames_max = floor(t_max_hr * 60 / time_int_min) + 1;   % inclusive

% ---- BUILD FIGURE -------------------------------------------------------
n_rows = 5;
n_cols = 2;

figure('Color','w','Position',[100 100 1100 1200]);

light_red = [1.0 0.4 0.4];
trace_grn = [0.2 0.9 0.0];

for r = 1:n_rows

    % ---------- Column 1 : WT ----------
    I_WT = WT_all(:, traj_WT(r));
    I_WT = I_WT(~isnan(I_WT));
    I_WT = I_WT(1:min(end, n_frames_max));
    t_WT = (0:numel(I_WT)-1) * time_int_min / 60;          % hours

    fit_WT = HMM_fit_fun(I_WT);
    fit_WT = fit_WT{1,1};

    subplot(n_rows, n_cols, (r-1)*n_cols + 1);
    plot(t_WT, I_WT, 'Color', trace_grn, 'LineWidth', 2); hold on;
    plot(t_WT, fit_WT, 'Color', light_red, 'LineWidth', 2);
    ylabel('MS2 intensity', 'FontSize', fs_label);
    xlim([0 t_max_hr]); ylim([0 y_max]);
    xticks(0:1:t_max_hr);
    set(gca, 'FontSize', fs_tick, 'Box', 'on');
    if r == 1, title('WT', 'FontSize', fs_title); end
    if r == n_rows, xlabel('Time (hours)', 'FontSize', fs_label); end

    % ---------- Column 2 : del OCT4/SOX2 ----------
    I_DEL = DEL_all(:, traj_DEL(r));
    I_DEL = I_DEL(~isnan(I_DEL));
    I_DEL = I_DEL(1:min(end, n_frames_max));
    t_DEL = (0:numel(I_DEL)-1) * time_int_min / 60;        % hours

    fit_DEL = HMM_fit_fun(I_DEL);
    fit_DEL = fit_DEL{1,1};

    subplot(n_rows, n_cols, (r-1)*n_cols + 2);
    plot(t_DEL, I_DEL, 'Color', trace_grn, 'LineWidth', 2); hold on;
    plot(t_DEL, fit_DEL, 'Color', light_red, 'LineWidth', 2);
    ylabel('MS2 intensity', 'FontSize', fs_label);
    xlim([0 t_max_hr]); ylim([0 y_max]);
    xticks(0:1:t_max_hr);
    set(gca, 'FontSize', fs_tick, 'Box', 'on');
    if r == 1, title('del OCT4/SOX2', 'FontSize', fs_title); end
    if r == n_rows, xlabel('Time (hours)', 'FontSize', fs_label); end
end

sgtitle('MS2 intensity with HMM 2-state fitting', 'FontSize', fs_title);

savefig(save_fig);

toc