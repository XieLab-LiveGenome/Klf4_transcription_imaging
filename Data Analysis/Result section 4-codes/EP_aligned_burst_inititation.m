%% Align E-P distance trajectories at burst initiation (HMM OFF -> ON)
%  HMM state is read directly from MS2_score(:,3) -- no re-fitting.
%  MS2 signal is plotted in green, E-P distance in magenta.

clear; clc; tic

%% ---------------- Parameters ----------------
folderPath = ['/Users/janaa/Desktop/MS2 transcription/Compiled data 7-20-2026/' ...
              'EP-MS2 imaging -3 color/Ctrl new'];

thr           = 1200;   % nm, reject a trajectory if any frame exceeds this (bad localization)
pre_thr       = 400;    % nm, require >=1 pre-transition frame closer than this (set Inf to disable)
frames_before = 21;     % frames kept before the ON frame
frames_after  = 9;      % frames kept from the ON frame onward
dt            = 20;      % s per frame (set to your interval; 1 = plot in frames)

%% ---------------- Collect trajectories ----------------
fileList = dir(fullfile(folderPath, '*.mat'));
nWin     = frames_before + frames_after;

MS2_all = [];   % one row per accepted trajectory
EP_all  = [];
nRej    = 0;
S       = struct('FileName', {}, 'N_accepted', {}, 'N_rejected', {});

for i = 1:numel(fileList)
    d = load(fullfile(folderPath, fileList(i).name));
    if ~isfield(d, 'MS2_score') || ~isfield(d, 'E_P_dist_3d'), continue; end

    keep = d.MS2_score(:,2) ~= 0;        % drop trailing zeros
    MS2  = d.MS2_score(keep, 2);
    hmm  = d.MS2_score(keep, 3);
    EP   = d.E_P_dist_3d(keep, 1);

    % ON = any state above the lowest one. Works whether column 3 is already
    % 0/1 or holds the HMM-fitted intensity levels. Override if needed,
    % e.g. on = hmm == 2;
    on = hmm > min(hmm) + eps;

    onsets = find(diff(on) == 1) + 1;    % index of the first ON frame of each burst

    min_off_gap = 7;                     % require this many OFF frames before onset

    nAcc = 0; nRejFile = 0;
    for k = 1:numel(onsets)
        w = (onsets(k) - frames_before) : (onsets(k) + frames_after - 1);
        if w(1) < 1 || w(end) > numel(MS2), continue; end

        % Reject if a previous burst ended within min_off_gap frames before onset
        gap_win = (onsets(k) - min_off_gap) : (onsets(k) - 1);
        if gap_win(1) < 1 || any(on(gap_win)), continue; end

        ep = EP(w);
        if any(isnan(ep)) || any(ep > thr), continue; end        % localization quality

        if any(ep(1:frames_before) < pre_thr)                    % pre-transition contact
            MS2_all = [MS2_all; MS2(w)'];  %#ok<AGROW>
            EP_all  = [EP_all;  ep'];      %#ok<AGROW>
            nAcc    = nAcc + 1;
        else
            nRejFile = nRejFile + 1;
        end
    end

    nRej = nRej + nRejFile;
    S(end+1) = struct('FileName', fileList(i).name, ...
                      'N_accepted', nAcc, 'N_rejected', nRejFile);  %#ok<AGROW>
    fprintf('%-40s  accepted: %3d   rejected: %3d\n', fileList(i).name, nAcc, nRejFile);
end

nTraj = size(MS2_all, 1);
assert(nTraj > 0, 'No trajectories passed the filters.');

%% ---------------- Per-trajectory MS2 normalization ----------------
MS2_max = max(MS2_all, [], 2);            % per-trajectory peak
MS2_max(MS2_max == 0) = 1;                % guard against divide-by-zero
MS2_all = MS2_all ./ MS2_max;             % each row now in [0, 1]

%% ---------------- Averages ----------------
t       = (-frames_before : frames_after - 1) * dt;
MS2_avg = mean(MS2_all, 1, 'omitnan');
EP_avg  = mean(EP_all,  1, 'omitnan');
MS2_sem = std(MS2_all, 0, 1, 'omitnan') / sqrt(nTraj);
EP_sem  = std(EP_all,  0, 1, 'omitnan') / sqrt(nTraj);

fprintf('\n=== SUMMARY ===\n');
disp(struct2table(S));
fprintf('Accepted: %d   Rejected: %d   (acceptance %.1f%%)\n', ...
        nTraj, nRej, 100 * nTraj / (nTraj + nRej));
fprintf('Window: %d frames (%d before, %d from ON)\n', nWin, frames_before, frames_after);

%% ---------------- Merged overlay ----------------
green   = [0.00 0.65 0.25];
magenta = [0.85 0.00 0.85];
shade   = @(y, e, c) fill([t fliplr(t)], [y - e, fliplr(y + e)], c, ...
                          'FaceAlpha', 0.25, 'EdgeColor', 'none');

figure('Position', [80 80 700 460]); hold on;

yyaxis left
shade(EP_avg, EP_sem, magenta);
plot(t, EP_avg, '-', 'Color', magenta, 'LineWidth', 2);
ylabel('E-P distance (nm)'); set(gca, 'YColor', magenta);

yyaxis right
shade(MS2_avg, MS2_sem, green);
plot(t, MS2_avg, '-', 'Color', green, 'LineWidth', 2);
ylabel('MS2 intensity (norm.)'); set(gca, 'YColor', green);

xline(0, 'k--', 'LineWidth', 1.2);
xlabel(ternary(dt == 1, 'Frame relative to burst onset', 'Time relative to burst onset (s)'));
title(sprintf('Burst initiation, n = %d trajectories', nTraj));
box off; set(gca, 'FontSize', 12);

%% ---------------- Save ----------------
save('burst_initiation_aligned.mat', 'MS2_all', 'EP_all', 'MS2_avg', 'EP_avg', ...
     'MS2_sem', 'EP_sem', 't', 'nTraj', 'nRej', 'frames_before', 'frames_after', ...
     'thr', 'pre_thr', 'S');
fprintf('\nSaved: burst_initiation_aligned.mat\n');
toc

%% ---------------- helper ----------------
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end