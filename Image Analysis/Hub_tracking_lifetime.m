% % % % BRD4_TRACKING_PIPELINE.m
% % %  Two-channel SR-SIM timecourse:
% % %    Channel 1 = SIM BRD4 signal      used for peak detection & lifetimes
% % %    Channel 2 = widefield BRD4       used for nuclear segmentation
% % %
% % %  Pipeline:
% % %    Step 1 : Load CZI (both channels) 
% % %    Step 2 : Nuclear segmentation on widefield, centroid-tracking
% % %    Step 3 : Detect intensity-defective frames (SIM channel)
% % %    Step 4 : Multiplicative intensity correction of defective frames
% % %    Step 5 : Condensate peak detection on SIM (pkfnd + cntrd)
% % %    Step 6 : Rigid-body registration (translation + gated rotation)
% % %    Step 7 : Hungarian linking (matchpairs, memory = 1)
% % %    Step 8 : Birth / death with NaN edge-censoring
% % %    Step 9 : Histogram + Kaplan-Meier survival
% % %    Step 10: Save .mat
% % %
% % %  Dependencies on path:
% % %    ReadImage6D2, cellpose_seg2, pkfnd, cntrd, matchpairs (R2019b+)

clc; clear; close all;

% %============================== PROJECT DIRECTORY SETUP=============================
here = fileparts(mfilename('fullpath'));   % Folder of this main script
addpath(fullfile(here, 'functions'));      % Folder with the functions
addpath(fullfile(here, 'bioformats'));     % Folder with the bioformats package for image loading
addpath(fullfile(here, 'HMM fitting'));    % Folder with the HMM fitting package
javaaddpath(fullfile(here, 'bioformats', 'bioformats_package.jar'));    % Add the JAVA path for loading bioformats
% %==============================%%==============================%%===================

% 
% % =====================================================================
% %                              PARAMETERS
% % =====================================================================
p = struct();
p.NUC_DIAMETER_PX  = 120;        % this is the value of the nucleus diameter fed under 0.3x downsampling
p.MIN_NUC_AREA     = 30000;
p.MAX_CENT_JUMP    = 80;
p.PEAK_SEPARATION  = 25;
p.PEAK_FEATSIZE    = 5;
p.PEAK_GAUSS_SIGMA = 1.0;
p.THR_MULT         = 3.0;
p.MAX_DISP_PX      = 45;
p.MEMORY           = 1;           % Buffer of 1 frame 
p.DEFECT_Z         = -2.0;
p.ROT_DEG_GATE     = 3.0;
p.ROT_ECC_GATE     = 0.60;
p.INTERVAL_SEC     = 120;

p.CH_SIM           = 1;     % channel index: SIM BRD4 (peak detection)
p.CH_WF            = 2;     % channel index: widefield (segmentation)

% Input paths
MIP_filename = '/Volumes/Aniket2/BRD4 single color 5-14-2026/SIM+WF/ctrl 1 2 minute _SIM_WF cell 2.czi';   % Change the file path approporiately here
xpixel       = 0.0313;
output_dir   = './output';
if ~exist(output_dir, 'dir'); mkdir(output_dir); end


%% ====================================================================
%                    STEP 1: LOAD CZI (BOTH CHANNELS)
% =====================================================================
fprintf('\n[1] Loading CZI via ReadImage6D2...\n');
scene       = 1;
MIP_out     = ReadImage6D2(MIP_filename, true, scene);
MIP_image6d = MIP_out{1};

% Expected 6D order from ReadImage6D2: (S, T, C, Z, Y, X)
% With scene=1 and MIP, S=Z=1. Squeeze the singleton S/Z dims only.
sz = size(MIP_image6d);
fprintf('    Raw image6D size: %s\n', mat2str(sz));

% Robust extraction by named dims — keeps the channel index explicit
% so we don't get bitten by squeeze() order.
if numel(sz) == 6
    % (S, T, C, Z, Y, X)
    frames_sim = double(squeeze(MIP_image6d(1, :, 1, p.CH_SIM, :, :)));
    frames_wf  = double(squeeze(MIP_image6d(1, :, 1, p.CH_WF, :, :)));
elseif numel(sz) == 4
    % Pre-squeezed (T, C, Y, X)
    frames_sim = double(squeeze(MIP_image6d(:, p.CH_SIM, :, :)));
    frames_wf  = double(squeeze(MIP_image6d(:, p.CH_WF,  :, :)));
else
    error('Unexpected MIP_image6d dim layout: %s', mat2str(sz));
end

[T, H, W] = size(frames_sim);
assert(isequal(size(frames_sim), size(frames_wf)), 'Channel size mismatch.');
fprintf('    %d frames of %d x %d  (%.3f um/px)\n', T, H, W, xpixel);
fprintf('    SIM (ch%d): mean=%.0f   Widefield (ch%d): mean=%.0f\n', ...
        p.CH_SIM, mean(frames_sim(:)), p.CH_WF, mean(frames_wf(:)));

% Keep original variable name for backward compatibility downstream
frames_brd4 = frames_sim;



%% ====================================================================
%       STEP 2: NUCLEAR SEGMENTATION ON WIDEFIELD + CENTROID TRACKING
% =====================================================================
fprintf('\n[2] Nuclear segmentation on widefield channel...\n');

masks         = false(T, H, W);
areas         = zeros(T, 1);
centroids     = zeros(T, 2);
seg_status    = strings(T, 1);
prev_centroid = [];
prev_mask     = false(H, W);

for t = 1:T
    img_t_wf  = squeeze(frames_wf(t,:,:));   % segmentation input

    % --- Attempt 1: direct segmentation on widefield -------------------
    [mask, cent] = pick_nucleus(img_t_wf, prev_centroid, ...
                                p.NUC_DIAMETER_PX, p.MIN_NUC_AREA, p.MAX_CENT_JUMP);
    status = "direct";

    % --- Attempt 2: average widefield (prev + curr) --------------------
    if isempty(mask) && t > 1
        img_prev_wf = squeeze(frames_wf(t-1,:,:));
        combined    = (img_t_wf + img_prev_wf) / 2;
        [mask, cent] = pick_nucleus(combined, prev_centroid, ...
                                    p.NUC_DIAMETER_PX, p.MIN_NUC_AREA, p.MAX_CENT_JUMP);
        status = "combined";
    end

    % --- Attempt 3: fall back to previous frame's mask -----------------
    if isempty(mask)
        if t == 1
            error('Frame 1 segmentation failed. Check MIN_NUC_AREA or cellpose parameters.');
        end
        mask   = prev_mask;
        cent   = prev_centroid;
        status = "fallback";
    end

    masks(t,:,:)   = mask;
    areas(t)       = sum(mask(:));
    centroids(t,:) = cent;
    seg_status(t)  = status;
    prev_centroid  = cent;
    prev_mask      = mask;

    fprintf('    t=%2d: A=%6d  centroid=(%.1f, %.1f)  [%s]\n', ...
            t, areas(t), cent(1), cent(2), status);

    % Live overlay on the widefield image (the QC that matters)
    imshow(img_t_wf, [], 'InitialMagnification', 'fit');
    hold on;
    visboundaries(mask, 'Color', 'c');
    plot(cent(2), cent(1), 'g+', 'MarkerSize', 14, 'LineWidth', 2);
    hold off;
    title(sprintf('Frame %d  [%s]  — widefield seg input', t, status));
    drawnow;
end

fprintf('    Areas : min=%d  max=%d  median=%d\n', ...
        min(areas), max(areas), round(median(areas)));
fprintf('    Status: direct=%d  combined=%d  fallback=%d\n', ...
        sum(seg_status == "direct"), ...
        sum(seg_status == "combined"), ...
        sum(seg_status == "fallback"));


%% ====================================================================
%          STEP 3: DETECT INTENSITY-DEFECTIVE FRAMES  (SIM CHANNEL)
% =====================================================================
fprintf('\n[3] Detecting defective SIM frames...\n');
metric99 = zeros(T, 1);
for t = 1:T
    img = squeeze(frames_sim(t,:,:));
    m   = squeeze(masks(t,:,:));
    if any(m(:))
        metric99(t) = prctile(img(m), 99);
    end
end

trend   = medfilt1(metric99, 5, 'truncate');
mad_val = median(abs(metric99 - trend)) * 1.4826;
z       = (metric99 - trend) / max(mad_val, 1e-6);
defective = z < p.DEFECT_Z;
fprintf('    Defective frames (z < %.1f): %s\n', p.DEFECT_Z, mat2str(find(defective)'));


%% ====================================================================
%            STEP 4: INTENSITY CORRECTION  (SIM CHANNEL)
% =====================================================================
fprintf('\n[4] Intensity correction of defective frames...\n');
good_idx      = find(~defective);
target_metric = interp1(good_idx, metric99(good_idx), (1:T)', 'linear', 'extrap');
scale_factors = ones(T, 1);
frames_corr   = frames_sim;                     % start from raw SIM

for t = find(defective)'
    scale_factors(t)   = target_metric(t) / max(metric99(t), 1e-6);
    frames_corr(t,:,:) = frames_sim(t,:,:) * scale_factors(t);
    fprintf('    t=%2d: x%.2f  (%.0f -> %.0f)\n', ...
            t, scale_factors(t), metric99(t), target_metric(t));
end


%% ====================================================================
%   STEP 5: PEAK DETECTION ON SIM  (pkfnd + cntrd, 3x mean nuclear thr)
% =====================================================================
fprintf('\n[5] Peak detection on SIM channel...\n');
peaks_per_frame = cell(T, 1);
mean_in         = zeros(T, 1);
n_peaks         = zeros(T, 1);

for t = 1:T
    img = squeeze(frames_corr(t,:,:));
    m   = squeeze(masks(t,:,:));
    [coords, mn] = detect_peaks_in_mask(img, m, p);
    peaks_per_frame{t} = coords;
    mean_in(t)         = mn;
    n_peaks(t)         = size(coords, 1);
end

for t = 1:T
    flag = ''; if defective(t); flag = '  <- defective'; end
    fprintf('    t=%2d: mean_in=%6.1f  thr=%6.1f  n_peaks=%3d%s\n', ...
            t, mean_in(t), p.THR_MULT * mean_in(t), n_peaks(t), flag);
end


%% ====================================================================
%           STEP 6: RIGID-BODY REGISTRATION
% =====================================================================
fprintf('\n[6] Rigid-body registration...\n');

orientations   = zeros(T, 1);
eccentricities = zeros(T, 1);
for t = 1:T
    m  = squeeze(masks(t,:,:));
    rp = regionprops(m, 'Orientation', 'Eccentricity');
    if ~isempty(rp)
        orientations(t)   = -deg2rad(rp(1).Orientation);
        eccentricities(t) = rp(1).Eccentricity;
    end
end

ref_centroid    = median(centroids, 1);
ref_orientation = median(orientations);
fprintf('    Ref centroid    : (%.2f, %.2f)\n', ref_centroid(1), ref_centroid(2));
fprintf('    Ref orientation : %+.2f deg\n', rad2deg(ref_orientation));

centroids_smooth = zeros(size(centroids));
for ax = 1:2
    centroids_smooth(:,ax) = medfilt1(centroids(:,ax), 5, 'truncate');
end
shifts = ref_centroid - centroids_smooth;

dthetas_raw = mod(orientations - ref_orientation + pi/2, pi) - pi/2;
gate = abs(rad2deg(dthetas_raw)) > p.ROT_DEG_GATE & ...
       eccentricities > p.ROT_ECC_GATE;
dthetas       = zeros(T, 1);
dthetas(gate) = dthetas_raw(gate);

peaks_ref = cell(T, 1);
for t = 1:T
    pts = peaks_per_frame{t};
    if isempty(pts); peaks_ref{t} = pts; continue; end
    pts_t = pts + shifts(t,:);
    if abs(dthetas(t)) > 1e-6
        c = cos(-dthetas(t)); s = sin(-dthetas(t));
        R = [c, -s; s, c];
        rel   = pts_t - ref_centroid;
        pts_t = rel * R' + ref_centroid;
    end
    peaks_ref{t} = pts_t;
end


%% ====================================================================
%       STEP 7: LINKING  (Hungarian via matchpairs, memory = 1)
% =====================================================================
fprintf('\n[7] Linking (search_range=%.1f px, memory=%d)...\n', ...
        p.MAX_DISP_PX, p.MEMORY);

tracks = struct('frames', {}, 'coords', {});
active = struct('tid', {}, 'last_coord', {}, 'missed', {});

for i = 1:size(peaks_ref{1}, 1)
    tracks(end+1).frames = 1;                                        
    tracks(end).coords   = peaks_ref{1}(i,:);
    active(end+1).tid        = numel(tracks);                         
    active(end).last_coord   = peaks_ref{1}(i,:);
    active(end).missed       = 0;
end

for t = 2:T
    cur = peaks_ref{t};
    n_a = numel(active);
    n_c = size(cur, 1);
    matched_a = false(n_a, 1);
    matched_c = false(n_c, 1);

    if n_a > 0 && n_c > 0
        A    = cat(1, active.last_coord);
        D    = pdist2(A, cur);
        BIG  = 1e6;
        cost = D;
        cost(D > p.MAX_DISP_PX) = BIG;

        M = matchpairs(cost, p.MAX_DISP_PX + 0.01);
        for k = 1:size(M, 1)
            r = M(k,1); c = M(k,2);
            if cost(r,c) < BIG
                tid = active(r).tid;
                tracks(tid).frames(end+1) = t;
                tracks(tid).coords(end+1,:) = cur(c,:);
                active(r).last_coord = cur(c,:);
                active(r).missed     = 0;
                matched_a(r) = true;
                matched_c(c) = true;
            end
        end
    end

    new_active = struct('tid', {}, 'last_coord', {}, 'missed', {});
    for c = 1:n_c
        if ~matched_c(c)
            tracks(end+1).frames = t;                                 
            tracks(end).coords   = cur(c,:);
            new_active(end+1).tid        = numel(tracks);             
            new_active(end).last_coord   = cur(c,:);
            new_active(end).missed       = 0;
        end
    end

    surviving = struct('tid', {}, 'last_coord', {}, 'missed', {});
    for r = 1:n_a
        if matched_a(r)
            surviving(end+1) = active(r);                             
        else
            active(r).missed = active(r).missed + 1;
            if active(r).missed <= p.MEMORY
                surviving(end+1) = active(r);                          
            end
        end
    end
    active = [surviving, new_active];
end

n_tracks = numel(tracks);
fprintf('    Total tracks: %d  (active at last frame: %d)\n', n_tracks, numel(active));


%% ====================================================================
%   STEP 8: BIRTH / DEATH FRAMES + LIFETIME  (NaN edge-censoring)
% =====================================================================
fprintf('\n[8] Computing birth/death and lifetimes...\n');

first_seen = zeros(n_tracks, 1);
last_seen  = zeros(n_tracks, 1);
duration   = zeros(n_tracks, 1);
n_det      = zeros(n_tracks, 1);

for k = 1:n_tracks
    first_seen(k) = tracks(k).frames(1);
    last_seen(k)  = tracks(k).frames(end);
    duration(k)   = last_seen(k) - first_seen(k) + 1;
    n_det(k)      = numel(tracks(k).frames);
end

duration_sec = duration * p.INTERVAL_SEC;
left_cens    = first_seen == 1;
right_cens   = last_seen  == T;

uncens                     = ~(left_cens | right_cens);
structural_lifetime_frames = duration(uncens);
structural_lifetime_sec    = duration_sec(uncens);

fprintf('    Total tracks               : %d\n', n_tracks);
fprintf('    Uncensored (birth + death) : %d\n', sum(uncens));
fprintf('    Left-censored only         : %d\n', sum(left_cens & ~right_cens));
fprintf('    Right-censored only        : %d\n', sum(~left_cens & right_cens));
fprintf('    Doubly censored            : %d\n', sum(left_cens & right_cens));
if any(uncens)
    fprintf('    Lifetime (sec): med=%.0f  mean=%.0f  max=%d\n', ...
            median(structural_lifetime_sec), mean(structural_lifetime_sec), ...
            max(structural_lifetime_sec));
end


%% ====================================================================
%                  STEP 9: FIGURES  (histogram + KM only)
% =====================================================================
fprintf('\n[9] Generating figures...\n');

dt_s = p.INTERVAL_SEC;

% --- Histogram of structural lifetimes (uncensored only) --------------
if any(uncens)
    fig_hist = figure('Position', [50 50 700 420], 'Name', 'Lifetime histogram');
    edges_sec = (0.5:1:max(structural_lifetime_frames)+0.5) * dt_s;
    histogram(structural_lifetime_sec, edges_sec, ...
              'FaceColor', [0.3 0.5 0.8], 'EdgeColor', 'k'); hold on;
    xline(median(structural_lifetime_sec), 'r--', ...
          sprintf('med = %.0f s', median(structural_lifetime_sec)), ...
          'LineWidth', 1.4);
    xline(mean(structural_lifetime_sec), 'g--', ...
          sprintf('mean = %.0f s', mean(structural_lifetime_sec)), ...
          'LineWidth', 1.4);
    xlabel('Structural lifetime (s)'); ylabel('Count');
    title(sprintf('Structural lifetime distribution  (n = %d uncensored, dt = %d s)', ...
                  numel(structural_lifetime_sec), dt_s));
    box off; set(gca, 'TickDir', 'out');
    exportgraphics(fig_hist, fullfile(output_dir, 'FIG_lifetime_histogram.png'), 'Resolution', 150);
end

% --- Kaplan-Meier (drop left-censored, right-censored = no event) -----
keep_km  = ~left_cens;
d_km_sec = duration_sec(keep_km);
event_km = ~right_cens(keep_km);

if any(event_km)
    [d_km_sorted, ord] = sort(d_km_sec);
    event_sorted = event_km(ord);
    unique_d = unique(d_km_sorted);
    S = 1; times = 0; survs = 1;
    for i = 1:numel(unique_d)
        di = unique_d(i);
        n_events = sum(event_sorted(d_km_sorted == di));
        n_risk   = sum(d_km_sorted >= di);
        if n_risk > 0; S = S * (1 - n_events / n_risk); end
        times(end+1,1) = di;  survs(end+1,1) = S;                      
    end

    fig_km = figure('Position', [50 50 700 480], 'Name', 'Kaplan-Meier');
    stairs(times, survs, 'b-', 'LineWidth', 2.2); hold on; grid on;
    yline(0.5, 'r:', 'LineWidth', 0.8);
    idx_med = find(survs <= 0.5, 1, 'first');
    if ~isempty(idx_med)
        xline(times(idx_med), 'r:', ...
              sprintf('median = %.0f s', times(idx_med)), 'LineWidth', 1.4);
    end
    xlabel('Time (s)'); ylabel('Survival S(t)'); ylim([0 1.05]);
    title(sprintf('Kaplan-Meier  (n = %d, %d events, %d right-censored)', ...
                  numel(d_km_sec), sum(event_km), sum(~event_km)));
    box off; set(gca, 'TickDir', 'out');
    exportgraphics(fig_km, fullfile(output_dir, 'FIG_kaplan_meier.png'), 'Resolution', 150);
end


%% ====================================================================
%                  STEP 10: SAVE .MAT
% =====================================================================
fprintf('\n[10] Saving results...\n');

[~, base_name, ~] = fileparts(MIP_filename);
save_path = fullfile(output_dir, ['brd4_tracking_results_' base_name '.mat']);

results = struct();
results.structural_lifetime_sec    = structural_lifetime_sec;
results.structural_lifetime_frames = structural_lifetime_frames;
results.tracks                     = tracks;
results.first_seen                 = first_seen;
results.last_seen                  = last_seen;
results.left_censored              = left_cens;
results.right_censored             = right_cens;
results.duration_sec               = duration_sec;
results.masks                      = masks;
results.centroids                  = centroids;
results.shifts_yx                  = shifts;
results.dthetas_rad                = dthetas;
results.defective_frames           = find(defective);
results.scale_factors              = scale_factors;
results.peaks_per_frame            = peaks_per_frame;
results.peaks_ref                  = peaks_ref;
results.seg_status                 = seg_status;
results.parameters                 = p;
results.xpixel_um                  = xpixel;
results.source_file                = MIP_filename;

save(save_path, '-struct', 'results', '-v7.3');
fprintf('    Saved to %s\n', save_path);
fprintf('\nDone.\n');


%% ====================================================================
%                       LOCAL FUNCTIONS
% =====================================================================

function [mask, centroid] = pick_nucleus(img, prev_centroid, ...
                                         diameter, min_area, max_jump)
% Run cellpose, filter by area, pick closest to prev centroid.
% Returns [] if nothing qualifies or best match exceeds max_jump.

    mask = []; centroid = [];

    out = cellpose_seg2(img, diameter);
    if isempty(out) || isempty(out{1}); return; end
    ind_masks = out{1};

    cand_masks = {};
    cand_cents = zeros(0, 2);
    cand_areas = [];
    for i = 1:numel(ind_masks)
        m_i = ind_masks{i};
        a_i = sum(m_i(:));
        if a_i < min_area; continue; end
        rp = regionprops(m_i, 'Centroid');
        cand_masks{end+1}    = m_i;                                    
        cand_cents(end+1,:)  = [rp(1).Centroid(2), rp(1).Centroid(1)]; 
        cand_areas(end+1)    = a_i;                                    
    end
    if isempty(cand_masks); return; end

    if isempty(prev_centroid)
        [~, idx] = max(cand_areas);
    else
        d = vecnorm(cand_cents - prev_centroid, 2, 2);
        [d_min, idx] = min(d);
        if d_min > max_jump; return; end
    end

    mask     = cand_masks{idx};
    centroid = cand_cents(idx,:);
end


function [coords, mean_in] = detect_peaks_in_mask(img, mask, p)
% pkfnd + cntrd inside the nuclear mask. Returns Nx2 [y x] sub-pixel.

    coords  = zeros(0, 2);
    mean_in = 0;
    if ~any(mask(:)); return; end

    mean_in  = mean(img(mask));
    thr      = p.THR_MULT * mean_in;
    img_sm   = imgaussfilt(img, p.PEAK_GAUSS_SIGMA);

    pk = pkfnd(img_sm, thr, p.PEAK_SEPARATION);
    if isempty(pk); return; end

    cnt = cntrd(img_sm, pk, p.PEAK_FEATSIZE);
    if isempty(cnt); return; end

    xy = cnt(:, 1:2);
    yi = max(min(round(xy(:,2)), size(mask,1)), 1);
    xi = max(min(round(xy(:,1)), size(mask,2)), 1);
    in_mask = mask(sub2ind(size(mask), yi, xi));
    coords  = [xy(in_mask,2), xy(in_mask,1)];
end