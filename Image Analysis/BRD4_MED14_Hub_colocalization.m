% % % % BRD4_MED14_hub_colcalization_PIPELINE.m

% % %  2-channel SR-SIM timecourse:
% % %    Channel 1 = MED14      labelled with mScarlet
% % %    Channel 2 = BRD4       labelled with HaloTag
% % %
% % %  Pipeline:
% % %    Step 1 : Load CZI (all channels) 
% % %    Step 2 : Nucleus segmentation using Cellpose
% % %    Step 3 : Condensate/Hub peak detection in each nucleus
% % %    Step 4 : MED14 hub to closest BRD4 hub/ BRD4 hub to closest BRD4 hub 3D distance calculations
% % %    Step 5 : Colcalization coefficient (Manders) calculation for individual hubs
% % %    Step 5 : Distance-based colocalized/isolated classification of hubs
% % %    Step 6:  Save .mat
% % %
% % %  Dependencies on path:
% % %  ReadImage6D2, cellpose_seg, pkfnd, cntrd, local_background_2Dg, radial_spot_size2D, gaussian_spot_size2D, feretDiameters


clear; clc; tic;

% %============================== PROJECT DIRECTORY SETUP=============================
here = fileparts(mfilename('fullpath'));   % folder this script lives in
addpath(fullfile(here, 'functions'));
addpath(fullfile(here, 'bioformats'));
javaaddpath(fullfile(here, 'bioformats', 'bioformats_package.jar'));
% %==============================%%==============================%%===================

%% ---- 1) Read Image -----------------------------------------------------
MIP_filename  = '/Volumes/Aniket2/4c MED14 dtag 4-22-2026/JQ1 1um/jq1 UM BRD4 MED14 1hr 3_SIM_Maximum intensity projection.czi';
save_filename = '/Volumes/xiel2lab/Aniket/brd4_med14_coloc_JQ1 4-22-26 image 3_v5.mat';

scene = 1;
MIP_out      = ReadImage6D(MIP_filename, true, scene);
MIP_image6d  = MIP_out{1};

xpixel = 0.0313;           % um/pixel

% ---- Channel assignments ----
channel_MED14 = 1;
channel_BRD4  = 2;

inputImage_MED14 = squeeze(MIP_image6d(1,1,1,channel_MED14,:,:));
inputImage_BRD4  = squeeze(MIP_image6d(1,1,1,channel_BRD4,:,:));

img_MED14 = double(inputImage_MED14);
img_BRD4  = double(inputImage_BRD4);

%% ---- 2) Nucleus segmentation -------------------------------------------
seg_OUT          = cellpose_seg(inputImage_BRD4, 150);
individual_masks = seg_OUT{1};
numNuclei        = numel(individual_masks);
fprintf('Detected %d nuclei.\n', numNuclei);

%% ---- 3) Parameters -----------------------------------------------------
pkfnd_sz        = 15;       % minimum separation between peaks (px)
dilate_margin   = 10;       % DILATE mask to include edge condensates (px)
thresh_mult     = 3;        % threshold = thresh_mult × nuclear mean
sz_cntrd        = 21;       % cntrd fitting window
R_search        = 20;       % radial-profile / background search radius (px)
R_manders       = 5;       % Manders ROI half-width (px)
coloc_radius_nm = 300;      % colocalization proximity cutoff (nm)
coloc_radius_px = coloc_radius_nm / (xpixel * 1000);

%% ========================================================================
%  HELPER: detect + size condensates in one channel
%  ========================================================================
function [out] = detect_and_size(img, individual_masks, numNuclei, ...
        pkfnd_sz, dilate_margin, thresh_mult, sz_cntrd, R_search, xpixel)

    C_in       = [];
    nucleus_id = [];

    for n = 1:numNuclei
        mask_n      = individual_masks{n};
        mask_dilated = imdilate(mask_n, strel('disk', dilate_margin));

        % Threshold from ORIGINAL (undilated) mask mean
        nuc_pixels = img(mask_n > 0);
        if isempty(nuc_pixels), continue; end

        nuc_mean     = mean(nuc_pixels);
        pkfnd_thresh = thresh_mult * nuc_mean;

        % Find peaks in dilated region
        img_masked = img .* mask_dilated;
        pk         = pkfnd(img_masked, pkfnd_thresh, pkfnd_sz);
        if isempty(pk), continue; end

        % Verify each peak is inside dilated mask
        valid = false(size(pk,1),1);
        for j = 1:size(pk,1)
            px = round(pk(j,1)); py = round(pk(j,2));
            if py>=1 && py<=size(mask_dilated,1) && px>=1 && px<=size(mask_dilated,2)
                valid(j) = mask_dilated(py,px) > 0;
            end
        end
        pk = pk(valid,:);
        if ~isempty(pk)
            C_in       = [C_in; pk];
            nucleus_id = [nucleus_id; n*ones(size(pk,1),1)];
            fprintf('  Nucleus %d: mean=%.0f, thresh=%.0f, %d spots\n', ...
                n, nuc_mean, pkfnd_thresh, size(pk,1));
        end
    end

    numPoints = size(C_in,1);
    if numPoints == 0
        out.C_in=[];  out.nucleus_id=[];  out.numPoints=0;
        out.cent=[];  out.FWHM=[];  out.FWHM_r=[];
        out.feret=[];  out.AR=[];  out.AR_feret=[];
        out.bg=[];  out.sigma_r=[];
        return;
    end

    % ---- Refine centroids ------------------------------------------------
    mx        = round(C_in);
    centroids = cntrd(img, mx, sz_cntrd);

    if size(centroids,1) < numPoints
        warning('cntrd rejected %d points near edges.', numPoints-size(centroids,1));
        numPoints = size(centroids,1);
        kept_idx  = zeros(numPoints,1);
        for i = 1:numPoints
            dists = sqrt((C_in(:,1)-centroids(i,1)).^2 + (C_in(:,2)-centroids(i,2)).^2);
            [~, kept_idx(i)] = min(dists);
        end
        C_in       = C_in(kept_idx,:);
        nucleus_id = nucleus_id(kept_idx);
    end

    % ---- Allocate --------------------------------------------------------
    cent      = zeros(numPoints,2);
    FWHM      = zeros(numPoints,1);
    AR        = zeros(numPoints,1);
    AR_feret  = zeros(numPoints,1);
    FWHM_r    = zeros(numPoints,1);
    feret     = zeros(numPoints,1);
    bg        = zeros(numPoints,1);
    sigma_r   = zeros(numPoints,1);
    [m, n_cols] = size(img);

    for i = 1:numPoints
        x0 = centroids(i,1);  y0 = centroids(i,2);
        center = [y0, x0];

        [loc_b, sigma_bg] = local_background_2Dg(img, center, R_search);
        bg(i)      = loc_b;
        sigma_r(i) = sigma_bg;

        x_min = max(1,      round(x0)-R_search);
        x_max = min(n_cols, round(x0)+R_search);
        y_min = max(1,      round(y0)-R_search);
        y_max = min(m,      round(y0)+R_search);
        ROI   = img(y_min:y_max, x_min:x_max) - loc_b;
        ROI(ROI<0) = 0;

        ROI_norm     = ROI / max(ROI(:));
        level        = graythresh(ROI_norm);
        bw           = imbinarize(ROI_norm, level);
        bw_label     = bwlabel(bw);
        center_label = bw_label(round(size(bw,1)/2), round(size(bw,2)/2));
        if center_label > 0, bw = (bw_label == center_label); end

        try
            stats = regionprops(bw, 'MaxFeretDiameter', 'MinFeretDiameter');
            if ~isempty(stats)
                feret(i)    = (stats(1).MaxFeretDiameter + stats(1).MinFeretDiameter)/2 * xpixel*1000;
                AR_feret(i) = stats(1).MaxFeretDiameter / stats(1).MinFeretDiameter;
            else
                feret(i)=NaN; AR_feret(i)=NaN;
            end
        catch
            [maxF, minF] = feretDiameters(bw);
            feret(i)    = (maxF+minF)/2 * xpixel*1000;
            AR_feret(i) = maxF/minF;
        end

        ROI_center = [round(size(ROI,1)/2), round(size(ROI,2)/2)];
        FWHM_r(i)  = radial_spot_size2D(ROI, ROI_center, R_search, xpixel);

        outp    = gaussian_spot_size2D(ROI, x0-x_min+1, y0-y_min+1, xpixel, R_search);
        FWHM(i) = outp(1);
        AR(i)   = outp(2);
        cent(i,1) = outp(3) + x_min - 1;
        cent(i,2) = outp(4) + y_min - 1;
    end

    out.C_in       = C_in;
    out.nucleus_id = nucleus_id;
    out.numPoints  = numPoints;
    out.cent       = cent;
    out.FWHM       = FWHM;
    out.FWHM_r     = FWHM_r;
    out.feret      = feret;
    out.AR         = AR;
    out.AR_feret   = AR_feret;
    out.bg         = bg;
    out.sigma_r    = sigma_r;
end

%% ---- 4) Run detection + sizing on BOTH channels -----------------------
fprintf('\n===== MED14 channel =====\n');
med = detect_and_size(img_MED14, individual_masks, numNuclei, ...
    pkfnd_sz, dilate_margin, thresh_mult, sz_cntrd, R_search, xpixel);

fprintf('\n===== BRD4 channel =====\n');
brd = detect_and_size(img_BRD4, individual_masks, numNuclei, ...
    pkfnd_sz, dilate_margin, thresh_mult, sz_cntrd, R_search, xpixel);

fprintf('\nTotal detected — MED14: %d   BRD4: %d\n', med.numPoints, brd.numPoints);

if med.numPoints == 0 || brd.numPoints == 0
    error('One channel has zero condensates. Adjust thresh_mult or pkfnd_sz.');
end

%% ========================================================================
%  5) COLOCALIZATION — proximity test (global)
%  ========================================================================
D = pdist2(med.cent, brd.cent);                     % nMED14 × nBRD4

[minDist_MED14, nearest_BRD4_idx] = min(D, [], 2);
coloc_MED14_flag  = minDist_MED14 <= coloc_radius_px;
coloc_MED14_count = sum(coloc_MED14_flag);
pct_MED14_coloc   = 100 * coloc_MED14_count / med.numPoints;

[minDist_BRD4, nearest_MED14_idx] = min(D', [], 2);
coloc_BRD4_flag  = minDist_BRD4 <= coloc_radius_px;
coloc_BRD4_count = sum(coloc_BRD4_flag);
pct_BRD4_coloc   = 100 * coloc_BRD4_count / brd.numPoints;

fprintf('\n--- MED14 → BRD4 (%.0f nm): %d/%d (%.1f%%)\n', ...
    coloc_radius_nm, coloc_MED14_count, med.numPoints, pct_MED14_coloc);
fprintf('--- BRD4 → MED14 (%.0f nm): %d/%d (%.1f%%)\n', ...
    coloc_radius_nm, coloc_BRD4_count, brd.numPoints, pct_BRD4_coloc);

dist_MED14_nm = minDist_MED14 * xpixel * 1000;
dist_BRD4_nm  = minDist_BRD4  * xpixel * 1000;

%% ========================================================================
%  6) PER-CONDENSATE MANDERS (background-subtracted)
%  ========================================================================

% ---- Precompute per-nucleus thresholds and means ------------------------
thr_MED14_per_nuc  = NaN(numNuclei,1);
thr_BRD4_per_nuc   = NaN(numNuclei,1);
mean_MED14_per_nuc = NaN(numNuclei,1);
mean_BRD4_per_nuc  = NaN(numNuclei,1);

for n = 1:numNuclei
    mask_n = individual_masks{n};
    pix_m  = img_MED14(mask_n > 0);
    pix_b  = img_BRD4(mask_n > 0);
    if ~isempty(pix_m)
        mean_MED14_per_nuc(n) = mean(pix_m);
        mean_BRD4_per_nuc(n)  = mean(pix_b);
        thr_MED14_per_nuc(n)  = 2*mean_MED14_per_nuc(n);  %%% for Manders calculation use a gentler threshold of just 2* mean nuclear intensity
        thr_BRD4_per_nuc(n)   = 2*mean_BRD4_per_nuc(n);
    end
end

[nRows, nCols] = size(img_MED14);

% ---- M1 per MED14 condensate -------------------------------------------
M1_per_cond = NaN(med.numPoints, 1);

for i = 1:med.numPoints
    nuc_id = med.nucleus_id(i);
    thr_b  = thr_BRD4_per_nuc(nuc_id);
    if isnan(thr_b), continue; end

    x0 = round(med.cent(i,1));  y0 = round(med.cent(i,2));
    r1 = max(1,y0-R_manders);  r2 = min(nRows,y0+R_manders);
    c1 = max(1,x0-R_manders);  c2 = min(nCols,x0+R_manders);

    roi_med  = img_MED14(r1:r2, c1:c2);
    roi_brd  = img_BRD4(r1:r2, c1:c2);
    roi_mask = individual_masks{nuc_id}(r1:r2, c1:c2);
    roi_med  = roi_med .* roi_mask;
    roi_brd  = roi_brd .* roi_mask;

    brd_above = roi_brd > thr_b;
    denom     = sum(roi_med(:));
    if denom > 0
        M1_per_cond(i) = sum(roi_med(brd_above)) / denom;
    end
end

% ---- M2 per BRD4 condensate --------------------------------------------
M2_per_cond = NaN(brd.numPoints, 1);

for i = 1:brd.numPoints
    nuc_id = brd.nucleus_id(i);
    thr_m  = thr_MED14_per_nuc(nuc_id);
    if isnan(thr_m), continue; end

    x0 = round(brd.cent(i,1));  y0 = round(brd.cent(i,2));
    r1 = max(1,y0-R_manders);  r2 = min(nRows,y0+R_manders);
    c1 = max(1,x0-R_manders);  c2 = min(nCols,x0+R_manders);

    roi_brd  = img_BRD4(r1:r2, c1:c2);
    roi_med  = img_MED14(r1:r2, c1:c2);
    roi_mask = individual_masks{nuc_id}(r1:r2, c1:c2);
    roi_brd  = roi_brd .* roi_mask;
    roi_med  = roi_med .* roi_mask;

    med_above = roi_med > thr_m;
    denom     = sum(roi_brd(:));
    if denom > 0
        M2_per_cond(i) = sum(roi_brd(med_above)) / denom;
    end
end

%% ========================================================================
%  7) PER-NUCLEUS: distance classification + Manders by population
%  ========================================================================

% ---- Allocate per-nucleus storage ---------------------------------------
BRD4_count_per_nuc  = zeros(numNuclei, 1);
MED14_count_per_nuc = zeros(numNuclei, 1);

M1_coloc_per_nuc    = NaN(numNuclei, 1);  % mean M1 of colocalized MED14 in each nucleus
M1_isol_per_nuc     = NaN(numNuclei, 1);  % mean M1 of isolated MED14 in each nucleus
M2_coloc_per_nuc    = NaN(numNuclei, 1);  % mean M2 of colocalized BRD4 in each nucleus
M2_isol_per_nuc     = NaN(numNuclei, 1);  % mean M2 of isolated BRD4 in each nucleus

coloc_MED14_count_per_nuc = zeros(numNuclei, 1);
isol_MED14_count_per_nuc  = zeros(numNuclei, 1);
coloc_BRD4_count_per_nuc  = zeros(numNuclei, 1);
isol_BRD4_count_per_nuc   = zeros(numNuclei, 1);

for n = 1:numNuclei

    % Indices of condensates belonging to this nucleus
    idx_m = find(med.nucleus_id == n);
    idx_b = find(brd.nucleus_id == n);

    MED14_count_per_nuc(n) = numel(idx_m);
    BRD4_count_per_nuc(n)  = numel(idx_b);

    % -- MED14 in this nucleus: classify coloc vs isolated -----------------
    if ~isempty(idx_m) && ~isempty(idx_b)
        D_nuc = pdist2(med.cent(idx_m,:), brd.cent(idx_b,:));
        minD_m = min(D_nuc, [], 2);
        is_coloc_m = minD_m <= coloc_radius_px;
    elseif ~isempty(idx_m)
        is_coloc_m = false(numel(idx_m), 1);
    else
        is_coloc_m = [];
    end

    coloc_MED14_count_per_nuc(n) = sum(is_coloc_m);
    isol_MED14_count_per_nuc(n)  = sum(~is_coloc_m);

    if any(is_coloc_m)
        M1_coloc_per_nuc(n) = mean(M1_per_cond(idx_m(is_coloc_m)), 'omitnan');
    end
    if any(~is_coloc_m)
        M1_isol_per_nuc(n) = mean(M1_per_cond(idx_m(~is_coloc_m)), 'omitnan');
    end

    % -- BRD4 in this nucleus: classify coloc vs isolated ------------------
    if ~isempty(idx_b) && ~isempty(idx_m)
        D_nuc_t = pdist2(brd.cent(idx_b,:), med.cent(idx_m,:));
        minD_b = min(D_nuc_t, [], 2);
        is_coloc_b = minD_b <= coloc_radius_px;
    elseif ~isempty(idx_b)
        is_coloc_b = false(numel(idx_b), 1);
    else
        is_coloc_b = [];
    end

    coloc_BRD4_count_per_nuc(n) = sum(is_coloc_b);
    isol_BRD4_count_per_nuc(n)  = sum(~is_coloc_b);

    if any(is_coloc_b)
        M2_coloc_per_nuc(n) = mean(M2_per_cond(idx_b(is_coloc_b)), 'omitnan');
    end
    if any(~is_coloc_b)
        M2_isol_per_nuc(n) = mean(M2_per_cond(idx_b(~is_coloc_b)), 'omitnan');
    end
end

% ---- Global split (all condensates) -------------------------------------
M1_coloc_all = M1_per_cond(coloc_MED14_flag);
M1_isol_all  = M1_per_cond(~coloc_MED14_flag);
M2_coloc_all = M2_per_cond(coloc_BRD4_flag);
M2_isol_all  = M2_per_cond(~coloc_BRD4_flag);

%% ========================================================================
%  8) PRINT SUMMARY
%  ========================================================================
fprintf('\n');
fprintf('=================================================================\n');
fprintf('           MED14 / BRD4 COLOCALIZATION + SIZING SUMMARY\n');
fprintf('=================================================================\n');
fprintf('  Pixel size              : %.4f um/px  (%.1f nm/px)\n', xpixel, xpixel*1000);
fprintf('  Proximity radius        : %.0f nm  (%.1f px)\n', coloc_radius_nm, coloc_radius_px);
fprintf('  Peak threshold          : %d x nuclear mean\n', thresh_mult);
fprintf('  Mask dilation           : %d px\n', dilate_margin);
fprintf('-----------------------------------------------------------------\n');
fprintf('  MED14 condensates       : %d\n', med.numPoints);
fprintf('    FWHM (2D Gauss)       : %.0f +/- %.0f nm\n', mean(med.FWHM), std(med.FWHM));
fprintf('    Feret diameter        : %.0f +/- %.0f nm\n', mean(med.feret,'omitnan'), std(med.feret,'omitnan'));
fprintf('-----------------------------------------------------------------\n');
fprintf('  BRD4  condensates       : %d\n', brd.numPoints);
fprintf('    FWHM (2D Gauss)       : %.0f +/- %.0f nm\n', mean(brd.FWHM), std(brd.FWHM));
fprintf('    Feret diameter        : %.0f +/- %.0f nm\n', mean(brd.feret,'omitnan'), std(brd.feret,'omitnan'));
fprintf('-----------------------------------------------------------------\n');
fprintf('  MED14 with nearby BRD4  : %d / %d  (%.1f%%)\n', coloc_MED14_count, med.numPoints, pct_MED14_coloc);
fprintf('  BRD4  with nearby MED14 : %d / %d  (%.1f%%)\n', coloc_BRD4_count,  brd.numPoints, pct_BRD4_coloc);
fprintf('-----------------------------------------------------------------\n');
fprintf('  Manders M1 coloc MED14  : %.3f +/- %.3f  (n=%d)\n', ...
    mean(M1_coloc_all,'omitnan'), std(M1_coloc_all,'omitnan'), sum(~isnan(M1_coloc_all)));
fprintf('  Manders M1 isol  MED14  : %.3f +/- %.3f  (n=%d)\n', ...
    mean(M1_isol_all,'omitnan'), std(M1_isol_all,'omitnan'), sum(~isnan(M1_isol_all)));
fprintf('  Manders M2 coloc BRD4   : %.3f +/- %.3f  (n=%d)\n', ...
    mean(M2_coloc_all,'omitnan'), std(M2_coloc_all,'omitnan'), sum(~isnan(M2_coloc_all)));
fprintf('  Manders M2 isol  BRD4   : %.3f +/- %.3f  (n=%d)\n', ...
    mean(M2_isol_all,'omitnan'), std(M2_isol_all,'omitnan'), sum(~isnan(M2_isol_all)));
fprintf('=================================================================\n');

% ---- Per-nucleus table --------------------------------------------------
fprintf('\n--- Per-nucleus breakdown ---\n');
fprintf('%-5s %6s %6s %6s %6s %6s %6s %10s %10s %10s %10s\n', ...
    'Nuc', '#MED', '#BRD', 'cMED', 'iMED', 'cBRD', 'iBRD', ...
    'M1_col', 'M1_iso', 'M2_col', 'M2_iso');
fprintf('%s\n', repmat('-',1,95));
for n = 1:numNuclei
    fprintf('%-5d %6d %6d %6d %6d %6d %6d %10.3f %10.3f %10.3f %10.3f\n', ...
        n, MED14_count_per_nuc(n), BRD4_count_per_nuc(n), ...
        coloc_MED14_count_per_nuc(n), isol_MED14_count_per_nuc(n), ...
        coloc_BRD4_count_per_nuc(n), isol_BRD4_count_per_nuc(n), ...
        M1_coloc_per_nuc(n), M1_isol_per_nuc(n), ...
        M2_coloc_per_nuc(n), M2_isol_per_nuc(n));
end

%% ---- 9) VISUALIZATION --------------------------------------------------
figure('Position',[50 50 1600 500]);

% Panel 1 — MED14
subplot(1,3,1);
imshow(uint16(65536 - inputImage_MED14), []);
hold on;
colors = lines(numNuclei);
for n = 1:numNuclei
    idx = (med.nucleus_id == n);
    if any(idx)
        plot(med.cent(idx,1), med.cent(idx,2), '*', 'Color', colors(n,:), 'MarkerSize', 6);
    end
    B = bwboundaries(individual_masks{n});
    for b = 1:numel(B)
        plot(B{b}(:,2), B{b}(:,1), '-', 'Color', colors(n,:), 'LineWidth', 2);
    end
end
title(sprintf('MED14: %d condensates', med.numPoints));

% Panel 2 — BRD4
subplot(1,3,2);
imshow(uint16(65536 - inputImage_BRD4), []);
hold on;
for n = 1:numNuclei
    idx = (brd.nucleus_id == n);
    if any(idx)
        plot(brd.cent(idx,1), brd.cent(idx,2), '*', 'Color', colors(n,:), 'MarkerSize', 6);
    end
    B = bwboundaries(individual_masks{n});
    for b = 1:numel(B)
        plot(B{b}(:,2), B{b}(:,1), '-', 'Color', colors(n,:), 'LineWidth', 2);
    end
end
title(sprintf('BRD4: %d condensates', brd.numPoints));

% Panel 3 — Colocalization overlay
subplot(1,3,3);
overlay = cat(3, mat2gray(img_BRD4), mat2gray(img_MED14), mat2gray(img_BRD4));
imshow(overlay);
hold on;
plot(med.cent(coloc_MED14_flag,1),  med.cent(coloc_MED14_flag,2), ...
    'wo', 'MarkerSize', 10, 'LineWidth', 1.5);
plot(med.cent(~coloc_MED14_flag,1), med.cent(~coloc_MED14_flag,2), ...
    'go', 'MarkerSize', 5, 'LineWidth', 0.5);
plot(brd.cent(coloc_BRD4_flag,1),  brd.cent(coloc_BRD4_flag,2), ...
    'w+', 'MarkerSize', 10, 'LineWidth', 1.5);
plot(brd.cent(~coloc_BRD4_flag,1), brd.cent(~coloc_BRD4_flag,2), ...
    'm+', 'MarkerSize', 5, 'LineWidth', 0.5);
title(sprintf('Coloc: MED14 %.0f%%, BRD4 %.0f%%', pct_MED14_coloc, pct_BRD4_coloc));
legend('MED14 coloc','MED14 solo','BRD4 coloc','BRD4 solo', ...
    'Location','southeast','TextColor','w');
sgtitle('MED14 / BRD4 Condensate Colocalization', 'FontSize', 14, 'FontWeight', 'bold');

%% ---- 10) Distance histograms ------------------------------------------
coloc_thresh_nm = coloc_radius_nm;

frac_MED14_left  = sum(dist_MED14_nm <= coloc_thresh_nm) / numel(dist_MED14_nm);
frac_MED14_right = 1 - frac_MED14_left;
frac_BRD4_left   = sum(dist_BRD4_nm  <= coloc_thresh_nm) / numel(dist_BRD4_nm);
frac_BRD4_right  = 1 - frac_BRD4_left;

fprintf('\n--- Area partition at %.0f nm ---\n', coloc_thresh_nm);
fprintf('  MED14→BRD4 :  left = %.3f  |  right = %.3f\n', frac_MED14_left, frac_MED14_right);
fprintf('  BRD4→MED14 :  left = %.3f  |  right = %.3f\n', frac_BRD4_left,  frac_BRD4_right);

figure('Position',[50 600 800 300]);
edges = 0:15:1500;

subplot(1,2,1);
histogram(dist_MED14_nm, edges, 'Normalization','probability', 'FaceColor',[0.2 0.8 0.2]);
hold on; xline(coloc_thresh_nm, 'r--', 'LineWidth', 1.5);
xlim([0 1500]); xlabel('Distance to nearest BRD4 (nm)'); ylabel('Probability');
title('MED14 \rightarrow nearest BRD4');

subplot(1,2,2);
histogram(dist_BRD4_nm, edges, 'Normalization','probability', 'FaceColor',[0.8 0.2 0.8]);
hold on; xline(coloc_thresh_nm, 'r--', 'LineWidth', 1.5);
xlim([0 1500]); xlabel('Distance to nearest MED14 (nm)'); ylabel('Probability');
title('BRD4 \rightarrow nearest MED14');

set(findall(gcf,'-property','FontSize'), 'FontSize', 14, 'FontName', 'Arial');

%% ---- 11) Manders histograms: coloc vs isolated -------------------------
figure('Position',[50 700 900 350]);

subplot(1,2,1);
histogram(M1_coloc_all, 0:0.05:1, 'Normalization','probability', ...
    'FaceColor',[0.2 0.8 0.2], 'FaceAlpha',0.8);
hold on;
histogram(M1_isol_all, 0:0.05:1, 'Normalization','probability', ...
    'FaceColor',[0.6 0.6 0.6], 'FaceAlpha',0.5);
xline(mean(M1_coloc_all,'omitnan'), 'g--', 'LineWidth', 1.5);
xline(mean(M1_isol_all,'omitnan'),  'k--', 'LineWidth', 1.5);
xlim([0 1]); xlabel('M1 (Manders)'); ylabel('Probability');
title('MED14 condensates');
legend(sprintf('Coloc (%.3f)', mean(M1_coloc_all,'omitnan')), ...
       sprintf('Isolated (%.3f)', mean(M1_isol_all,'omitnan')), 'Location','northwest');

subplot(1,2,2);
histogram(M2_coloc_all, 0:0.05:1, 'Normalization','probability', ...
    'FaceColor',[0.8 0.2 0.8], 'FaceAlpha',0.8);
hold on;
histogram(M2_isol_all, 0:0.05:1, 'Normalization','probability', ...
    'FaceColor',[0.6 0.6 0.6], 'FaceAlpha',0.5);
xline(mean(M2_coloc_all,'omitnan'), 'm--', 'LineWidth', 1.5);
xline(mean(M2_isol_all,'omitnan'),  'k--', 'LineWidth', 1.5);
xlim([0 1]); xlabel('M2 (Manders)'); ylabel('Probability');
title('BRD4 condensates');
legend(sprintf('Coloc (%.3f)', mean(M2_coloc_all,'omitnan')), ...
       sprintf('Isolated (%.3f)', mean(M2_isol_all,'omitnan')), 'Location','northwest');

sgtitle(sprintf('Manders: Colocalized vs Isolated (%.0f nm cutoff)', coloc_radius_nm), ...
    'FontSize', 14, 'FontWeight','bold');
set(findall(gcf,'-property','FontSize'), 'FontSize', 14, 'FontName', 'Arial');

%% ========================================================================
%  12) SAVE — numbered to match your list
%  ========================================================================

results = struct();

% 1) BRD4 condensate count per nucleus
results.BRD4_count_per_nuc         = BRD4_count_per_nuc;

% 2) MED14 condensate count per nucleus
results.MED14_count_per_nuc        = MED14_count_per_nuc;

% 3) BRD4 Feret diameter — all condensates in image
results.BRD4_feret_all_nm          = brd.feret;

% 4) MED14 Feret diameter — all condensates in image
results.MED14_feret_all_nm         = med.feret;

% 5) Pairwise distances: all BRD4 → nearest MED14 (nm)
results.dist_BRD4_to_MED14_nm     = dist_BRD4_nm;

% 6) Pairwise distances: all MED14 → nearest BRD4 (nm)
results.dist_MED14_to_BRD4_nm     = dist_MED14_nm;

% 7) Manders M2 of all colocalized BRD4 in image
results.M2_coloc_BRD4_all         = M2_coloc_all;

% 8) Manders M2 of colocalized BRD4 — averaged per nucleus
results.M2_coloc_BRD4_per_nuc     = M2_coloc_per_nuc;

% 9) Manders M2 of all isolated BRD4 in image
results.M2_isol_BRD4_all          = M2_isol_all;

% 10) Manders M2 of isolated BRD4 — averaged per nucleus
results.M2_isol_BRD4_per_nuc      = M2_isol_per_nuc;

% 11) Manders M1 of all colocalized MED14 in image
results.M1_coloc_MED14_all        = M1_coloc_all;

% 12) Manders M1 of colocalized MED14 — averaged per nucleus
results.M1_coloc_MED14_per_nuc    = M1_coloc_per_nuc;

% 13) Manders M1 of all isolated MED14 in image
results.M1_isol_MED14_all         = M1_isol_all;

% 14) Manders M1 of isolated MED14 — averaged per nucleus
results.M1_isol_MED14_per_nuc     = M1_isol_per_nuc;

% 15) Colocalized fraction of BRD4 per nucleus
results.BRD4_coloc_frac_per_nuc    = coloc_BRD4_count_per_nuc ./ BRD4_count_per_nuc;

% 16) Isolated fraction of BRD4 per nucleus
results.BRD4_isol_frac_per_nuc     = isol_BRD4_count_per_nuc ./ BRD4_count_per_nuc;

% 17) Colocalized fraction of MED14 per nucleus
results.MED14_coloc_frac_per_nuc   = coloc_MED14_count_per_nuc ./ MED14_count_per_nuc;

% 18) Isolated fraction of MED14 per nucleus
results.MED14_isol_frac_per_nuc    = isol_MED14_count_per_nuc ./ MED14_count_per_nuc;

% ---- Extras (full structs for downstream use) ----
results.med                        = med;
results.brd                        = brd;
results.coloc_radius_nm            = coloc_radius_nm;
results.numNuclei                  = numNuclei;

save(save_filename, 'results');
fprintf('\nSaved to %s\n', save_filename);
toc