% % % % BRD4_MED14_hub_colcalization_PIPELINE.m

% % %  4-channel SR-SIM timecourse:
% % %    Channel 1 = MS2       used for transcription state(bursting/inactive)
% % %    Channel 2 = Enhancer  used for enhancer 3D location
% % %    Channel 3 = MED14     labelled with mScarlet
% % %    Channel 4 = BRD4      labelled with HaloTag
% % %
% % %  Pipeline:
% % %    Step 1 : Load CZI (all channels) 
% % %    Step 2 : Nucleus segmentation using Cellpose
% % %    Step 3 : Refine enhancer 3D position via 3D Gaussian fit
% % %    Step 4 : MS2 spot refined 3D position starting from enhancer co-ordinates and spot intensity extraction
% % %    Step 5 : Find top-3 (primary/secondary/tertiary) nearest BRD4 condensates/hubs in each nucleus
% % %    Step 6 : MED14 hub to closest BRD4 hub/ BRD4 hub to closest BRD4 hub 3D distance calculations
% % %    Step 7 : For each BRD4 condensate/hub, do a focused local MED14 search and flag colocalization (yes/no)
% % %    Step 8 : Save .mat
% % %
% % %  Dependencies on path:
% % %  ReadImage6D2, cellpose_seg, point_location2D, single_zstack, fit_Gaussian3D, condensate_search_v6, pkfnd, cntrd, imagemask, maskavg

clear; clc; tic;

% %============================== PROJECT DIRECTORY SETUP=============================
here = fileparts(mfilename('fullpath'));   % folder this script lives in
addpath(fullfile(here, 'functions'));
addpath(fullfile(here, 'bioformats'));
javaaddpath(fullfile(here, 'bioformats', 'bioformats_package.jar'));
% %==============================%%==============================%%===================

%% ---- 1) File paths and key parameters ----------------------------------
Input_zstack  = '/Volumes/Aniket2/BRD4 med14 4c 5-7-2026/BRD4 MED14 ctrl 1_SIM.czi';
MIP_filename  = '/Volumes/Aniket2/BRD4 med14 4c 5-7-2026/BRD4 MED14 ctrl 1_SIM_Maximum intensity projection.czi';
save_filename = 'BRD4 MED14 5-7-26 ctrl 1_4color v1.mat';

scene = 1;

% --- Channel assignments ---
channel_MS2   = 1;
channel_enh   = 2;
channel_MED14 = 3;
channel_BRD4  = 4;

% --- Approximate enhancer locations (x_pixel, y_pixel) ---
C_in = [
969.1054313	2335.207668
816.5814696	1011.277955
1801.341853	535.2076677
1893.674121	1005.591054
1839.680511	1471.629393
700.5750799	330.1916933
590.543131	1878.690096
1201.629393	1251.629393 ];     % Approximate pixel co-ordinates of enhancer (Verify each cell has all 4 channel labels in ImageJ)


% --- Imaging / pixel parameters ---
xpixel = 0.0313;     % um/pixel
ypixel = 0.0313;

% --- Detection / fitting parameters ---
R_fit               = 8;     % 3D Gaussian fitting radius (px)
mult                = 3;      % BRD4/MED14 detection threshold multiplier (x nuc mean)
R_search_c          = 10;     % BRD4 condensate initial search radius (auto-expands)
R_manders           = 5;      % Manders ROI half-width (px)
coloc_radius_nm     = 300;    % BRD4-MED14 proximity cutoff (nm)
dilate_margin       = 10;     % nucleus dilation for BRD4 condensate search; 0 = off
manders_thresh_mult = 2;      % gentler threshold (x nuc mean) for Manders denominator
fallback_radius_px  = 300;    % synthetic-ROI radius when locus is outside all
                              
m_enh = 0.1;

%% ---- 2) Read images ----------------------------------------------------
MIP_out      = ReadImage6D2(MIP_filename, true, scene);
stack_out    = ReadImage6D2(Input_zstack, true, scene);
metadata     = stack_out{2};
full_stack   = stack_out{1};
MIP_image6d  = MIP_out{1};

zpixel  = metadata.ScaleZ;
spacing = [xpixel ypixel zpixel];
z_slice = metadata.SizeZ;


% MIP per channel
img_MS2_mip   = squeeze(MIP_image6d(1,1,1,channel_MS2,:,:));
img_enh_mip   = squeeze(MIP_image6d(1,1,1,channel_enh,:,:));
img_MED14_mip = squeeze(MIP_image6d(1,1,1,channel_MED14,:,:));
img_BRD4_mip  = squeeze(MIP_image6d(1,1,1,channel_BRD4,:,:));

img_MED14_d = double(img_MED14_mip);
img_BRD4_d  = double(img_BRD4_mip);

I_cond = img_MED14_d + img_BRD4_d;

[nRows, nCols] = size(img_BRD4_mip);

vis             = size(C_in,1);
coloc_radius_px = coloc_radius_nm / (xpixel * 1000);

% Pre-extract z-stacks per channel (single timepoint k=1)
stack_enh   = single_zstack(full_stack, 1, 1, channel_enh);
stack_MS2   = single_zstack(full_stack, 1, 1, channel_MS2);
stack_MED14 = single_zstack(full_stack, 1, 1, channel_MED14);
stack_BRD4  = single_zstack(full_stack, 1, 1, channel_BRD4);

%% ---- 3) Allocate output structures -------------------------------------
% Enhancer
enh_xyz  = zeros(vis, 3);
C_cent_e = zeros(vis, 3);

% MS2 / TSS
TSS_xyz   = NaN(vis, 3);
MS2_score = zeros(vis, 2);     % col 1: spot flag (1 yes, -1 no); col 2: intensity

% BRD4 (3 condensates per locus)
BRD4_xyz  = NaN(vis, 3, 3);    % vis x xyz x {1,2,3}
BRD4_size = NaN(vis, 3);
BRD4_int  = NaN(vis, 3);
BRD4_AR   = NaN(vis, 3);

% MED14 colocalization per BRD4 condensate
MED14_nearest_xyz  = NaN(vis, 3, 3);
MED14_nearest_d3D  = NaN(vis, 3);     % nm
MED14_nearest_d2D  = NaN(vis, 3);     % nm
MED14_coloc_flag   = false(vis, 3);
M2_BRD4_per_cond   = NaN(vis, 3);     % Manders M2 (BRD4 -> MED14)

% Distances (nm)
E_TSS_dist    = NaN(vis, 2);          % col1 3D, col2 2D
E_BRD4_dist   = NaN(vis, 2, 3);
TSS_BRD4_dist = NaN(vis, 2, 3);

% Bookkeeping
% host_nuc_id semantics:  >0 = real Cellpose nucleus id
%                          0 = fallback (synthetic ROI; no real host)
%                        NaN = never reached (shouldn't occur)
host_nuc_id = NaN(vis, 1);

%% ---- 4) Cellpose segmentation (on BRD4 channel) ------------------------
% seg_OUT          = cellpose_seg(img_BRD4_mip, 150);
seg_OUT          = cellpose_seg(I_cond, 150);
individual_masks = seg_OUT{1};
nuc_label_image  = seg_OUT{1,3};
numNuclei        = numel(individual_masks);
fprintf('Detected %d nuclei.\n', numNuclei);

% Per-nucleus mean intensities for BRD4 and MED14
nuc_mean_BRD4  = NaN(numNuclei, 1);
nuc_mean_MED14 = NaN(numNuclei, 1);
for n = 1:numNuclei
    mask_n = individual_masks{n};
    pix_b  = double(img_BRD4_mip(mask_n > 0));
    pix_m  = double(img_MED14_mip(mask_n > 0));
    if ~isempty(pix_b)
        nuc_mean_BRD4(n)  = mean(pix_b);
        nuc_mean_MED14(n) = mean(pix_m);
    end
end

% --- Global fallback means: used when a locus is outside all nuclei -----
all_nuc_mask = false(nRows, nCols);
for n = 1:numNuclei
    all_nuc_mask = all_nuc_mask | individual_masks{n};
end
if any(all_nuc_mask(:))
    global_mean_BRD4  = mean(img_BRD4_d(all_nuc_mask));
    global_mean_MED14 = mean(img_MED14_d(all_nuc_mask));
else
    global_mean_BRD4  = mean(img_BRD4_d(:));
    global_mean_MED14 = mean(img_MED14_d(:));
end
fprintf('Global fallback means -> BRD4: %.1f | MED14: %.1f\n', ...
        global_mean_BRD4, global_mean_MED14);

% --- Visualize Cellpose segmentation ---
figure;
imshow(uint16(65536 - img_BRD4_mip), []);
hold on;
colors = lines(numNuclei);
for n = 1:numNuclei
    B = bwboundaries(individual_masks{n});
    for b = 1:numel(B)
        plot(B{b}(:,2), B{b}(:,1), '-', 'Color', colors(n,:), 'LineWidth', 1.5);
    end
    props = regionprops(individual_masks{n}, 'Centroid');
    if ~isempty(props)
        text(props(1).Centroid(1), props(1).Centroid(2), sprintf('%d', n), ...
             'Color', colors(n,:), 'FontSize', 12, 'FontWeight', 'bold', ...
             'HorizontalAlignment', 'center');
    end
end
for j = 1:vis
    plot(C_in(j,1), C_in(j,2), 'r+', 'MarkerSize', 10, 'LineWidth', 1.5);
    text(C_in(j,1)+15, C_in(j,2), sprintf('e%d', j), 'Color', 'r', 'FontSize', 10);
end
hold off;
title(sprintf('Cellpose segmentation: %d nuclei (BRD4 channel) + enhancer seeds', numNuclei));

%% ---- 5) Refine enhancer 3D positions -----------------------------------
for j = 1:vis
    % guess_e = [C_in(j,1), C_in(j,2), round(z_slice/2)];
    [out_1, mip_c1] = spotPosition3D(img_enh_mip, stack_enh, C_in(j,1), C_in(j,2), z_slice, R_fit, m_enh, spacing);

    enh_xyz(j, :) = [out_1(1), out_1(2), out_1(3)];
    % [x_ref_e, y_ref_e, z_ref_e, ~] = fit_Gaussian3D(stack_enh, guess_e, R_fit, spacing);

    x_ref_e = enh_xyz(j,1)/xpixel;
    y_ref_e  = enh_xyz(j,2)/ypixel;
    z_ref_e = enh_xyz(j,3)/zpixel + 1;
    C_cent_e(j,:) = round([x_ref_e, y_ref_e, z_ref_e]);
end

% --- visualize refined enhancer positions ---
positions_enh = zeros(vis, 3);
labels_enh    = cell(1, vis);
for j = 1:vis
    positions_enh(j,:) = [C_cent_e(j,1), C_cent_e(j,2), 5];
    labels_enh{j}      = sprintf('e%d', j);
end
figure;
RGB_enh = insertObjectAnnotation(uint16(img_enh_mip), 'circle', positions_enh, labels_enh);
imshow(RGB_enh, []);
title('Refined enhancer locations');

%% ---- 6) MS2 / TSS extraction -------------------------------------------
for j = 1:vis
    masked_img_MS2 = imagemask(img_MS2_mip, C_cent_e(j,1), C_cent_e(j,2), 30);
    t_MS2          = 0.05 * max(img_MS2_mip(:));
    pk_MS2         = pkfnd(masked_img_MS2, t_MS2, 5);

    if size(pk_MS2,1) >= 1
        I_spots = zeros(size(pk_MS2,1), 1);
        for i = 1:size(pk_MS2,1)
            I_spots(i) = img_MS2_mip(pk_MS2(i,2), pk_MS2(i,1));
        end
        [~, idx_brightest] = max(I_spots);
        TSS_x = pk_MS2(idx_brightest, 1);
        TSS_y = pk_MS2(idx_brightest, 2);
        z_m   = enh_xyz(j,3) / zpixel + 1;

        guess_m = [TSS_x, TSS_y, z_m];
        [x_ref_m, y_ref_m, z_ref_m, ~] = fit_Gaussian3D(stack_MS2, guess_m, R_fit, spacing);

        TSS_xyz(j,1) = x_ref_m * xpixel;
        TSS_xyz(j,2) = y_ref_m * ypixel;
        TSS_xyz(j,3) = (z_ref_m - 1) * zpixel;

        MS2_score(j,1) = 1;
        MS2_score(j,2) = maskavg(img_MS2_mip, round(x_ref_m), round(y_ref_m), 5);
    else
        TSS_xyz(j,:)   = NaN;
        MS2_score(j,1) = -1;
        MS2_score(j,2) = maskavg(img_MS2_mip, C_cent_e(j,1), C_cent_e(j,2), 5);
    end
end

%% ---- 7) Per-locus: host nucleus -> top-3 BRD4 -> focused MED14 search --
all_BRD4_positions = [];
all_BRD4_labels    = {};

R_coloc_px = ceil(coloc_radius_nm / (xpixel * 1000));   % ~10 px for 300 nm

for j = 1:vis
    x_loc = C_cent_e(j,1);
    y_loc = C_cent_e(j,2);

    nuc_id     = point_location2D(nuc_label_image, x_loc, y_loc);
    in_nucleus = ~(isnan(nuc_id) || nuc_id == 0 || nuc_id > numNuclei);

    if in_nucleus
        host_nuc_id(j) = nuc_id;
        nuc_mask       = individual_masks{nuc_id};
        t_cond_BRD4    = nuc_mean_BRD4(nuc_id);
        t_cond_MED14   = nuc_mean_MED14(nuc_id);
    else
        fprintf(['Locus %d: outside all Cellpose nuclei - ' ...
                 'using global fallback (synthetic ROI around locus)\n'], j);
        host_nuc_id(j) = 0;    % sentinel: 0 = fallback, >0 = real nucleus
        [Xg, Yg]       = meshgrid(1:nCols, 1:nRows);
        nuc_mask       = (Xg - x_loc).^2 + (Yg - y_loc).^2 <= fallback_radius_px^2;
        t_cond_BRD4    = global_mean_BRD4;
        t_cond_MED14   = global_mean_MED14;
    end

    % --- BRD4 top-3 within host nucleus (or fallback ROI) ---------------
    if dilate_margin > 0
        mask_use = imdilate(nuc_mask, strel('disk', dilate_margin));
    else
        mask_use = nuc_mask;
    end

    masked_BRD4 = img_BRD4_d .* double(mask_use);

    BW  = imbinarize(masked_BRD4, mult * t_cond_BRD4 / 2);
    BW2 = bwareafilt(BW, [10 inf]);
    B3  = masked_BRD4 .* BW2;

    enh_pos = enh_xyz(j,:);

    OUT = condensate_search_v6(stack_BRD4, B3, enh_pos, spacing, ...
                               t_cond_BRD4, mult, R_search_c);
    % auto-expand search radius if first attempt fails
    for k = 2:20
        if isnan(OUT{1}(1,1))
            OUT = condensate_search_v6(stack_BRD4, B3, enh_pos, spacing, ...
                                       t_cond_BRD4, mult, k * R_search_c);
        end
    end

    for c = 1:3
        BRD4_xyz(j,1,c) = OUT{c}(1,1);
        BRD4_xyz(j,2,c) = OUT{c}(1,2);
        BRD4_xyz(j,3,c) = OUT{c}(1,3);
        BRD4_size(j,c)  = OUT{c}(1,4);
        BRD4_int(j,c)   = OUT{c}(1,5);
        BRD4_AR(j,c)    = OUT{c}(1,6);
    end

    % --- MED14 colocalization for each BRD4 (focused local search) ------
    pkfnd_thresh_med = mult * t_cond_MED14;
    thr_med14        = manders_thresh_mult * t_cond_MED14;

    for c = 1:3
        if isnan(BRD4_xyz(j,1,c)), continue; end

        x0 = round(BRD4_xyz(j,1,c) / xpixel);
        y0 = round(BRD4_xyz(j,2,c) / ypixel);
        if x0 < 1 || x0 > nCols || y0 < 1 || y0 > nRows, continue; end

        % --- Local MED14 search window around the BRD4 centroid ---------
        r1 = max(1, y0 - R_coloc_px);  r2 = min(nRows, y0 + R_coloc_px);
        c1 = max(1, x0 - R_coloc_px);  c2 = min(nCols, x0 + R_coloc_px);

        local_med14 = img_MED14_d(r1:r2, c1:c2) .* double(nuc_mask(r1:r2, c1:c2));
        pk_med      = pkfnd(local_med14, pkfnd_thresh_med, 5);

        if ~isempty(pk_med)
            % convert local px coords -> global px coords
            pk_global_x = pk_med(:,1) + c1 - 1;
            pk_global_y = pk_med(:,2) + r1 - 1;
            nCand       = numel(pk_global_x);

            cand_xyz = NaN(nCand, 3);
            cand_d3D = NaN(nCand, 1);
            for q = 1:nCand
                z_guess = round(BRD4_xyz(j,3,c) / zpixel) + 1;
                guess_m = [pk_global_x(q), pk_global_y(q), z_guess];
                try
                    [x_r, y_r, z_r, ~] = fit_Gaussian3D(stack_MED14, guess_m, R_fit, spacing);
                    if ~isnan(x_r)
                        cand_xyz(q,:) = [x_r*xpixel, y_r*ypixel, (z_r-1)*zpixel];
                    else
                        cand_xyz(q,:) = [pk_global_x(q)*xpixel, pk_global_y(q)*ypixel, ...
                                         BRD4_xyz(j,3,c)];
                    end
                catch
                    cand_xyz(q,:) = [pk_global_x(q)*xpixel, pk_global_y(q)*ypixel, ...
                                     BRD4_xyz(j,3,c)];
                end

                d3 = cand_xyz(q,:) - [BRD4_xyz(j,1,c) BRD4_xyz(j,2,c) BRD4_xyz(j,3,c)];
                cand_d3D(q) = norm(d3) * 1000;
            end

            [min_d3D, idx_min] = min(cand_d3D);
            d2 = cand_xyz(idx_min,1:2) - [BRD4_xyz(j,1,c) BRD4_xyz(j,2,c)];

            MED14_nearest_xyz(j,:,c) = cand_xyz(idx_min,:);
            MED14_nearest_d3D(j,c)   = min_d3D;
            MED14_nearest_d2D(j,c)   = norm(d2) * 1000;
            MED14_coloc_flag(j,c)    = (min_d3D <= coloc_radius_nm);
        else
            % no MED14 peak within the search window -> not colocalized
            MED14_coloc_flag(j,c) = false;
            % MED14_nearest_xyz / d3D / d2D stay as NaN
        end

        % --- Manders M2 (BRD4 -> MED14) in 5-px ROI around BRD4 ---------
        rr1 = max(1, y0 - R_manders);  rr2 = min(nRows, y0 + R_manders);
        cc1 = max(1, x0 - R_manders);  cc2 = min(nCols, x0 + R_manders);

        roi_brd  = img_BRD4_d(rr1:rr2, cc1:cc2);
        roi_med  = img_MED14_d(rr1:rr2, cc1:cc2);
        roi_mask = double(nuc_mask(rr1:rr2, cc1:cc2));     % UN-dilated mask
        roi_brd  = roi_brd .* roi_mask;
        roi_med  = roi_med .* roi_mask;

        med_above = roi_med > thr_med14;
        denom     = sum(roi_brd(:));
        if denom > 0
            M2_BRD4_per_cond(j,c) = sum(roi_brd(med_above)) / denom;
        end

        all_BRD4_positions = [all_BRD4_positions; ...
            BRD4_xyz(j,1,c)/xpixel  BRD4_xyz(j,2,c)/ypixel  5]; %#ok<AGROW>
        all_BRD4_labels{end+1} = sprintf('e%d-B%d', j, c); %#ok<SAGROW>
    end
end

% --- visualize all top-3 BRD4 condensates ---
if ~isempty(all_BRD4_positions)
    figure;
    I16bit2  = uint16(65536 - img_BRD4_mip);
    RGB_BRD4 = insertObjectAnnotation(I16bit2, 'circle', ...
                                      all_BRD4_positions, all_BRD4_labels);
    imshow(RGB_BRD4, []);
    title('Top-3 BRD4 condensates per locus (1/2/3)');
end

%% ---- 8) Compute all distances ------------------------------------------
for j = 1:vis
    if MS2_score(j,1) == 1
        d3 = enh_xyz(j,:) - TSS_xyz(j,:);
        E_TSS_dist(j,1) = norm(d3)        * 1000;
        E_TSS_dist(j,2) = norm(d3(1:2))   * 1000;
    end

    for c = 1:3
        if isnan(BRD4_xyz(j,1,c)), continue; end
        B = [BRD4_xyz(j,1,c)  BRD4_xyz(j,2,c)  BRD4_xyz(j,3,c)];

        d3 = enh_xyz(j,:) - B;
        E_BRD4_dist(j,1,c) = norm(d3)      * 1000;
        E_BRD4_dist(j,2,c) = norm(d3(1:2)) * 1000;

        if MS2_score(j,1) == 1
            d3 = TSS_xyz(j,:) - B;
            TSS_BRD4_dist(j,1,c) = norm(d3)      * 1000;
            TSS_BRD4_dist(j,2,c) = norm(d3(1:2)) * 1000;
        end
    end
end

%% ---- 9) Save -----------------------------------------------------------
results = struct();

results.C_in              = C_in;
results.host_nuc_id       = host_nuc_id;        % >0 real nuc, 0 fallback
results.numNuclei         = numNuclei;

results.enh_xyz           = enh_xyz;            % vis x 3 (um)
results.TSS_xyz           = TSS_xyz;            % vis x 3 (um), NaN if no spot
results.MS2_score         = MS2_score;          % col1 spot flag, col2 intensity

results.BRD4_xyz          = BRD4_xyz;           % vis x xyz x {1,2,3}, um
results.BRD4_size         = BRD4_size;
results.BRD4_int          = BRD4_int;
results.BRD4_AR           = BRD4_AR;

results.MED14_nearest_xyz = MED14_nearest_xyz;  % vis x xyz x {1,2,3}
results.MED14_nearest_d3D = MED14_nearest_d3D;  % nm (NaN if no peak found)
results.MED14_nearest_d2D = MED14_nearest_d2D;  % nm
results.MED14_coloc_flag  = MED14_coloc_flag;   % logical, vis x 3
results.M2_BRD4_per_cond  = M2_BRD4_per_cond;

results.E_TSS_dist        = E_TSS_dist;         % vis x 2 (3D, 2D) nm
results.E_BRD4_dist       = E_BRD4_dist;        % vis x 2 x 3
results.TSS_BRD4_dist     = TSS_BRD4_dist;      % vis x 2 x 3

results.nuc_mean_BRD4     = nuc_mean_BRD4;
results.nuc_mean_MED14    = nuc_mean_MED14;
results.global_mean_BRD4  = global_mean_BRD4;
results.global_mean_MED14 = global_mean_MED14;

results.params.xpixel              = xpixel;
results.params.ypixel              = ypixel;
results.params.zpixel              = zpixel;
results.params.mult                = mult;
results.params.coloc_radius_nm     = coloc_radius_nm;
results.params.R_manders           = R_manders;
results.params.R_fit               = R_fit;
results.params.R_search_c          = R_search_c;
results.params.dilate_margin       = dilate_margin;
results.params.manders_thresh_mult = manders_thresh_mult;
results.params.fallback_radius_px  = fallback_radius_px;

save(save_filename, 'results');
fprintf('Saved to %s\n', save_filename);
toc