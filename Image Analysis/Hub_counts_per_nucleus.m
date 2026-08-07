% % % % Hub_count_per nucleus_calculation_PIPELINE.m

% % %  Single/Multi-channel SR-SIM images: Batch processing of multiple images
% % %
% % %  Pipeline:
% % %    Step 1 : Load individual image CZI (all channels) 
% % %    Step 2 : Nucleus segmentation using Cellpose
% % %    Step 3 : Filter incomplete nuclei at the edges
% % %    Step 4 : Puncta detection using DBSCAN
% % %    Step 5 : Pool all hub counts and save .mat
% % %
% % %  Dependencies on path:
% % %  ReadImage6D2, cellpose_seg, pkfnd, cntrd, local_background_2Dg, radial_spot_size2D, gaussian_spot_size2D, feretDiameters

clear
clc
tic

%%============================== PROJECT DIRECTORY SETUP=============================
here = fileparts(mfilename('fullpath'));   % folder this script lives in
addpath(fullfile(here, 'functions'));
addpath(fullfile(here, 'bioformats'));
javaaddpath(fullfile(here, 'bioformats', 'bioformats_package.jar'));
%%==============================%%==============================%%===================


k=1; %%%%%single images

save_filename='/Volumes/xiel2lab/Aniket/BRD4 MS2 enh 10-14-25/Single images/brd4 ctrl 14_SIM condensate intensity analysis V4.mat';
MIP_filename  = '/Volumes/xiel2lab/Aniket/Brd4 drug 2-25-2026/Ctrl/Ctrl_MIPs/ctrl 11 100 ms_SIM_Maximum intensity projection.czi';

channel_cond = 3;      % 642 channel

scene=1;    % set manually

% % % READ IMAGE%%%=================================================

MIP_out = ReadImage6D2(MIP_filename, true, scene);
MIP_image6d = MIP_out{1};

metadata = MIP_out{2};

xpixel = 0.0313;
ypixel = 0.0313;


%==============MIP =============

inputImage_cond=squeeze(MIP_image6d(1,k,1,channel_cond,:,:));

% =================NUCLEUS SEGMENTATION =====================

seg_OUT = cellpose_seg(inputImage_cond,150);  % roughly segment the nuclei

%%================% %%================

% % % --- Access individual masks ---
individual_masks = seg_OUT{1};   % cell array containing all masks, one binary mask per nucleus
% --- Remove border-touching nuclei ---
[rows, cols] = size(individual_masks{1});  % image dimensions

% Set your threshold and convert it to pixels
% Assuming 'xpixel' is your pixel size in um/pixel (e.g., 0.1)
distance_threshold_um = 5;
max_border_pixels = distance_threshold_um / xpixel;

% Preallocate 'keep' as an array of true values
keep = true(1, length(individual_masks));

for k = 1:length(individual_masks)
    m = individual_masks{k};

    % Sum the true pixels along each edge
    % We index the left and right edges from 2 to end-1 to prevent
    % double-counting the corner pixels if a nucleus sits exactly in a corner.
    top_touch    = sum(m(1, :));
    bottom_touch = sum(m(end, :));
    left_touch   = sum(m(2:end-1, 1));
    right_touch  = sum(m(2:end-1, end));

    % Total number of border-touching pixels
    total_touch_pixels = top_touch + bottom_touch + left_touch + right_touch;

    % Reject ONLY if the total touch length is greater than the threshold
    if total_touch_pixels > max_border_pixels
        keep(k) = false;
    end
end

% Filter the masks
individual_masks = individual_masks(keep);
fprintf('Kept %d / %d nuclei (removed %d border nuclei)\n', ...
    sum(keep), length(keep), sum(~keep));


figure;

% % % % % Colored label overlay on original image
label_clean = zeros(rows, cols);
for k = 1:length(individual_masks)
    label_clean(individual_masks{k}) = k;
end

% Show original with colored overlay
imshow(inputImage_cond, []);
hold on;
overlay = labeloverlay(imadjust(im2uint8(mat2gray(inputImage_cond))), label_clean, ...
    'Colormap', lines(length(individual_masks)), ...
    'Transparency', 0.5);
imshow(overlay, 'InitialMagnification','Fit');
title(sprintf('%d nuclei after border removal', length(individual_masks)));

num_nuclei = length(individual_masks);

% Epsilon: roughly the expected radius of a punctum in pixels
% MinPts: roughly the number of pixels in the smallest valid punctum
eps_val = 8;
min_pixels = 40;

% ---------------- Preallocate ----------------
all_condensate_counts        = zeros(num_nuclei, 1);
all_condensate_Areas         = cell(num_nuclei, 1);
all_condensate_centroids     = cell(num_nuclei, 1);
all_condensate_peak_raw      = cell(num_nuclei, 1);   % raw peak intensity per condensate
all_condensate_mean_raw      = cell(num_nuclei, 1);   % raw mean intensity per condensate
all_condensate_peak_int_norm = cell(num_nuclei, 1);   % peak / nucleus mean  <-- your metric
all_condensate_mean_int_norm = cell(num_nuclei, 1);   % mean / nucleus mean  (partition coeff.)

nucleus_mean_int         = nan(num_nuclei, 1);   % raw mean over whole nucleus
nucleus_nucleoplasm_int  = nan(num_nuclei, 1);   % raw mean excluding condensate pixels
nucleus_median_peak_norm = nan(num_nuclei, 1);   % per-nucleus summary of the ratio

raw = double(inputImage_cond);   % <-- every intensity measurement comes from here

for j = 1:num_nuclei
    mask_j = logical(individual_masks{j});

    % ---------- DETECTION: processed image only ----------
    cond_I   = raw .* mask_j;
    % cond_I_M = SNR_adjust(cond_I, 0.1, 0.5);
    % intensity_thresh = mean(nonzeros(cond_I_M(:))) + 2*std(nonzeros(cond_I_M(:)));
    intensity_thresh = mean(nonzeros(cond_I(:))) + 2*std(nonzeros(cond_I(:)));

    [labels, stats, mask] = detectPunctaDBSCAN(cond_I, 'Epsilon', eps_val, ...
        'TophatRadius', 30, 'Threshold', intensity_thresh, ...
        'MinArea', min_pixels);

 
    labels_img = bwlabel(logical(mask) & mask_j);

    n_bw = max(labels_img(:));
    if n_bw ~= numel(stats)
        warning('Nucleus %d: bwlabel found %d regions, detector reported %d.', ...
                j, n_bw, numel(stats));
    end

    % ---------- PHOTOMETRY: raw image only ----------
    nucleus_mean_int(j) = mean(raw(mask_j));

    cond_pix        = mask_j & (labels_img > 0);
    nucleoplasm_pix = mask_j & ~cond_pix;
    if any(nucleoplasm_pix(:))
        nucleus_nucleoplasm_int(j) = mean(raw(nucleoplasm_pix));
    end

    rp = regionprops(labels_img, raw, 'Area', 'Centroid', 'MaxIntensity', 'MeanIntensity');
    if ~isempty(rp)
        rp = rp([rp.Area] > 0);
    end

    all_condensate_counts(j) = numel(rp);

    if ~isempty(rp)
        all_condensate_Areas{j}     = [rp.Area]' * (xpixel * ypixel);   % um^2
        all_condensate_centroids{j} = cat(1, rp.Centroid);              % [x y] px

        peak_raw = [rp.MaxIntensity]';
        mean_raw = [rp.MeanIntensity]';

        all_condensate_peak_raw{j}      = peak_raw;
        all_condensate_mean_raw{j}      = mean_raw;
        all_condensate_peak_int_norm{j} = peak_raw ./ nucleus_mean_int(j);
        all_condensate_mean_int_norm{j} = mean_raw ./ nucleus_mean_int(j);

        nucleus_median_peak_norm(j) = median(all_condensate_peak_int_norm{j});
    else
        all_condensate_Areas{j}         = [];
        all_condensate_centroids{j}     = [];
        all_condensate_peak_raw{j}      = [];
        all_condensate_mean_raw{j}      = [];
        all_condensate_peak_int_norm{j} = [];
        all_condensate_mean_int_norm{j} = [];
    end
end

% ================= POOL ACROSS ALL NUCLEI =================
pooled_peak_int_norm = vertcat(all_condensate_peak_int_norm{:});
pooled_mean_int_norm = vertcat(all_condensate_mean_int_norm{:});
pooled_peak_raw      = vertcat(all_condensate_peak_raw{:});
pooled_mean_raw      = vertcat(all_condensate_mean_raw{:});
pooled_areas         = vertcat(all_condensate_Areas{:});

% nucleus tag for each pooled condensate (keeps per-nucleus structure for stats)
pooled_nucleus_id = repelem((1:num_nuclei)', all_condensate_counts);

fprintf('\n%d condensates pooled across %d nuclei\n', ...
    numel(pooled_peak_int_norm), num_nuclei);
fprintf('Peak / nucleus-mean: mean %.2f, median %.2f, IQR [%.2f  %.2f]\n', ...
    mean(pooled_peak_int_norm), median(pooled_peak_int_norm), ...
    prctile(pooled_peak_int_norm,25), prctile(pooled_peak_int_norm,75));

figure;
histogram(pooled_peak_int_norm, 60, 'Normalization','probability');
xlabel('Peak hub intensity / nucleus mean intensity');
ylabel('Fraction of condensates');
title(sprintf('Pooled condensate enrichment (n = %d)', numel(pooled_peak_int_norm)));

save(save_filename, ...
    'all_condensate_counts','all_condensate_Areas','all_condensate_centroids', ...
    'all_condensate_peak_raw','all_condensate_mean_raw', ...
    'all_condensate_peak_int_norm','all_condensate_mean_int_norm', ...
    'nucleus_mean_int','nucleus_nucleoplasm_int','nucleus_median_peak_norm', ...
    'pooled_peak_int_norm','pooled_mean_int_norm','pooled_peak_raw','pooled_mean_raw', ...
    'pooled_areas','pooled_nucleus_id','individual_masks')

toc
