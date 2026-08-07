% % % % Hub_size_calculation_PIPELINE.m

% % %  Single/Multi-channel SR-SIM images: Batch processing of multiple images
% % %
% % %  Pipeline:
% % %    Step 1 : Setup image directory for batch processing
% % %    Step 2 : Load individual image CZI (all channels) 
% % %    Step 3 : Nucleus segmentation using Cellpose
% % %    Step 4 : Global condensate/Hub peak detection in each nucleus
% % %    Step 5 : Hub centroid refinement using Gaussian fitting
% % %    Step 6 : Estimate Local intensity background radially around centroid
% % %    Step 7 : Calculate Feret's diameter of eachg hub
% % %    Step 8:  Pool all hub measurements and save .mat
% % %
% % %  Dependencies on path:
% % %  ReadImage6D2, cellpose_seg, pkfnd, cntrd, local_background_2Dg, radial_spot_size2D, gaussian_spot_size2D, feretDiameters

clear; clc; tic;

%%============================== PROJECT DIRECTORY SETUP=============================
here = fileparts(mfilename('fullpath'));   % folder this script lives in
addpath(fullfile(here, 'functions'));
addpath(fullfile(here, 'bioformats'));
javaaddpath(fullfile(here, 'bioformats', 'bioformats_package.jar'));
%%==============================%%==============================%%===================

folderPath = '/Volumes/xiel2lab/Aniket/4c BRD4 MED14 dual/2 color /';
fileList = dir(fullfile(folderPath, '*.czi'));

for f = 1:length(fileList)
    %% --- 1) Build file paths and save name ---
    MIP_filename = fullfile(folderPath, fileList(f).name);
    [filepath, baseName, ~] = fileparts(MIP_filename);

    fprintf('Processing file %d/%d: %s\n', f, length(fileList), fileList(f).name);


scene = 1;

%% 2) Read individual CZI files
MIP_out = ReadImage6D2(MIP_filename, true, scene);
MIP_image6d = MIP_out{1};

xpixel = 0.0313;       % um/pixel
channel_cond = 1;      % condensate channel

inputImage_cond = squeeze(MIP_image6d(1,1,1,channel_cond,:,:));
img = double(inputImage_cond);

%% ---- 3) Nucleus segmentation ----
seg_OUT = cellpose_seg(inputImage_cond, 150);
individual_masks = seg_OUT{1};
numNuclei = numel(individual_masks);
fprintf('Detected %d nuclei.\n', numNuclei);

%% ---- 4) Find condensates/hubs within each nucleus using pkfnd ----
pkfnd_sz = 25;           % minimum separation between peaks (pixels)
edge_margin = 10;       % exclude spots near nucleus edge (pixels)
thresh_mult = 3;        % threshold = thresh_mult × mean nuclear intensity

C_in = [];
nucleus_id = [];

for n = 1:numNuclei
    mask_n = individual_masks{n};

    % Erode mask to avoid edge artifacts
    mask_eroded = imerode(mask_n, strel('disk', edge_margin));

    % Compute adaptive threshold from mean intensity inside this nucleus
    nuc_pixels = img(mask_eroded > 0);
    if isempty(nuc_pixels)
        continue;
    end
    nuc_mean = mean(nuc_pixels);
    pkfnd_thresh = thresh_mult * nuc_mean;

    % Apply mask and find peaks
    img_masked = img .* mask_eroded;
    pk = pkfnd(img_masked, pkfnd_thresh, pkfnd_sz);

    if isempty(pk)
        continue;
    end

    % Verify each peak is inside the eroded mask
    valid = false(size(pk, 1), 1);
    for j = 1:size(pk, 1)
        px = round(pk(j, 1));
        py = round(pk(j, 2));
        if py >= 1 && py <= size(mask_eroded, 1) && ...
           px >= 1 && px <= size(mask_eroded, 2)
            valid(j) = mask_eroded(py, px) > 0;
        end
    end
    pk = pk(valid, :);

    if ~isempty(pk)
        C_in = [C_in; pk];
        nucleus_id = [nucleus_id; n * ones(size(pk, 1), 1)];
        fprintf('  Nucleus %d: mean = %.0f, thresh = %.0f, found %d spots\n', ...
            n, nuc_mean, pkfnd_thresh, size(pk, 1));
    end
end

numPoints = size(C_in, 1);
fprintf('Found %d condensate candidates across %d nuclei.\n', numPoints, numNuclei);

if numPoints == 0
    error('No condensates detected. Try lowering thresh_mult (currently %.1f).', thresh_mult);
end

%% ---- 5) Refine centroids (Crocker-Grier cntrd → 2D Gaussian) and measure ----
sz_cntrd = 21;
R_search = 20;

mx = round(C_in);
centroids = cntrd(img, mx, sz_cntrd);

% Handle cntrd edge rejection
if size(centroids, 1) < numPoints
    warning('cntrd rejected %d points near image edges.', numPoints - size(centroids, 1));
    numPoints = size(centroids, 1);
    kept_idx = zeros(numPoints, 1);
    for i = 1:numPoints
        dists = sqrt((C_in(:,1) - centroids(i,1)).^2 + (C_in(:,2) - centroids(i,2)).^2);
        [~, kept_idx(i)] = min(dists);
    end
    C_in = C_in(kept_idx, :);
    nucleus_id = nucleus_id(kept_idx);
end

% Allocate outputs
condensate_cent    = zeros(numPoints, 2);
condensate_FWHM    = zeros(numPoints, 1);
condensate_AR      = zeros(numPoints, 1);
condensate_AR_feret= zeros(numPoints, 1);
condensate_FWHM_r  = zeros(numPoints, 1);
condensate_feret   = zeros(numPoints, 1);
condensate_bg      = zeros(numPoints, 1);
condensate_sigma_r = zeros(numPoints, 1);

[m, n_cols] = size(img);

for i = 1:numPoints
    x0 = centroids(i, 1);
    y0 = centroids(i, 2);
    center = [y0, x0];

    %% ---- 6) Radial profile and local background ----
    [loc_b, sigma_bg] = local_background_2Dg(img, center, R_search);
    condensate_bg(i) = loc_b;
    condensate_sigma_r(i) = sigma_bg;

    %% ---- 7) Background-subtracted ROI and Feret's diameter ----
    x_min = max(1, round(x0) - R_search);
    x_max = min(n_cols, round(x0) + R_search);
    y_min = max(1, round(y0) - R_search);
    y_max = min(m, round(y0) + R_search);

    ROI = img(y_min:y_max, x_min:x_max) - loc_b;
    ROI(ROI < 0) = 0;

    % Binarize
    ROI_norm = ROI / max(ROI(:));
    level = graythresh(ROI_norm);
    bw = imbinarize(ROI_norm, level);

    % Keep only central object
    bw_label = bwlabel(bw);
    center_label = bw_label(round(size(bw,1)/2), round(size(bw,2)/2));
    if center_label > 0
        bw = (bw_label == center_label);
    end

    % Feret's diameter
    try
        stats = regionprops(bw, 'MaxFeretDiameter', 'MinFeretDiameter');
        if ~isempty(stats)
            feret_avg = (stats(1).MaxFeretDiameter + stats(1).MinFeretDiameter) / 2;
            condensate_feret(i) = feret_avg * xpixel * 1000;
            condensate_AR_feret(i) = (stats(1).MaxFeretDiameter)/(stats(1).MinFeretDiameter);
        else
            condensate_feret(i) = NaN;
            condensate_AR_feret(i) = NaN;
        end
    catch
        [maxF, minF] = feretDiameters(bw);
        condensate_feret(i) = (maxF + minF) / 2 * xpixel * 1000;
        condensate_AR_feret(i) = maxF/minF;
    end

    %% ---- FWHM from radial profile ----
    ROI_center = [round(size(ROI,1)/2), round(size(ROI,2)/2)];
    condensate_FWHM_r(i) = radial_spot_size2D(ROI, ROI_center, R_search, xpixel);

    %% ---- 2D Gaussian fit ----
    outp = gaussian_spot_size2D(ROI, ...
        x0 - x_min + 1, y0 - y_min + 1, xpixel, R_search);
    condensate_FWHM(i) = outp(1);
    condensate_AR(i)   = outp(2);
    condensate_cent(i,1) = outp(3) + x_min - 1;
    condensate_cent(i,2) = outp(4) + y_min - 1;
end

%% ---- Visualization ----
figure;
imshow(uint16(65536 - inputImage_cond), []);
hold on;

colors = lines(numNuclei);
for n = 1:numNuclei
    idx = (nucleus_id == n);
    if any(idx)
        plot(condensate_cent(idx,1), condensate_cent(idx,2), '*', ...
            'Color', colors(n,:), 'MarkerSize', 6);
    end
end

for n = 1:numNuclei
    B = bwboundaries(individual_masks{n});
    for b = 1:numel(B)
        plot(B{b}(:,2), B{b}(:,1), '-', 'Color', colors(n,:), 'LineWidth', 3);
    end
end
hold off;
title(sprintf('%d condensates in %d nuclei (thresh = %d x mean)', ...
    numPoints, numNuclei, thresh_mult));

%% ---- 8) Print summary table ----
fprintf('\n%-6s %6s %10s %10s %10s %10s %10s\n', ...
    'Spot', 'Nuc', 'FWHM_2D', 'FWHM_rad', 'Feret', 'AR', 'Bkg');
fprintf('%-6s %6s %10s %10s %10s %10s %10s\n', ...
    '', '', '(nm)', '(nm)', '(nm)', '', '(counts)');
fprintf('%s\n', repmat('-', 1, 68));
for i = 1:numPoints
    fprintf('%-6d %6d %10.1f %10.1f %10.1f %10.2f %10.1f\n', ...
        i, nucleus_id(i), condensate_FWHM(i), condensate_FWHM_r(i), ...
        condensate_feret(i), condensate_AR(i), condensate_bg(i));
end

% ----  Per-nucleus summary ----
fprintf('\n--- Per-nucleus summary ---\n');
fprintf('%-6s %6s %12s %12s %12s\n', 'Nuc', 'Count', 'FWHM_2D', 'Feret', 'AR');
fprintf('%-6s %6s %12s %12s %12s\n', '', '', 'mean+-sd', 'mean+-sd', 'mean+-sd');
fprintf('%s\n', repmat('-', 1, 55));
for n = 1:numNuclei
    idx = (nucleus_id == n);
    if sum(idx) == 0, continue; end
    fprintf('%-6d %6d %5.0f +- %-5.0f %5.0f +- %-5.0f %4.2f +- %-4.2f\n', ...
        n, sum(idx), ...
        mean(condensate_FWHM(idx)), std(condensate_FWHM(idx)), ...
        mean(condensate_feret(idx)), std(condensate_feret(idx)), ...
        mean(condensate_AR(idx)), std(condensate_AR(idx)));
end


% ---- Save ----
save_filename = fullfile(filepath, [baseName, '.mat']);

save(save_filename, 'C_in', 'condensate_cent', 'nucleus_id', 'numNuclei',...
    'condensate_FWHM', 'condensate_FWHM_r', 'condensate_feret', ...
    'condensate_AR', 'condensate_AR_feret','condensate_bg', 'condensate_sigma_r');

fprintf('\nSaved to %s\n', save_filename);

end

toc