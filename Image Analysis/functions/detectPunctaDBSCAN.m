function [labels, stats, mask] = detectPunctaDBSCAN(img, varargin)
% detectPunctaDBSCAN  Detect fluorescent puncta/condensates via DBSCAN.
%
% USAGE:
%   [labels, stats, mask] = detectPunctaDBSCAN(img)
%   [labels, stats, mask] = detectPunctaDBSCAN(img, 'Epsilon', 3, ...)
%
% INPUTS:
%   img  - 2D MIP image (any numeric type; converted to double internally)
%
% NAME-VALUE PARAMETERS:
%   Epsilon       - DBSCAN neighbourhood radius in pixels      (default: 3)
%   MinPts        - Minimum points to seed a cluster            (default: 5)
%   Threshold     - Intensity threshold ('auto' | scalar)       (default: 'auto')
%   SigmaBlur     - Gaussian pre-filter sigma in pixels         (default: 1.0)
%   TophatRadius  - Disk radius for rolling-ball background     (default: 15)
%                   subtraction. Set 0 to skip.
%   MinArea       - Discard puncta smaller than this (px)       (default: 3)
%   MaxArea       - Discard puncta larger  than this (px)       (default: Inf)
%   UseIntensity  - Weight DBSCAN coords by intensity (logical) (default: true)
%   IntensityWeight - Scaling factor for the intensity axis     (default: 0.5)
%   Visualize     - Show diagnostic figure (logical)            (default: false)

    %% ---- Parse inputs -------------------------------------------------------
    p = inputParser;
    addRequired(p, 'img', @isnumeric);
    addParameter(p, 'Epsilon',         3,      @isnumeric);
    addParameter(p, 'MinPts',          5,      @isnumeric);
    addParameter(p, 'Threshold',       'auto');
    addParameter(p, 'SigmaBlur',       1.0,    @isnumeric);
    addParameter(p, 'TophatRadius',    30,     @isnumeric);
    addParameter(p, 'MinArea',         10,      @isnumeric);
    addParameter(p, 'MaxArea',         Inf,    @isnumeric);
    addParameter(p, 'UseIntensity',    true,   @islogical);
    addParameter(p, 'IntensityWeight', 0.5,    @isnumeric);
    addParameter(p, 'Visualize',       false,  @islogical);
    parse(p, img, varargin{:});
    opts = p.Results;

    %% ---- 1. Preprocessing ---------------------------------------------------
    img = double(img);

    % img = subtractBackground(img,250);        % sigma = 100 (default)

    % 1a. Background subtraction via morphological top-hat
    if opts.TophatRadius > 0
        se  = strel('disk', opts.TophatRadius);
        img = imtophat(img, se);
    end

    % 1b. Gentle Gaussian blur to suppress shot noise while preserving puncta
    if opts.SigmaBlur > 0
        img = imgaussfilt(img, opts.SigmaBlur);
    end

    % Normalise to [0 1] for consistent downstream handling
    imgMin = min(img(:));
    imgMax = max(img(:));
    if imgMax - imgMin < eps
        warning('Image has no dynamic range after preprocessing.');
        labels = []; stats = []; mask = zeros(size(img));
        return;
    end
    imgNorm = (img - imgMin) / (imgMax - imgMin);

    %% ---- 2. Adaptive thresholding -------------------------------------------
    if ischar(opts.Threshold) || isstring(opts.Threshold)
        % Robust auto-threshold: mean + 2*std of non-zero pixels,
        % cross-checked against Otsu on the preprocessed image.
        vals       = imgNorm(imgNorm > 0);
        statThresh = mean(vals) + 2 * std(vals);
        otsuThresh = graythresh(imgNorm);          % Otsu's method
        threshold  = max(statThresh, otsuThresh);  % take the more conservative
    else
        % User-supplied threshold: rescale to match normalised image
        threshold = (double(opts.Threshold) - imgMin) / (imgMax - imgMin);
    end

    binaryMask = imgNorm > threshold;

    % Clean up: remove isolated single-pixel noise, fill tiny holes
    binaryMask = bwareaopen(binaryMask, max(1, floor(opts.MinPts / 2)));
    binaryMask = imclose(binaryMask, strel('disk', 1));

    [rows, cols] = find(binaryMask);
    if isempty(rows)
        warning('No pixels survive thresholding (T = %.4f). Consider lowering the threshold.', threshold);
        labels = []; stats = []; mask = zeros(size(img));
        return;
    end

    %% ---- 3. Build feature matrix for DBSCAN ---------------------------------
    % Vanilla DBSCAN on (x,y) alone treats every above-threshold pixel equally.
    % Adding a scaled intensity axis lets the algorithm separate nearby but
    % distinct puncta whose intensity profiles don't overlap.
    intensities = imgNorm(sub2ind(size(imgNorm), rows, cols));

    if opts.UseIntensity
        % Scale intensity axis so that it meaningfully contributes to distance
        % without overwhelming spatial proximity.
        X = [double(cols), double(rows), intensities * opts.IntensityWeight * opts.Epsilon];
    else
        X = [double(cols), double(rows)];
    end

    %% ---- 4. DBSCAN clustering -----------------------------------------------
    labels = dbscan(X, opts.Epsilon, opts.MinPts);

    %% ---- 5. Build labelled mask & relabel contiguously ----------------------
    mask = zeros(size(img));
    validIdx       = labels > 0;                            % drop noise (label -1)
    validRows      = rows(validIdx);
    validCols      = cols(validIdx);
    validLabels    = labels(validIdx);

    % Relabel to 1:N (DBSCAN may skip integers)
    [uniqueIDs, ~, remapped] = unique(validLabels);
    nClusters = numel(uniqueIDs);

    linInd = sub2ind(size(img), validRows, validCols);
    mask(linInd) = remapped;

    %% ---- 6. Size-based filtering --------------------------------------------
    convergeMask = false;
    finalMap      = 1:nClusters;   % maps old cluster id -> new (0 = removed)

    for k = 1:nClusters
        area = nnz(mask == k);
        if area < opts.MinArea || area > opts.MaxArea
            mask(mask == k)  = 0;
            finalMap(k)      = 0;
            convergeMask     = true;
        end
    end

    % Relabel again after removal so IDs stay contiguous
    if convergeMask
        props   = regionprops(logical(mask), 'PixelIdxList');
        newMask = zeros(size(img));
        for k = 1:numel(props)
            newMask(props(k).PixelIdxList) = k;
        end
        mask = newMask;
    end

    %% ---- 7. Extract puncta statistics (on original preprocessed intensities) -
    stats = regionprops(mask, img, ...
        'Centroid', 'WeightedCentroid', 'Area', ...
        'MeanIntensity', 'MaxIntensity', 'PixelValues', ...
        'BoundingBox', 'Eccentricity');

    % Append integrated intensity (sum of pixel values) — useful for
    % quantifying total condensate signal.
    for k = 1:numel(stats)
        stats(k).IntegratedIntensity = sum(stats(k).PixelValues);
    end

    %% ---- 8. Update label vector to match final mask -------------------------
    finalLabels = zeros(size(labels));
    for i = 1:numel(rows)
        finalLabels(i) = mask(rows(i), cols(i));
    end
    labels = finalLabels;

    %% ---- 9. Optional diagnostic visualisation --------------------------------
    if opts.Visualize
        convergeFig(img, imgNorm, binaryMask, mask, stats, threshold);
    end
end

%% =========================================================================
%  HELPER: diagnostic figure
%  =========================================================================
function convergeFig(imgRaw, imgNorm, binaryMask, mask, stats, thresh)
    figure('Name', 'Puncta Detection Diagnostics', 'NumberTitle', 'off', ...
           'Position', [100 100 1400 400]);

    % Panel 1 — preprocessed image with threshold contour
    ax1 = subplot(1,3,1);
    imshow(imgNorm, [], 'Parent', ax1); hold(ax1, 'on');
    contour(ax1, imgNorm, [thresh thresh], 'r', 'LineWidth', 0.8);
    title(ax1, sprintf('Preprocessed + Threshold (T=%.3f)', thresh));

    % Panel 2 — binary mask fed into DBSCAN
    ax2 = subplot(1,3,2);
    imshow(binaryMask, [], 'Parent', ax2);
    title(ax2, sprintf('Thresholded Mask (%d px)', nnz(binaryMask)));

    % Panel 3 — final labelled puncta with centroids
    ax3 = subplot(1,3,3);
    nPuncta = max(mask(:));
    if nPuncta > 0
        rgb = label2rgb(mask, 'jet', 'k', 'shuffle');
        imshow(rgb, 'Parent', ax3); hold(ax3, 'on');
        centroids = cat(1, stats.WeightedCentroid);
        plot(ax3, centroids(:,1), centroids(:,2), 'w+', 'MarkerSize', 8, 'LineWidth', 1.5);
    else
        imshow(zeros(size(mask)), 'Parent', ax3);
    end
    title(ax3, sprintf('Detected Puncta: %d', nPuncta));
end