%% Alternative Puncta Detection Methods
% This script provides different approaches for separating puncta from background
% Use this to test which method works best for your data

function [puncta_mask, background_intensity] = detect_puncta_alternative(nuclear_image, nucleus_mask, method)
% Detects puncta within a nucleus using various methods
%
% Inputs:
%   nuclear_image - intensity image (full image)
%   nucleus_mask - binary mask of single nucleus
%   method - string specifying detection method
%
% Outputs:
%   puncta_mask - binary mask of detected puncta
%   background_intensity - mean background intensity

    % Extract nuclear region
    nuclear_values = nuclear_image(nucleus_mask);
    nuclear_region = nuclear_image .* double(nucleus_mask);
    
    switch method
        case 'percentile'
            % Method 1: Top percentile threshold
            threshold_percentile = 75; % Adjust: 70-90 typical
            threshold = prctile(nuclear_values, threshold_percentile);
            puncta_mask = (nuclear_region > threshold) & nucleus_mask;
            
        case 'std_based'
            % Method 2: Mean + N*STD
            mean_int = mean(nuclear_values);
            std_int = std(nuclear_values);
            n_std = 2.5; % Adjust: 2-3 typical
            threshold = mean_int + n_std * std_int;
            puncta_mask = (nuclear_region > threshold) & nucleus_mask;
            
        case 'triangle'
            % Method 3: Triangle thresholding (good for bimodal distributions)
            hist_counts = histcounts(nuclear_values, 256);
            threshold_normalized = triangle_threshold(hist_counts);
            min_val = min(nuclear_values);
            max_val = max(nuclear_values);
            threshold = min_val + threshold_normalized * (max_val - min_val);
            puncta_mask = (nuclear_region > threshold) & nucleus_mask;
            
        case 'local_maxima'
            % Method 4: Local maxima detection (good for distinct puncta)
            % Smooth image slightly
            smoothed = imgaussfilt(nuclear_region, 1);
            
            % Find regional maxima
            regional_max = imregionalmax(smoothed .* double(nucleus_mask));
            
            % Dilate maxima to capture puncta area
            se = strel('disk', 2);
            puncta_mask = imdilate(regional_max, se) & nucleus_mask;
            
            % Threshold to keep only bright maxima
            threshold = prctile(nuclear_values, 70);
            puncta_mask = puncta_mask & (nuclear_region > threshold);
            
        case 'top_hat'
            % Method 5: Top-hat filtering (enhances bright spots)
            se = strel('disk', 5); % Adjust size based on expected puncta size
            tophat = imtophat(nuclear_region, se);
            
            % Threshold the top-hat image
            threshold = graythresh(tophat(nucleus_mask)) * max(tophat(:));
            puncta_mask = (tophat > threshold) & nucleus_mask;
            
        case 'laplacian'
            % Method 6: Laplacian of Gaussian (blob detection)
            sigma = 1.5; % Adjust based on puncta size
            h = fspecial('log', round(6*sigma), sigma);
            log_filtered = imfilter(double(nuclear_region), h, 'replicate');
            
            % Negative values indicate bright blobs
            log_filtered = -log_filtered;
            threshold = prctile(log_filtered(nucleus_mask), 95);
            puncta_mask = (log_filtered > threshold) & nucleus_mask;
            
        case 'adaptive'
            % Method 7: Adaptive thresholding
            % Only works on nuclear region, so extract it first
            stats = regionprops(nucleus_mask, 'BoundingBox');
            if ~isempty(stats)
                bbox = round(stats(1).BoundingBox);
                % Ensure bbox is within image bounds
                x1 = max(1, bbox(1));
                y1 = max(1, bbox(2));
                x2 = min(size(nuclear_region,2), bbox(1)+bbox(3));
                y2 = min(size(nuclear_region,1), bbox(2)+bbox(4));
                
                nuclear_crop = nuclear_region(y1:y2, x1:x2);
                mask_crop = nucleus_mask(y1:y2, x1:x2);
                
                % Adaptive threshold
                T = adaptthresh(mat2gray(nuclear_crop), 0.4); % Adjust sensitivity
                puncta_crop = imbinarize(mat2gray(nuclear_crop), T) & mask_crop;
                
                % Place back in full image
                puncta_mask = false(size(nuclear_region));
                puncta_mask(y1:y2, x1:x2) = puncta_crop;
            else
                puncta_mask = false(size(nuclear_region));
            end
            
        case 'kmeans'
            % Method 8: K-means clustering (2 clusters: background vs puncta)
            if sum(nucleus_mask(:)) > 50
                [~, centers] = kmeans(nuclear_values, 2, 'MaxIter', 100);
                % Higher center is puncta
                puncta_cluster = centers == max(centers);
                cluster_idx = kmeans(nuclear_values, 2, 'MaxIter', 100);
                
                % Create mask
                puncta_mask = false(size(nuclear_region));
                puncta_mask(nucleus_mask) = (cluster_idx == find(puncta_cluster));
            else
                puncta_mask = false(size(nuclear_region));
            end
            
        otherwise
            error('Unknown method: %s', method);
    end
    
    % Clean up puncta mask
    puncta_mask = bwareaopen(puncta_mask, 3); % Remove tiny objects
    
    % Calculate background (nucleus excluding puncta)
    background_mask = nucleus_mask & ~puncta_mask;
    if any(background_mask(:))
        background_intensity = mean(nuclear_image(background_mask));
    else
        background_intensity = mean(nuclear_values); % Fallback
    end
end

function threshold = triangle_threshold(histogram)
    % Triangle thresholding algorithm
    % Find peak
    [~, peak_idx] = max(histogram);
    
    % Find last non-zero bin
    last_idx = find(histogram > 0, 1, 'last');
    
    % Calculate distances from line connecting peak to end
    x_values = peak_idx:last_idx;
    line_y = linspace(histogram(peak_idx), histogram(last_idx), length(x_values));
    distances = abs(histogram(x_values) - line_y);
    
    % Find maximum distance
    [~, max_dist_idx] = max(distances);
    threshold = (x_values(max_dist_idx) - 1) / 255; % Normalize to [0,1]
end

%% Example usage comparing methods
% Run this section after loading your data

% Assuming you have: inputImage_cond, nuclear_labels, num_nuclei

methods_to_test = {'percentile', 'std_based', 'local_maxima', 'top_hat'};
nuc_idx = 1; % Test on first nucleus

figure('Name', 'Method Comparison', 'Position', [100 100 1200 800]);
for m = 1:length(methods_to_test)
    current_method = methods_to_test{m};
    current_nucleus_mask = (nuclear_labels == nuc_idx);
    
    [puncta_mask, bg_int] = detect_puncta_alternative(inputImage_cond, ...
        current_nucleus_mask, current_method);
    
    subplot(2, 4, m);
    nuclear_brd4 = inputImage_cond .* double(current_nucleus_mask);
    imshow(labeloverlay(mat2gray(nuclear_brd4), puncta_mask, 'Colormap', 'jet'));
    title(sprintf('%s\nBG: %.1f, Puncta: %d', current_method, ...
        bg_int, max(bwlabel(puncta_mask))));
    
    subplot(2, 4, m+4);
    nuclear_values = inputImage_cond(current_nucleus_mask);
    histogram(nuclear_values, 50, 'FaceAlpha', 0.5);
    hold on;
    xline(bg_int, 'r--', 'LineWidth', 2);
    title('Intensity Distribution');
    xlabel('Intensity');
    ylabel('Pixels');
end

fprintf('Method comparison complete. Choose the method that best separates puncta.\n');
