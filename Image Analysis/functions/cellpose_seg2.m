function out = cellpose_seg2(inputImage_cond, cell_D)
    
% --- IMPROVED PREPROCESSING FOR FAINT SIGNAL ---
    img_double = double(inputImage_cond);
    
    % 1. Percentile Clipping (Cap the bright hubs)
    p_high = prctile(img_double(:), 98); 
    img_clipped = min(img_double, p_high);
    
    % 2. Initial Normalization to 0-1 range
    img_norm = (img_clipped - min(img_clipped(:))) ./ (max(img_clipped(:)) - min(img_clipped(:)));
    
    % 3. Aggressive Histogram Stretching 
    % This ignores the bottom 5% of pixels (background noise) and stretches 
    % the rest to fill the 0-1 range, maximizing the contrast of the faint signal.
    limits = stretchlim(img_norm, [0.05 0.99]); 
    img_stretched = imadjust(img_norm, limits, []);
    
    % 4. Median Filtering (Erase puncta)
    img_smoothed = medfilt2(img_stretched, [5 5]); 
    
    % 5. Stronger Gamma Correction
    % Dropping from 0.5 to 0.4 to pull up those very dark mid-tones even more
    img_gamma = imadjust(img_smoothed, [], [], 0.4);
    
    % Convert to uint16
    img_uint16 = uint16(img_gamma * 65535);

    % --- DOWNSAMPLING ---
    scale_factor = 0.3; % Reduce size by 70%
    img_resized = imresize(img_uint16, scale_factor, 'bicubic');

    % --- Cellpose Segmentation (Python Environment) ---
    cp_models = py.importlib.import_module('cellpose.models');
    np        = py.importlib.import_module('numpy');

    % Convert downsampled image to numpy array
    img_np = np.array(img_resized);

    % Load nuclei model (CPU only)
model = cp_models.Cellpose(pyargs('gpu', false, 'model_type', 'nuclei'));

% % % % Run segmentation with cell_D
results = model.eval(img_np, pyargs(...
    'diameter',  cell_D, ...
    'channels',  py.list({0, 0}), ...
    'augment',   false));

    % Extract resized mask
    masks_resized = uint16(double(py.numpy.array(results{1})));

% --- UPSAMPLING & FILTERING ---
    original_size = size(img_uint16);
    
    % CRITICAL: Use 'nearest' interpolation so mask label IDs are preserved!
    masks = imresize(masks_resized, original_size, 'nearest');
    
    % 1. FILTER BY SIZE (Without destroying unique IDs)
    % Measure the area of every unique cell ID
    stats = regionprops(masks, 'Area');
    
    % Find which Cell IDs meet your pixel threshold
    valid_ids = find([stats.Area] >= 10000); 
    
    % Keep only the pixels that belong to those valid IDs
    masks_filtered = masks;
    masks_filtered(~ismember(masks, valid_ids)) = 0;
    
    % 2. EXPAND THE MASKS (Using Distance Transform)
    dil_radius = 10; % Set your desired expansion in pixels here
    
    % Calculate distance from background to the nearest cell pixel
    % D holds the distance, nearest_idx holds the linear index of the closest cell pixel
    [D, nearest_idx] = bwdist(masks_filtered > 0);
    
    % Define the expansion zone: pixels within the radius that are currently background
    expansion_zone = (D <= dil_radius) & (masks_filtered == 0);
    
    % Safely push the labels outward into the expansion zone
    masks_filtered(expansion_zone) = masks_filtered(nearest_idx(expansion_zone));

    % 3. RE-NUMBER LABELS (To ensure they are strictly sequential 1 to N)
    unique_ids = unique(masks_filtered(masks_filtered > 0));
    final_masks = zeros(size(masks_filtered), 'uint16');
    
    for i = 1:length(unique_ids)
        final_masks(masks_filtered == unique_ids(i)) = i;
    end
    
    masks_filtered = final_masks;
    numObjects = max(masks_filtered(:));
    fprintf('Found %d nuclei\n', numObjects);
 
    if numObjects == 0
        warning('No nuclei found — check thresholds or area filter.');
        out{1} = {};
        out{2} = [];
        out{3} = [];
        return
    end
    
    % Use the dimensions of your final label matrix
    [rows, cols] = size(masks_filtered);                  
    
    % Initialize arrays to hold the outputs
    individual_masks = cell(numObjects, 1);
    nuclear_intensities = zeros(numObjects, 1);
    
    % regionprops automatically groups pixels by their unique integer ID
    props = regionprops(masks_filtered, 'PixelIdxList');
    for i = 1:numObjects
        % 1) Build the individual binary mask for object 'i'
        mask = false(rows, cols);
        mask(props(i).PixelIdxList) = true;
        individual_masks{i} = mask;
        
        % 2) Calculate nuclear intensity in this mask using the original image
        values = double(inputImage_cond(props(i).PixelIdxList));
        % Mean intensity (excluding exact zeros)
        nuclear_intensities(i) = mean(values(values ~= 0));   
    end
    
    % 3) The combined image with all masks and unique IDs
    combined_image = masks_filtered;
    
    % Assign to the final output array
    out{1} = individual_masks;
    out{2} = nuclear_intensities;
    out{3} = combined_image;
end