function out = cellpose_seg(inputImage_cond,cell_D)

M = median(inputImage_cond(:));

% Normalize and adjust
img_uint16 = uint16(inputImage_cond / max(inputImage_cond(:)) * 65535);

thresholds = M;
gammas = 0.6;
saturations = min(5*M,65535);

img_uint16 = adjust_multichannel_image(img_uint16, thresholds, gammas, saturations); % improve brightness by saturating signal

% Replace medfilt2 with imgaussfilt for smoother flow-field calculations
img_uint16 = imgaussfilt(img_uint16, 3);

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

% Run segmentation with auto-diameter (0)
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
    
    % Find which Cell IDs meet your 3000 pixel threshold
    valid_ids = find([stats.Area] >= 10000); 
    
    % Keep only the pixels that belong to those valid IDs
    masks_filtered = masks;
    masks_filtered(~ismember(masks, valid_ids)) = 0;

    % 2. GENTLE dilation (Without merging touching cells)
    dil_radius = 15; 
    se = strel('disk', dil_radius);
    
    trimmed_boundaries = imdilate(masks_filtered > 0, se);
    masks_filtered(~trimmed_boundaries) = 0;

    unique_ids = unique(masks_filtered(masks_filtered > 0));

    final_masks = zeros(size(masks_filtered), 'uint16');

    for i = 1:length(unique_ids)
        final_masks(masks_filtered == unique_ids(i)) = i;
    end
    
    masks_filtered = final_masks;

numObjects = max(masks_filtered(:));

fprintf('Found %d nuclei\n', numObjects);

    numObjects = max(masks_filtered(:));
 
    
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
    % This is much faster than running (masks_filtered == i) inside a loop
    props = regionprops(masks_filtered, 'PixelIdxList');

    for i = 1:numObjects
        % 1) Build the individual binary mask for object 'i'
        mask = false(rows, cols);
        mask(props(i).PixelIdxList) = true;
        individual_masks{i} = mask;
        
        % 2) Calculate nuclear intensity in this mask using the original image
        values = double(inputImage_cond(props(i).PixelIdxList));
        % Mean intensity (excluding exact zeros, per your original logic)
        nuclear_intensities(i) = mean(values(values ~= 0));   
    end

    % 3) The combined image with all masks and unique IDs
    % masks_filtered is already exactly this: a label matrix where 
    % background = 0, and cells = 1, 2, 3... N.
    combined_image = masks_filtered;


    % Assign to the final output array
    out{1} = individual_masks;
    out{2} = nuclear_intensities;
    out{3} = combined_image;

end