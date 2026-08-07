function out = cellpose_seg_MS2(inputImage_cond, cell_D, preproc, snr_mul)
% CELLPOSE_SEG  Segment nuclei from SIM images with weak diffuse nuclear signal.
%
%   Preprocessing: large Gaussian blur → percentile norm → gamma → CLAHE
%   The blur dilutes bright MS2 puncta while preserving the diffuse nuclear halo.
%
%   Inputs:
%       inputImage_cond - 2D image (any numeric type)
%       cell_D          - expected nucleus diameter in downsampled pixels
%       preproc         - (optional) 'none' (default) or 'snr' to run SNR_inc2
%                         before the standard preprocessing
%       snr_mul         - (optional) mul passed to SNR_inc2 (default 0.5)
%
%   Outputs:
%       out{1} - cell array of individual binary masks
%       out{2} - mean nuclear intensity per object
%       out{3} - combined label matrix (uint16)

    if nargin < 3 || isempty(preproc), preproc = 'none'; end
    if nargin < 4 || isempty(snr_mul),  snr_mul = 0.5;    end

    img = double(inputImage_cond);

    % === STEP 0 (OPTIONAL): SNR_inc2 CONTRAST BOOST ===
    % Brightens pixels above snr_mul*max and dims the rest, then median filters.
    % Off by default: on weak diffuse nuclear signal it mainly amplifies the
    % MS2 puncta, which STEP 1 then blurs away again.
    if strcmpi(preproc, 'snr')
        img = SNR_inc2(img, snr_mul);
    end

    % === STEP 1: LARGE GAUSSIAN BLUR ===
    % MS2 spots (~5 px) get diluted into background; diffuse nuclear 
    % signal (~200+ px) is preserved. No morphological opening needed.
    img_smooth = imgaussfilt(img, 30);

    % === STEP 2: PERCENTILE NORMALIZATION ===
    % Use 1st-99.5th to avoid any residual bright pixels setting ceiling
    p_low  = prctile(img_smooth(:), 1);
    p_high = prctile(img_smooth(:), 99.5);
    img_norm = (img_smooth - p_low) / max(p_high - p_low, eps);
    img_norm = max(0, min(1, img_norm));

    % === STEP 3: AGGRESSIVE GAMMA ===
    % Gamma 0.3 boosts dim nuclear signal into the upper range
    img_gamma = img_norm .^ 0.3;

    % === STEP 4: CLAHE ===
    % Local adaptive contrast stretches dim nuclei to fill local dynamic range
    img_uint16 = uint16(img_gamma * 65535);
    img_uint16 = adapthisteq(img_uint16, ...
        'NumTiles',     [8 8], ...
        'ClipLimit',    0.03, ...
        'Distribution', 'rayleigh');

    % === STEP 5: LIGHT SMOOTH FOR CELLPOSE FLOW FIELDS ===
    img_uint16 = imgaussfilt(img_uint16, 3);

    % --- DOWNSAMPLING ---
    scale_factor = 0.3;
    img_resized = imresize(img_uint16, scale_factor, 'bicubic');

    % --- CELLPOSE SEGMENTATION (Python Environment) ---
cp_models = py.importlib.import_module('cellpose.models');
np        = py.importlib.import_module('numpy');
img_np    = np.array(img_resized);

cp_ver   = char(py.importlib.metadata.version('cellpose'));
cp_major = sscanf(cp_ver, '%d', 1);

if cp_major < 4
    model   = cp_models.Cellpose(pyargs('gpu', false, 'model_type', 'nuclei'));
    results = model.eval(img_np, pyargs('diameter', cell_D, ...
                                        'channels', py.list({0, 0}), ...
                                        'augment',  false));
else
    warning('cellpose_seg_MS2:v4', ...
        ['Cellpose %s detected (Cellpose-SAM). Masks may differ from ' ...
         'the v3 results used for published analyses.'], cp_ver);
    model   = cp_models.CellposeModel(pyargs('gpu', false));
    results = model.eval(img_np, pyargs('diameter', cell_D));
end

masks_resized = uint16(double(py.numpy.array(results{1})));

    % --- UPSAMPLING & FILTERING ---
    original_size = size(img_uint16);
    masks = imresize(masks_resized, original_size, 'nearest');

    % 1. FILTER BY SIZE
    stats = regionprops(masks, 'Area');
    valid_ids = find([stats.Area] >= 10000);
    masks_filtered = masks;
    masks_filtered(~ismember(masks, valid_ids)) = 0;

    % 2. GENTLE DILATION (without merging touching cells)
    dil_radius = 15;
    se = strel('disk', dil_radius);
    trimmed_boundaries = imdilate(masks_filtered > 0, se);
    masks_filtered(~trimmed_boundaries) = 0;

    % 3. RELABEL TO CONTIGUOUS IDs
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

    % --- BUILD OUTPUT MASKS AND INTENSITIES ---
    [rows, cols] = size(masks_filtered);
    individual_masks = cell(numObjects, 1);
    nuclear_intensities = zeros(numObjects, 1);

    props = regionprops(masks_filtered, 'PixelIdxList');

    for i = 1:numObjects
        mask = false(rows, cols);
        mask(props(i).PixelIdxList) = true;
        individual_masks{i} = mask;

        values = double(inputImage_cond(props(i).PixelIdxList));
        nuclear_intensities(i) = mean(values(values ~= 0));
    end

    combined_image = masks_filtered;

    out{1} = individual_masks;
    out{2} = nuclear_intensities;
    out{3} = combined_image;

end