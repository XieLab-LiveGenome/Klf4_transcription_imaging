function nuc = track_nuclei_MS2(MIP_image6d, channel1, C_in, varargin)
%TRACK_NUCLEI_MS2  Consistent nucleus detection, frame-by-frame registration
%                  and TS-to-nucleus assignment.
%
%   nuc = TRACK_NUCLEI_MS2(MIP_image6d, channel1, C_in)
%   nuc = TRACK_NUCLEI_MS2(..., 'Name', value, ...)
%   nuc = TRACK_NUCLEI_MS2(..., optsStruct)
%
%   Step 1  block-averages the movie and runs Cellpose on each block, linking
%           blocks into a global nucleus registry (Hungarian matching).
%   Step 2  keeps nuclei present in at least min_presence of the blocks.
%   Step 3  segments every frame and relabels it with the stable global IDs.
%   Step 4  assigns each TS in C_in to a stable nucleus (pass [] to skip).
%
%   Options (defaults in brackets):
%       block_size    [10]    frames per block
%       min_presence  [0.75]  fraction of blocks a nucleus must appear in
%       cell_D_seg    [160]   cellpose diameter
%       max_drift_px  [200]   max centroid drift for matching
%       alpha_ema     [0.3]   centroid EMA smoothing
%       verbose       [true]  print the consistency report and progress
%
%   Output struct nuc:
%       .stable_ids       global IDs of the consistently detected nuclei
%       .n_stable         number of stable nuclei
%       .linked_masks     (m x n x T) uint16 label matrix with stable IDs
%       .centroids        (T x 2 x n_stable) centroid [y, x] per frame
%       .shifts           (T x 2 x n_stable) cumulative [dy, dx] from frame 1
%       .gid_to_stable    containers.Map from global ID to stable index
%       .presence         (n_global x n_blocks) detection table
%       .spot_gid         (vis x 1) global nucleus ID per TS
%       .spot_stable      (vis x 1) stable index per TS (0 if unassigned)
%       .individual_masks binary masks from the last block segmentation
%       .size             [m_img n_img T]

if nargin < 3, C_in = []; end

opts = struct('block_size',10, 'min_presence',0.75, 'cell_D_seg',160, ...
              'max_drift_px',200, 'alpha_ema',0.3, 'verbose',true);
opts = local_parse_opts(opts, varargin);

total_frames   = size(MIP_image6d, 2);
[m_img, n_img] = size(squeeze(MIP_image6d(1,1,1,channel1,:,:)));
n_blocks       = ceil(total_frames / opts.block_size);

%% === BLOCK SEGMENTATION -> GLOBAL NUCLEUS REGISTRY ===
global_centroids = [];
global_last_seen = [];
global_n         = 0;

presence      = false(100, n_blocks);   % preallocate generously
block_seg     = cell(n_blocks, 1);      % label matrices
block_id_maps = cell(n_blocks, 1);      % local-to-global ID maps
block_cents   = cell(n_blocks, 1);      % centroids
seg_OUT       = {};

for b = 1:n_blocks

    f_start    = (b - 1) * opts.block_size + 1;
    f_end      = min(b * opts.block_size, total_frames);
    n_in_block = f_end - f_start + 1;

    % --- time average ---
    I_sum = zeros(m_img, n_img);
    for f = f_start:f_end
        I_sum = I_sum + double(squeeze(MIP_image6d(1, f, 1, channel1, :, :)));
    end
    I_avg = I_sum / n_in_block;

    % --- segment ---
    seg_OUT    = cellpose_seg_MS2(I_avg, opts.cell_D_seg);
    blk_labels = seg_OUT{3};
    n_nuc_blk  = max(blk_labels(:));

    if n_nuc_blk == 0
        block_seg{b}     = zeros(m_img, n_img, 'uint16');
        block_id_maps{b} = [];
        block_cents{b}   = [];
        continue
    end

    blk_props     = regionprops(blk_labels, 'Centroid');
    blk_centroids = zeros(n_nuc_blk, 2);
    for i = 1:n_nuc_blk
        if ~isempty(blk_props(i).Centroid)
            blk_centroids(i,:) = blk_props(i).Centroid;
        end
    end

    % --- match to the global registry (Hungarian) ---
    id_map_blk = zeros(n_nuc_blk, 1);

    if global_n == 0
        for i = 1:n_nuc_blk
            global_n      = global_n + 1;
            id_map_blk(i) = global_n;
        end
        global_centroids = blk_centroids;
        global_last_seen = ones(n_nuc_blk, 1) * b;
    else
        cost = ones(n_nuc_blk, global_n) * opts.max_drift_px * 2;
        for i = 1:n_nuc_blk
            for jj = 1:global_n
                d       = sqrt(sum((blk_centroids(i,:) - global_centroids(jj,:)).^2));
                penalty = (b - global_last_seen(jj)) * 10;
                cost(i, jj) = d + penalty;
            end
        end

        M = matchpairs(cost, opts.max_drift_px);

        for row = 1:size(M, 1)
            blk_id  = M(row, 1);
            glob_id = M(row, 2);
            id_map_blk(blk_id) = glob_id;
            global_centroids(glob_id,:) = ...
                (1 - opts.alpha_ema) * global_centroids(glob_id,:) + ...
                 opts.alpha_ema * blk_centroids(blk_id,:);
            global_last_seen(glob_id) = b;
        end

        unmatched = find(id_map_blk == 0);
        for u = 1:length(unmatched)
            global_n = global_n + 1;
            id_map_blk(unmatched(u)) = global_n;
            global_centroids(global_n,:) = blk_centroids(unmatched(u),:);
            global_last_seen(global_n)   = b;
        end
    end

    if global_n > size(presence, 1)
        presence(end+1:global_n, :) = false;
    end

    for i = 1:n_nuc_blk
        presence(id_map_blk(i), b) = true;
    end

    block_seg{b}     = blk_labels;
    block_id_maps{b} = id_map_blk;
    block_cents{b}   = blk_centroids;
end

presence = presence(1:global_n, :);

%% === IDENTIFY STABLE NUCLEI ===
frac_present = sum(presence, 2) / n_blocks;
stable_ids   = find(frac_present >= opts.min_presence);
n_stable     = length(stable_ids);

if opts.verbose
    fprintf('\n=== NUCLEUS CONSISTENCY REPORT ===\n');
    fprintf('Total unique nuclei detected: %d\n', global_n);
    fprintf('Minimum presence threshold:   %.0f%% (%d/%d blocks)\n', ...
        opts.min_presence*100, ceil(opts.min_presence*n_blocks), n_blocks);
    fprintf('Stable nuclei: %d  (IDs: %s)\n', n_stable, num2str(stable_ids'));
    fprintf('\nPer-nucleus presence:\n');
    for i = 1:global_n
        marker = '  ';
        if ismember(i, stable_ids), marker = '>>'; end
        fprintf('  %s Nuc %2d: %2d/%d blocks (%.0f%%)\n', ...
            marker, i, sum(presence(i,:)), n_blocks, frac_present(i)*100);
    end
end

if n_stable == 0
    warning('track_nuclei_MS2:noStableNuclei', ...
        'No nucleus met the presence threshold - lower min_presence or check cell_D_seg.');
end

%% === FRAME-BY-FRAME REGISTRATION ===
if opts.verbose
    fprintf('\nRunning frame-by-frame registration for %d stable nuclei...\n', n_stable);
end

linked_masks  = zeros(m_img, n_img, total_frames, 'uint16');
nuc_centroids = NaN(total_frames, 2, n_stable);   % [y, x]

% reference centroid per frame, interpolated from the block centroids
ref_cent_per_frame = zeros(total_frames, 2, n_stable);

for s = 1:n_stable
    gid = stable_ids(s);

    block_frames = [];
    block_xy     = [];

    for b = 1:n_blocks
        if presence(gid, b) && ~isempty(block_id_maps{b})
            local_idx = find(block_id_maps{b} == gid);
            if ~isempty(local_idx)
                f_mid = round(((b-1)*opts.block_size + 1 + ...
                    min(b*opts.block_size, total_frames)) / 2);
                block_frames(end+1) = f_mid;                      %#ok<AGROW>
                block_xy(end+1,:)   = block_cents{b}(local_idx,:); %#ok<AGROW>
            end
        end
    end

    if length(block_frames) >= 2
        ref_cent_per_frame(:, 1, s) = interp1(block_frames, block_xy(:,1), ...
            1:total_frames, 'pchip', 'extrap');
        ref_cent_per_frame(:, 2, s) = interp1(block_frames, block_xy(:,2), ...
            1:total_frames, 'pchip', 'extrap');
    elseif length(block_frames) == 1
        ref_cent_per_frame(:, 1, s) = block_xy(1, 1);
        ref_cent_per_frame(:, 2, s) = block_xy(1, 2);
    end
end

for k = 1:total_frames

    I_raw = squeeze(MIP_image6d(1, k, 1, channel1, :, :));

    seg_k        = cellpose_seg_MS2(I_raw, opts.cell_D_seg);
    frame_labels = seg_k{3};
    frame_ids    = unique(frame_labels(frame_labels > 0));
    frame_props  = regionprops(frame_labels, 'Centroid');

    frame_linked = zeros(m_img, n_img, 'uint16');

    for s = 1:n_stable
        gid = stable_ids(s);
        rcx = ref_cent_per_frame(k, 1, s);  % expected x
        rcy = ref_cent_per_frame(k, 2, s);  % expected y

        ry = max(1, min(round(rcy), m_img));
        rx = max(1, min(round(rcx), n_img));

        matched_id = frame_labels(ry, rx);

        % fallback: nearest centroid search
        if matched_id == 0
            best_dist = Inf;
            best_id   = 0;
            for jj = 1:length(frame_ids)
                cid = frame_ids(jj);
                if cid <= length(frame_props) && ~isempty(frame_props(cid).Centroid)
                    d = sqrt(sum((frame_props(cid).Centroid - [rcx rcy]).^2));
                    if d < best_dist && d < opts.max_drift_px
                        best_dist = d;
                        best_id   = cid;
                    end
                end
            end
            matched_id = best_id;
        end

        if matched_id > 0
            frame_linked(frame_labels == matched_id) = gid;

            if matched_id <= length(frame_props) && ~isempty(frame_props(matched_id).Centroid)
                nuc_centroids(k, 1, s) = frame_props(matched_id).Centroid(2);  % y (row)
                nuc_centroids(k, 2, s) = frame_props(matched_id).Centroid(1);  % x (col)
            end
        end
    end

    linked_masks(:,:,k) = frame_linked;

    if opts.verbose && (mod(k, 20) == 0 || k == total_frames)
        fprintf('  Frame %d/%d done\n', k, total_frames);
    end
end

% --- fill NaN gaps in the centroids ---
for s = 1:n_stable
    for dim = 1:2
        v     = nuc_centroids(:, dim, s);
        valid = find(~isnan(v));
        if length(valid) >= 2
            nuc_centroids(:, dim, s) = interp1(valid, v(valid), ...
                1:total_frames, 'linear', 'extrap');
        elseif length(valid) == 1
            nuc_centroids(:, dim, s) = v(valid);
        end
    end
end

% --- cumulative rigid-body shifts relative to frame 1 ---
nuc_shifts = zeros(total_frames, 2, n_stable);   % [dy, dx]
for s = 1:n_stable
    for k = 2:total_frames
        nuc_shifts(k, 1, s) = nuc_centroids(k, 1, s) - nuc_centroids(1, 1, s);
        nuc_shifts(k, 2, s) = nuc_centroids(k, 2, s) - nuc_centroids(1, 2, s);
    end
end

% --- lookup: stable index from global ID ---
gid_to_stable = containers.Map('KeyType', 'int32', 'ValueType', 'int32');
for s = 1:n_stable
    gid_to_stable(stable_ids(s)) = s;
end

%% === ASSIGN TS SPOTS TO STABLE NUCLEI ===
vis             = size(C_in, 1);
spot_nuc_gid    = zeros(vis, 1);
spot_nuc_stable = zeros(vis, 1);

for p = 1:vis
    cx = round(C_in(p, 1));
    cy = round(C_in(p, 2));
    cy = max(1, min(cy, m_img));
    cx = max(1, min(cx, n_img));

    gid = linked_masks(cy, cx, 1);

    if gid == 0
        best_dist = Inf;
        for s = 1:n_stable
            d = sqrt((nuc_centroids(1,2,s) - cx)^2 + (nuc_centroids(1,1,s) - cy)^2);
            if d < best_dist
                best_dist = d;
                gid       = stable_ids(s);
            end
        end
    end

    spot_nuc_gid(p) = gid;
    if gid > 0 && isKey(gid_to_stable, gid)
        spot_nuc_stable(p) = gid_to_stable(gid);
    else
        warning('track_nuclei_MS2:spotOutsideNucleus', ...
            'Spot %d at [%.0f, %.0f] not inside a stable nucleus (gid=%d)', p, cx, cy, gid);
        spot_nuc_stable(p) = 0;
    end
    if opts.verbose
        fprintf('Spot %d -> nucleus %d (stable index %d)\n', p, gid, spot_nuc_stable(p));
    end
end

if opts.verbose
    fprintf('\n=== REGISTRATION COMPLETE ===\n');
    fprintf('Stable nuclei: %d (IDs: %s)\n', n_stable, num2str(stable_ids'));
    fprintf('Frames registered: %d\n', total_frames);
    fprintf('Spots assigned: %d/%d to stable nuclei\n', sum(spot_nuc_stable > 0), vis);
end

%% === PACK OUTPUT ===
nuc.stable_ids       = stable_ids;
nuc.n_stable         = n_stable;
nuc.linked_masks     = linked_masks;
nuc.centroids        = nuc_centroids;
nuc.shifts           = nuc_shifts;
nuc.gid_to_stable    = gid_to_stable;
nuc.presence         = presence;
nuc.spot_gid         = spot_nuc_gid;
nuc.spot_stable      = spot_nuc_stable;
nuc.block_seg        = block_seg;
nuc.block_id_maps    = block_id_maps;
nuc.size             = [m_img n_img total_frames];
nuc.opts             = opts;
if ~isempty(seg_OUT)
    nuc.individual_masks = seg_OUT{1,1};
else
    nuc.individual_masks = {};
end

end


% =====================================================================
function opts = local_parse_opts(opts, args)
%LOCAL_PARSE_OPTS  Accept either a single options struct or name/value pairs.

if isempty(args), return; end

if numel(args) == 1 && isstruct(args{1})
    f = fieldnames(args{1});
    for i = 1:numel(f)
        opts.(f{i}) = args{1}.(f{i});
    end
    return
end

if mod(numel(args), 2) ~= 0
    error('Options must be name/value pairs or a single struct.');
end
for i = 1:2:numel(args)
    if ~isfield(opts, args{i})
        error('Unknown option "%s".', args{i});
    end
    opts.(args{i}) = args{i+1};
end

end
