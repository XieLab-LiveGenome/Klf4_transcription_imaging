function nuc = nuc_motion_MS2(MIP_image6d, channel1, C_in, varargin)
%NUC_MOTION_MS2  Per-nucleus rigid translation from BLOCK AVERAGES by image
%                registration (phase correlation) - no per-frame segmentation.
%
%   nuc = NUC_MOTION_MS2(MIP_image6d, channel1, C_in)
%   nuc = NUC_MOTION_MS2(..., 'Name', value, ...)
%
%   Why this exists: taking the centroid of a Cellpose mask as the nucleus
%   position makes the motion estimate inherit every shape fluctuation of the
%   segmentation. On weak diffuse nuclear signal the mask area swings by ~25%
%   frame to frame, which moves the centroid by tens of pixels even when the
%   nucleus has not moved. Registration compares IMAGE CONTENT instead, so it
%   does not care where the mask boundary landed.
%
%   How it works:
%     1. Segment ONCE on the whole-movie time average (clean, high SNR) to get
%        the nucleus ROIs and to assign each TS in C_in to a nucleus.
%     2. Average the movie in blocks of block_size frames.
%     3. For each nucleus, phase-correlate its ROI in block b against block 1
%        -> cumulative [dy dx] per block.
%     4. Interpolate the block shifts onto every frame.
%
%   Options (defaults in brackets):
%       block_size  [10]     frames per block
%       pad         [100]    ROI padding around the nucleus bbox, px
%       sigma       [4]      Gaussian blur before registration, px. Keep this
%                            SMALL - heavy blur (e.g. 30) biases the shift
%                            badly (~2 px error vs ~0.1 px at sigma 4).
%       max_step    [150]    reject a block shift that jumps more than this
%                            from the previous block (holds previous instead)
%       interp      ['pchip'] interpolation of block shifts onto frames
%       preproc     ['none'] 'none' or 'snr' (apply SNR_inc2 before segmenting)
%       snr_mul     [0.5]    mul passed to SNR_inc2
%       cell_D_seg  [160]    cellpose diameter for the one-off segmentation
%       verbose     [true]
%
%   Output struct nuc - same field names track_TS_spots expects, so this is a
%   drop-in replacement for track_nuclei_MS2 when you only need motion:
%       .shifts        (T x 2 x n_stable) cumulative [dy dx] from frame 1
%       .centroids     (T x 2 x n_stable) ROI centre + shift, [y x]
%       .spot_gid      (vis x 1) nucleus ID per TS
%       .spot_stable   (vis x 1) stable index per TS
%       .stable_ids .n_stable .linked_masks (single reference mask, not per frame)
%       .block_shifts  (n_blocks x 2 x n_stable) raw per-block estimate
%       .block_mids    (n_blocks x 1) mid-frame of each block
%       .resp          (n_blocks x n_stable) correlation peak height (quality)

opts = struct('block_size',10, 'pad',100, 'sigma',4, 'max_step',150, ...
              'interp','pchip', 'preproc','none', 'snr_mul',0.5, ...
              'cell_D_seg',160, 'verbose',true);
opts = local_parse_opts(opts, varargin);

T              = size(MIP_image6d, 2);
[m_img, n_img] = size(squeeze(MIP_image6d(1,1,1,channel1,:,:)));

%% ===== 1. ONE-OFF SEGMENTATION ON THE TIME AVERAGE =====
if opts.verbose, fprintf('\n=== NUCLEUS MOTION (block registration) ===\n'); end

I_avg_all = zeros(m_img, n_img);
for f = 1:T
    I_avg_all = I_avg_all + double(squeeze(MIP_image6d(1,f,1,channel1,:,:)));
end
I_avg_all = I_avg_all / T;

I_seg = I_avg_all;
if strcmpi(opts.preproc,'snr')
    I_seg = SNR_inc2(I_seg, opts.snr_mul);
end

seg_OUT   = cellpose_seg_MS2(I_seg, opts.cell_D_seg);
ref_labels = seg_OUT{3};
n_stable   = double(max(ref_labels(:)));

if n_stable == 0
    error('nuc_motion_MS2:noNuclei', ...
        'No nuclei segmented on the time average - check cell_D_seg / preproc.');
end
stable_ids = (1:n_stable)';

props = regionprops(ref_labels, 'Centroid', 'BoundingBox');
rois  = zeros(n_stable, 4);   % [r0 r1 c0 c1]
cent0 = zeros(n_stable, 2);   % [y x]
for s = 1:n_stable
    bb = props(s).BoundingBox;                 % [x y w h]
    c0 = max(1,      floor(bb(1)) - opts.pad);
    r0 = max(1,      floor(bb(2)) - opts.pad);
    c1 = min(n_img,  ceil(bb(1)+bb(3)) + opts.pad);
    r1 = min(m_img,  ceil(bb(2)+bb(4)) + opts.pad);
    rois(s,:)  = [r0 r1 c0 c1];
    cent0(s,:) = [props(s).Centroid(2) props(s).Centroid(1)];
    if opts.verbose
        fprintf('  nucleus %d: centroid (%.0f, %.0f), ROI %dx%d px\n', ...
            s, cent0(s,1), cent0(s,2), r1-r0+1, c1-c0+1);
    end
end

%% ===== 2-3. BLOCK AVERAGES -> PER-NUCLEUS PHASE CORRELATION =====
n_blocks     = ceil(T / opts.block_size);
block_mids   = zeros(n_blocks,1);
block_shifts = zeros(n_blocks, 2, n_stable);
resp         = zeros(n_blocks, n_stable);

ref_roi = cell(n_stable,1);      % reference ROI = first block average

for b = 1:n_blocks
    f0 = (b-1)*opts.block_size + 1;
    f1 = min(b*opts.block_size, T);
    block_mids(b) = (f0 + f1) / 2;

    I_blk = zeros(m_img, n_img);
    for f = f0:f1
        I_blk = I_blk + double(squeeze(MIP_image6d(1,f,1,channel1,:,:)));
    end
    I_blk = I_blk / (f1 - f0 + 1);

    for s = 1:n_stable
        r = rois(s,:);
        cur = imgaussfilt(I_blk(r(1):r(2), r(3):r(4)), opts.sigma);

        if b == 1
            ref_roi{s}          = cur;
            block_shifts(b,:,s) = [0 0];
            resp(b,s)           = 1;
            continue
        end

        [dy, dx, pk] = local_phasecorr(ref_roi{s}, cur);

        % reject implausible jumps: hold the previous block instead
        prev = squeeze(block_shifts(b-1,:,s));
        if hypot(dy - prev(1), dx - prev(2)) > opts.max_step
            if opts.verbose
                fprintf('  ! block %d nuc %d: rejected jump (%.0f, %.0f) px, holding previous\n', ...
                    b, s, dy, dx);
            end
            dy = prev(1); dx = prev(2);
        end

        block_shifts(b,:,s) = [dy dx];
        resp(b,s)           = pk;
    end

    if opts.verbose
        fprintf('  block %d/%d (frames %d-%d) done\n', b, n_blocks, f0, f1);
    end
end

%% ===== 4. INTERPOLATE BLOCK SHIFTS ONTO EVERY FRAME =====
nuc_shifts    = zeros(T, 2, n_stable);
nuc_centroids = zeros(T, 2, n_stable);

for s = 1:n_stable
    for k = 1:2
        if n_blocks >= 2
            nuc_shifts(:,k,s) = interp1(block_mids, block_shifts(:,k,s), ...
                                        (1:T)', opts.interp, 'extrap');
        else
            nuc_shifts(:,k,s) = block_shifts(1,k,s);
        end
    end
    % express relative to frame 1 so .shifts means "moved since frame 1"
    nuc_shifts(:,1,s) = nuc_shifts(:,1,s) - nuc_shifts(1,1,s);
    nuc_shifts(:,2,s) = nuc_shifts(:,2,s) - nuc_shifts(1,2,s);

    nuc_centroids(:,1,s) = cent0(s,1) + nuc_shifts(:,1,s);
    nuc_centroids(:,2,s) = cent0(s,2) + nuc_shifts(:,2,s);
end

%% ===== ASSIGN TS SPOTS TO NUCLEI (on the reference mask) =====
vis             = size(C_in,1);
spot_nuc_gid    = zeros(vis,1);
spot_nuc_stable = zeros(vis,1);

for p = 1:vis
    cx = max(1, min(round(C_in(p,1)), n_img));
    cy = max(1, min(round(C_in(p,2)), m_img));
    gid = double(ref_labels(cy, cx));

    if gid == 0                              % nearest nucleus centroid
        best = Inf;
        for s = 1:n_stable
            d = hypot(cent0(s,2)-cx, cent0(s,1)-cy);
            if d < best, best = d; gid = s; end
        end
        if opts.verbose
            fprintf('  spot %d outside every mask - assigned nearest nucleus %d\n', p, gid);
        end
    end

    spot_nuc_gid(p)    = gid;
    spot_nuc_stable(p) = gid;                % 1:1 here, ids are contiguous
    if opts.verbose
        fprintf('Spot %d -> nucleus %d\n', p, gid);
    end
end

if opts.verbose
    fprintf('Blocks: %d (size %d) | nuclei: %d | frames: %d\n', ...
        n_blocks, opts.block_size, n_stable, T);
    for s = 1:n_stable
        tot = hypot(nuc_shifts(end,1,s), nuc_shifts(end,2,s));
        fprintf('  nucleus %d net drift frame 1 -> %d: %.1f px\n', s, T, tot);
    end
end

%% ===== PACK OUTPUT (field names match track_nuclei_MS2) =====
nuc.shifts       = nuc_shifts;
nuc.centroids    = nuc_centroids;
nuc.spot_gid     = spot_nuc_gid;
nuc.spot_stable  = spot_nuc_stable;
nuc.stable_ids   = stable_ids;
nuc.n_stable     = n_stable;
nuc.linked_masks = ref_labels;          % single reference mask, not per frame
nuc.block_shifts = block_shifts;
nuc.block_mids   = block_mids;
nuc.rois         = rois;
nuc.resp         = resp;
nuc.size         = [m_img n_img T];
nuc.opts         = opts;
nuc.method       = 'block_registration';

end


% =====================================================================
function [dy, dx, pk] = local_phasecorr(fixed, moving)
%LOCAL_PHASECORR  Sub-pixel translation of MOVING relative to FIXED.
%   Positive dy/dx = content moved toward larger row/column index.
%   Plain FFT phase correlation - no toolbox version dependency.

fixed  = double(fixed);
moving = double(moving);

% zero-mean + separable Hann window to suppress edge effects
fixed  = fixed  - mean(fixed(:));
moving = moving - mean(moving(:));
[r, c] = size(fixed);
w = hann_local(r) * hann_local(c)';
fixed  = fixed  .* w;
moving = moving .* w;

F = fft2(fixed);
M = fft2(moving);
R = conj(F) .* M;                       % peak lands at +shift of MOVING
R = R ./ max(abs(R), eps);
corr = real(ifft2(R));

[pk, idx]  = max(corr(:));
[py, px]   = ind2sub([r c], idx);

% parabolic sub-pixel refinement
dy = local_subpix(corr, py, px, 1, r);
dx = local_subpix(corr, py, px, 2, c);

% wrap to signed shift
if dy > r/2, dy = dy - r; end
if dx > c/2, dx = dx - c; end

end


function d = local_subpix(corr, py, px, dim, n)
ip = @(i) mod(i-1, n) + 1;
if dim == 1
    a = corr(ip(py-1), px);  b = corr(py, px);  cc = corr(ip(py+1), px);
    base = py - 1;
else
    a = corr(py, ip(px-1));  b = corr(py, px);  cc = corr(py, ip(px+1));
    base = px - 1;
end
den = (a - 2*b + cc);
if den == 0
    d = base;
else
    d = base + 0.5*(a - cc)/den;
end
end


function w = hann_local(n)
w = 0.5*(1 - cos(2*pi*(0:n-1)'/(n-1)));
end


% =====================================================================
function opts = local_parse_opts(opts, args)
if isempty(args), return; end
if numel(args) == 1 && isstruct(args{1})
    f = fieldnames(args{1});
    for i = 1:numel(f), opts.(f{i}) = args{1}.(f{i}); end
    return
end
if mod(numel(args), 2) ~= 0
    error('Options must be name/value pairs or a single struct.');
end
for i = 1:2:numel(args)
    if ~isfield(opts, args{i}), error('Unknown option "%s".', args{i}); end
    opts.(args{i}) = args{i+1};
end
end
