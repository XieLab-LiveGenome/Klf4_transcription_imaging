function trk = track_TS_spots(MIP_image6d, channel1, C_in, nuc, motion_correct, varargin)
%TRACK_TS_SPOTS  Track MS2 transcription sites frame-by-frame, optionally
%                correcting the search seed for host-nucleus motion.
%
%   trk = TRACK_TS_SPOTS(MIP_image6d, channel1, C_in, nuc, motion_correct)
%   trk = TRACK_TS_SPOTS(..., 'Name', value, ...)
%   trk = TRACK_TS_SPOTS(..., optsStruct)
%
%   motion_correct = false : seed every frame's search from the previous
%                            fitted position - independent tracking, exactly
%                            as before. nuc may be [] in this case.
%   motion_correct = true  : before searching frame k, shift the seed by how
%                            far the host nucleus moved between frame k-1 and
%                            frame k, and carry that same drift into frames
%                            where no spot is found. Requires a populated nuc
%                            struct from track_nuclei_MS2 (needs nuc.shifts and
%                            nuc.spot_stable).
%
%   In both modes track_spot2D re-detects the actual brightest spot within the
%   search radius, and the intensity is sampled at the FITTED position. Turning
%   correction on therefore changes WHERE the search is centred, not WHAT is
%   measured - it keeps the search window riding along with a drifting nucleus
%   so a fast-moving or transiently dark spot is not lost.
%
%   Options (defaults in brackets):
%       start_f [1]     first frame (nucleus shifts assume registration on 1:zs,
%                       so motion_correct requires start_f == 1)
%       Rp      [80]    search radius after a hit, in px
%       Rn      [150]   search radius after a miss, in px
%       m_MS2   [0.3]   threshold multiplier passed to track_spot2D
%       R_c     [5]     radius of the intensity mask (maskavg)
%       show    [true]  draw the per-frame tracking overlay
%       fig     []      figure handle to reuse (created if empty and show=true)
%
%   Output struct trk:
%       .C_cent    (zs+1 x 2 x vis)  [x y]; row k+1 holds the frame-k result
%       .I_t       (zs x vis)        raw intensity at the fitted spot
%       .sd        (zs x vis)        per-frame intensity SD
%       .R_spot_n  (zs x vis)        search radius used each frame
%       .spot_size (zs x vis)        fitted spot size (kept for compatibility)
%       .motion_correct              the flag actually used
%       .corrected_spots             logical (vis x 1): was correction applied

opts = struct('start_f',1, 'Rp',80, 'Rn',150, 'm_MS2',0.3, 'R_c',5, ...
              'show',true, 'fig',[]);
opts = local_parse_opts(opts, varargin);

zs  = size(MIP_image6d, 2);
vis = size(C_in, 1);

% --- validate the motion-correction inputs up front ---
use_shift = false(vis,1);
if motion_correct
    if isempty(nuc) || ~isfield(nuc,'shifts') || isempty(nuc.shifts)
        error('track_TS_spots:noNucleusData', ...
            ['motion_correct = true but no nucleus data was supplied. ' ...
             'Run track_nuclei_MS2 first (run_nuclei = true) and pass its output.']);
    end
    if opts.start_f ~= 1
        error('track_TS_spots:startFrameMismatch', ...
            ['motion_correct = true requires start_f = 1 so the nucleus shifts ' ...
             '(registered on frames 1:zs) line up with the tracked frames.']);
    end
    nuc_shifts = nuc.shifts;                 % (zs x 2 x n_stable) [dy dx]
    spot_stable = nuc.spot_stable;           % (vis x 1) stable index per spot
    for j = 1:vis
        if j <= numel(spot_stable) && spot_stable(j) > 0
            use_shift(j) = true;
        else
            warning('track_TS_spots:spotUnassigned', ...
                'Spot %d has no host nucleus - it is tracked WITHOUT correction.', j);
        end
    end
end

% --- preallocate ---
R_spot_n = zeros(zs, vis);
R_spot_n(1,:) = opts.Rp;
I_t       = zeros(zs, vis);
sd        = zeros(zs, vis);
spot_size = zeros(zs, vis);
C_cent    = zeros(zs+1, 2, vis);
for p = 1:vis
    C_cent(1,1,p) = C_in(p,1);
    C_cent(1,2,p) = C_in(p,2);
end

% --- fixed 16-bit display limits from the first frame ---
cmax = 65535;
if opts.show
    I_ref   = double(squeeze(MIP_image6d(1,opts.start_f,1,channel1,:,:)));
    v_ref   = sort(I_ref(:));
    disp_lo = v_ref(max(1, round(0.010*numel(v_ref))));
    disp_hi = v_ref(max(1, round(0.999*numel(v_ref))));
    if disp_hi <= disp_lo, disp_hi = disp_lo + 1; end
    col_marker = [cmax 0 0];
    col_text   = [cmax cmax cmax];
    if isempty(opts.fig)
        opts.fig = figure('Name','TS tracking','NumberTitle','off');
    else
        figure(opts.fig);
    end
end

%%===== cycle through all movie frames to track TSS=========
for j = 1:vis

    s = 0;
    if motion_correct && use_shift(j)
        s = nuc.spot_stable(j);
    end

    for k = 1:zs

        inputImageMS2_n = squeeze(MIP_image6d(1,k,1,channel1,:,:));

        % --- seed = previous position, optionally nudged by nucleus step ---
        seed_x = C_cent(k,1,j);
        seed_y = C_cent(k,2,j);
        if s > 0 && k >= 2
            ddy = nuc_shifts(k,1,s) - nuc_shifts(k-1,1,s);   % dy step frame k-1 -> k
            ddx = nuc_shifts(k,2,s) - nuc_shifts(k-1,2,s);   % dx step
            seed_x = seed_x + ddx;
            seed_y = seed_y + ddy;
        end

        out_E = track_spot2D(inputImageMS2_n, seed_x, seed_y, R_spot_n(k,j), opts.m_MS2);

        if out_E(3) >= 1                       % at least 1 bright spot found
            C_cent(k+1,1,j) = out_E(1);
            C_cent(k+1,2,j) = out_E(2);

            I1d_n = double(inputImageMS2_n);
            sd(k,j) = std(I1d_n(:));           % frame intensity SD

            R_spot_n(k+1,j) = opts.Rp;         % normal search radius next frame
        else                                   % out_E(3) == 0, nothing found
            C_cent(k+1,1,j) = seed_x;          % hold position (already drifted)
            C_cent(k+1,2,j) = seed_y;

            I1d_n = double(inputImageMS2_n);
            sd(k,j) = std(I1d_n(:));           % frame intensity SD
            
            R_spot_n(k+1,j) = opts.Rn;         % widen the search next frame
        end

        % --- visualization: fitted TSS on inverted 16-bit frame ---
        if opts.show
            I_disp = (double(inputImageMS2_n) - disp_lo) / (disp_hi - disp_lo);
            I_disp = min(max(I_disp, 0), 1);
            RGB    = repmat(uint16(cmax*(1 - I_disp)), [1 1 3]);   % dark spots on white

            position = [C_cent(k+1,1,j) C_cent(k+1,2,j) 6];
            RGB = insertObjectAnnotation(RGB, "circle", position, 'TSS', ...
                    'Color', col_marker, 'TextColor', col_text, 'LineWidth', 2);

            imshow(RGB, 'InitialMagnification', 'fit')
            if motion_correct && s > 0
                mc = 'MC on';
            else
                mc = 'MC off';
            end
            title(sprintf('spot %d/%d | frame %d/%d | %s', j, vis, k, zs, mc));
            drawnow
        end

        % --- TS spot intensity in a constant R_c mask ---
        I_t(k,j) = maskavg(inputImageMS2_n, C_cent(k+1,1,j), C_cent(k+1,2,j), opts.R_c);
    end
end

trk.C_cent          = C_cent;
trk.I_t             = I_t;
trk.sd              = sd;
trk.R_spot_n        = R_spot_n;
trk.spot_size       = spot_size;
trk.motion_correct  = motion_correct;
trk.corrected_spots = use_shift;
trk.opts            = opts;

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
