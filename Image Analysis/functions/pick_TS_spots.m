function [C_in, I_pick] = pick_TS_spots(MIP_image6d, channel1, varargin)
%PICK_TS_SPOTS  Click transcription sites on a block average of the first frames.
%
%   [C_in, I_pick] = PICK_TS_SPOTS(MIP_image6d, channel1)
%   [C_in, I_pick] = PICK_TS_SPOTS(..., 'Name', value, ...)
%   [C_in, I_pick] = PICK_TS_SPOTS(..., optsStruct)
%
%   Averages the first frames of the TSS channel, shows the image, and lets
%   you click every TS. Zoom / pan first, then press ENTER in the command
%   window to start picking. BACKSPACE undoes the last point, ENTER finishes.
%
%   Options (defaults in brackets):
%       start_f   [1]     first frame of the average
%       n_frames  [10]    number of frames averaged
%       snap_R    [8]     snap each click to the brightest pixel within this
%                         radius, in px (0 = keep the raw click)
%       invert    [true]  display inverted (dark spots on white)
%       manual    []      if non-empty, skip picking and return it as C_in
%
%   Outputs:
%       C_in   - (n x 2) [x y] positions, same convention as the tracking code
%       I_pick - the averaged image used for picking

opts = struct('start_f',1, 'n_frames',10, 'snap_R',8, 'invert',true, 'manual',[]);
opts = local_parse_opts(opts, varargin);

% --- build the block average ---
T  = size(MIP_image6d, 2);
f1 = max(1, opts.start_f);
f2 = min(f1 + opts.n_frames - 1, T);

I_pick = zeros(size(squeeze(MIP_image6d(1,1,1,channel1,:,:))));
for f = f1:f2
    I_pick = I_pick + double(squeeze(MIP_image6d(1,f,1,channel1,:,:)));
end
I_pick = I_pick / (f2 - f1 + 1);

[m_pick, n_pick] = size(I_pick);

% --- bypass picking if positions were supplied ---
if ~isempty(opts.manual)
    C_in = opts.manual;
    return
end

% --- robust display limits (no Statistics Toolbox needed) ---
v_pick = sort(I_pick(:));
pk_lo  = v_pick(max(1, round(0.010*numel(v_pick))));
pk_hi  = v_pick(max(1, round(0.999*numel(v_pick))));
if pk_hi <= pk_lo, pk_hi = pk_lo + 1; end

hFig = figure('Name','Pick TS spots','NumberTitle','off');
if opts.invert
    imshow(pk_hi - I_pick, [0 pk_hi - pk_lo], 'InitialMagnification', 'fit');
else
    imshow(I_pick, [pk_lo pk_hi], 'InitialMagnification', 'fit');
end
title(sprintf('Average of frames %d-%d  -  click each TS, ENTER when done', f1, f2));

fprintf('\n=== PICK TS SPOTS ===\n');
fprintf('  Averaged frames %d-%d of channel %d\n', f1, f2, channel1);
fprintf('  1. zoom / pan in the figure if needed\n');
fprintf('  2. press ENTER here to start picking\n');
fprintf('  3. click each transcription site (BACKSPACE undoes the last point)\n');
fprintf('  4. press ENTER in the figure when finished\n\n');
input('Press ENTER when ready to start picking...','s');

if exist('getpts','file')
    [xi, yi] = getpts(hFig);
else
    [xi, yi] = ginput;          % fallback if getpts is unavailable
end

if isempty(xi)
    error('pick_TS_spots:noSelection', 'No TS spots were selected.');
end

% --- drop accidental duplicates (a double-click can add the same point twice) ---
P    = [xi(:) yi(:)];
keep = true(size(P,1),1);
for a = 2:size(P,1)
    if any(sqrt(sum((P(1:a-1,:) - P(a,:)).^2, 2)) < 2)
        keep(a) = false;
    end
end
C_in = P(keep,:);               % column 1 = x, column 2 = y

% --- snap each click to the local intensity maximum (on the RAW average) ---
if opts.snap_R > 0
    for p = 1:size(C_in,1)
        cx = round(C_in(p,1));
        cy = round(C_in(p,2));
        x1 = max(1, cx - opts.snap_R);  x2 = min(n_pick, cx + opts.snap_R);
        y1 = max(1, cy - opts.snap_R);  y2 = min(m_pick, cy + opts.snap_R);
        sub = I_pick(y1:y2, x1:x2);
        [~, idx] = max(sub(:));
        [ry, rx] = ind2sub(size(sub), idx);
        C_in(p,1) = x1 + rx - 1;
        C_in(p,2) = y1 + ry - 1;
    end
end

% --- show what was picked ---
figure(hFig); hold on
plot(C_in(:,1), C_in(:,2), 'ro', 'MarkerSize', 12, 'LineWidth', 1.5);
for p = 1:size(C_in,1)
    text(C_in(p,1)+12, C_in(p,2), sprintf('%d', p), ...
        'Color', 'y', 'FontSize', 12, 'FontWeight', 'bold');
end
hold off
title(sprintf('%d TS spots selected', size(C_in,1)));

fprintf('\n%d TS spots selected. To reuse them later, paste:\n\n', size(C_in,1));
fprintf('C_in_manual = [ ...\n');
fprintf('    %8.1f %8.1f\n', C_in');
fprintf('    ];\n\n');

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
