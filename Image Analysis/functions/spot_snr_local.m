function s = spot_snr_local(I, cx, cy, searchR, sig)
% SPOT_SNR_LOCAL  Robust single-spot detection + SNR in one 2D frame.
%
%   s = spot_snr_local(I, cx, cy, searchR, sig)
%
%   I        : 2D image (any numeric type)
%   cx, cy   : approximate centre to search around (x = column, y = row).
%              Pass NaN for both to do a GLOBAL search (brightest spot in frame).
%   searchR  : peak-search radius in px (only used when cx,cy are given).
%   sig      : Gaussian smoothing sigma (px) for peak finding (default 1).
%
%   Returns struct s with fields:
%       x, y   : refined peak location (px, x=col / y=row)
%       peak   : mean intensity over the spot core (raw image)
%       bg     : robust background median (raw image)
%       noise  : robust background sigma (1.4826*MAD, no toolbox needed)
%       snr    : (peak - bg) / noise
%       fmax   : max of the smoothed frame (hot-pixel-safe)
%
% Tunables (kept as local constants so they live next to the maths):
%   coreR  - radius of the disk used to estimate spot signal
%   exclR  - pixels within this radius of the peak are excluded from background
coreR = 2;
exclR = 8;          % should comfortably exceed the SIM PSF/array-spot extent

I = double(I);
[ny, nx] = size(I);
if nargin < 5 || isempty(sig), sig = 1; end

Is = imgaussfilt(I, sig);                 % light smoothing for robust peak finding
[X, Y] = meshgrid(1:nx, 1:ny);

% --- restrict the peak search if a centre was supplied ---
if nargin >= 4 && ~isempty(searchR) && all(~isnan([cx cy]))
    win = (X - cx).^2 + (Y - cy).^2 <= searchR^2;
else
    win = true(ny, nx);                   % global search
end
Iw = Is; Iw(~win) = -inf;
[~, idx]  = max(Iw(:));
[py, px]  = ind2sub([ny, nx], idx);       % py = row (y), px = col (x)

% --- signal: mean over the spot core on the RAW image ---
core = (X - px).^2 + (Y - py).^2 <= coreR^2;
peak = mean(I(core));

% --- background: robust stats outside an exclusion disk ---
bgMask = (X - px).^2 + (Y - py).^2 > exclR^2;
bv     = I(bgMask);
bg     = median(bv);
noise  = 1.4826 * median(abs(bv - bg));   % robust sigma, no Stats toolbox needed
if noise <= 0, noise = std(bv) + eps; end

s.x = px;  s.y = py;
s.peak = peak;  s.bg = bg;  s.noise = noise;
s.snr  = (peak - bg) / noise;
s.fmax = max(Is(:));
end
