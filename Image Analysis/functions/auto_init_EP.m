function init = auto_init_EP(MIP, ch_enh, ch_prom, start_f, total_f, opt)
% AUTO_INIT_EP  Auto-initialisation for E-P tracking on a ~1-pair crop.
%
%   init = auto_init_EP(MIP_image6d, ch_enh, ch_prom, start_f, total_f, opt)
%
% Separates two concerns that used to be tangled:
%   * TRACKING ORDER  -> the brighter locus is the reference (tracked first,
%                        seeds the dimmer one). Purely for robustness.
%   * OUTPUT IDENTITY -> fixed by your labelling. ch_enh is ALWAYS the enhancer
%                        (SNAP-JF552 / RED); ch_prom is ALWAYS the promoter
%                        (Halo-JF646 / 642). The caller uses init.ref_is_enh to
%                        remap the tracked arrays back to these fields.
%
% INPUTS
%   MIP      : 6D MIP array [S T Z C Y X] from ReadImage6D2 (Z=1 for a MIP)
%   ch_enh   : enhancer channel  (JF552 / RED)   -- biological, fixed
%   ch_prom  : promoter channel  (JF646 / 642)   -- biological, fixed
%   start_f  : first frame to consider
%   total_f  : total frames available (metadata.SizeT)
%   opt      : optional struct (any omitted field uses the defaults below)
%
% OUTPUT struct 'init'
%   channel_ref, channel_sec : tracking order (ref = brighter, tracked first).
%                              Feed these into your loop's channel_enh/channel_prom.
%   m_ref, m_sec             : denoising m for the ref/sec tracking blocks.
%   C_in                     : [x y] seed for the reference (the first block).
%   ref_is_enh               : true if the brighter (reference) locus is the
%                              enhancer -> no remap needed at save time.
%   channel_enh, channel_prom: echoed biological channels (provenance).
%   seed_enh, seed_prom      : refined [x y] start positions per locus.
%   snr_enh, snr_prom        : start-frame SNR per locus.
%   start_f, stop_f          : usable frame range (stop = last frame before
%                              either locus vanishes).
%   snr_trace                : [T x 2] SNR per frame, columns [enhancer promoter].
%   notes                    : one-line summary.

if nargin < 6, opt = struct(); end

% ---------- defaults (override any via opt) ----------
def.searchR_seed   = 30;    % start-frame search radius if approx_center given
def.searchR_scan   = 12;    % radius to follow each spot frame-to-frame (px)
def.smooth_sigma   = 1;     % Gaussian sigma for peak finding (px)
def.snr_present    = 4;     % SNR below this = "absent" in that frame
def.stop_buffer    = 2;     % memory: tolerate up to this many consecutive absent
                            % frames; a longer gap ends the trajectory
def.m_min          = 0.20;  % m for a bright/clean spot
def.m_max          = 0.45;  % m for a dim/high-background spot
def.snr_hi         = 12;    % SNR mapped to m_min
def.snr_lo         = 3;     % SNR mapped to m_max
def.approx_center  = [];    % optional [x y] to disambiguate a busy crop
def.plot           = false; % show the SNR-vs-frame diagnostic
fn = fieldnames(def);
for i = 1:numel(fn)
    if ~isfield(opt, fn{i}) || isempty(opt.(fn{i})), opt.(fn{i}) = def.(fn{i}); end
end

T       = size(MIP, 2);
total_f = min(total_f, T);
frame   = @(k, c) squeeze(MIP(1, k, 1, c, :, :));
m_of    = @(snr) opt.m_max + (opt.m_min - opt.m_max) * ...
                 min(max((snr - opt.snr_lo) / (opt.snr_hi - opt.snr_lo), 0), 1);

% ---------- detect each locus in the start frame ----------
if ~isempty(opt.approx_center)
    c0 = opt.approx_center;  Rs = opt.searchR_seed;
else
    c0 = [NaN NaN];          Rs = [];     % global: brightest punctum per channel
end
sE = spot_snr_local(frame(start_f, ch_enh ), c0(1), c0(2), Rs, opt.smooth_sigma);
sP = spot_snr_local(frame(start_f, ch_prom), c0(1), c0(2), Rs, opt.smooth_sigma);

seed_enh  = [sE.x, sE.y];
seed_prom = [sP.x, sP.y];

% ---------- tracking order: brighter locus is the reference ----------
ref_is_enh = (sE.snr >= sP.snr);
if ref_is_enh
    channel_ref = ch_enh;   channel_sec = ch_prom;  sRef = sE;  C_in = seed_enh;
else
    channel_ref = ch_prom;  channel_sec = ch_enh;   sRef = sP;  C_in = seed_prom;
end
m_ref = m_of(sRef.snr);
if ref_is_enh, m_sec = m_of(sP.snr); else, m_sec = m_of(sE.snr); end

% ---------- follow both loci across ALL frames; record per-frame presence ----------
% Each frame we first search around the last good position; if that comes back
% weak we re-acquire globally (the crop holds ~1 punctum per channel). This stops
% a brief dip from freezing the guess and permanently killing the track.
biochan   = [ch_enh, ch_prom];
g         = [seed_enh; seed_prom];      % running guesses (row1=enh, row2=prom)
snr_trace = nan(total_f, 2);
present   = false(total_f, 2);
for k = start_f:total_f
    for c = 1:2
        I = frame(k, biochan(c));
        s = spot_snr_local(I, g(c,1), g(c,2), opt.searchR_scan, opt.smooth_sigma);
        if s.snr < opt.snr_present                       % re-acquire if it moved
            sg = spot_snr_local(I, NaN, NaN, [], opt.smooth_sigma);
            if sg.snr > s.snr, s = sg; end
        end
        snr_trace(k, c) = s.snr;
        if s.snr >= opt.snr_present
            g(c, :) = [s.x, s.y];
            present(k, c) = true;
        end
    end
end

% ---------- stop = last frame before either locus is gone past the buffer ----------
endK = [start_f start_f];               % per-channel last usable frame [enh prom]
for c = 1:2
    miss = 0;
    for k = start_f:total_f
        if present(k, c)
            endK(c) = k;  miss = 0;
        else
            miss = miss + 1;
            if miss > opt.stop_buffer    % gap longer than the buffer -> ends here
                break;                   % endK(c) holds the last present frame
            end
        end
    end
end
stop_f = max(min(endK), start_f);

if sE.snr < opt.snr_present || sP.snr < opt.snr_present
    warning('auto_init_EP: weak spot at start frame %d (SNR enh=%.1f, prom=%.1f).', ...
            start_f, sE.snr, sP.snr);
end

% ---------- pack + report ----------
init.channel_ref  = channel_ref;   init.channel_sec  = channel_sec;
init.m_ref        = m_ref;         init.m_sec        = m_sec;
init.C_in         = C_in;          init.ref_is_enh   = ref_is_enh;
init.channel_enh  = ch_enh;        init.channel_prom = ch_prom;
init.seed_enh     = seed_enh;      init.seed_prom    = seed_prom;
init.snr_enh      = sE.snr;        init.snr_prom     = sP.snr;
init.start_f      = start_f;       init.stop_f       = stop_f;
init.snr_trace    = snr_trace;
refName = 'enhancer'; if ~ref_is_enh, refName = 'promoter'; end
init.notes = sprintf(['ref = %s (ch%d, SNR %.1f, m %.2f) tracked first | ' ...
    'enhancer ch%d SNR %.1f | promoter ch%d SNR %.1f | seed [%d %d] | ' ...
    'frames %d-%d (enh ends %d, prom ends %d) | remap-needed: %d'], ...
    refName, channel_ref, sRef.snr, m_ref, ch_enh, sE.snr, ch_prom, sP.snr, ...
    C_in(1), C_in(2), start_f, stop_f, endK(1), endK(2), ~ref_is_enh);

if opt.plot
    figure('Name', 'auto_init_EP'); hold on;
    plot(start_f:total_f, snr_trace(start_f:total_f, 1), '-o', 'DisplayName', 'enhancer (JF552)');
    plot(start_f:total_f, snr_trace(start_f:total_f, 2), '-s', 'DisplayName', 'promoter (JF646)');
    yline(opt.snr_present, '--', 'SNR threshold');
    xline(stop_f, 'r-', 'stop\_f');
    xlabel('frame'); ylabel('SNR'); legend; title('Spot SNR vs frame'); hold off;
end
end
