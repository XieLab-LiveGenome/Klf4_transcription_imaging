% % % % MS2_TRACKING_PIPELINE.m

% % %  Single-channel SR-SIM timecourse:
% % %    Channel 1 = MS2        used for transcription state(bursting/inactive)
% % %
% % %  Pipeline:
% % %    Step 1 : Load CZI 
% % %    Step 2 : Nucleus segmentation and motion extraction
% % %    Step 3 : MS2 spot tracking starting from seed co-ordinates and using nucleus motion-correction
% % %    Step 4 : MS2 spot intensity extraction over time
% % %    Step 5 : HMM fitting of MS2 spot intensity for transcription state

% % %    Step 6 : Characterization of burst parameters from inferred HMM fitting
% % %    Step 7:  Save .mat
% % %
% % %  Dependencies on path:
% % %  ReadImage6D2, pick_TS_spots, nuc_motion_MS2, track_nuclei_MS2, track_TS_spots, burst_analysis_MS2

clear
clc
tic
%%============================== PROJECT DIRECTORY SETUP=============================
here = fileparts(mfilename('fullpath'));   % Folder of this main script
addpath(genpath(fullfile(here, 'functions')));   % Folder with the functions
addpath(genpath(fullfile(here, 'HMM fitting')));  % Folder with the HMM fitting package
addpath(fullfile(here, 'bioformats'));     % Folder with the bioformats package for image loading
javaaddpath(fullfile(here, 'bioformats', 'bioformats_package.jar'));    % Add the JAVA path for loading bioformats
%%==============================%%==============================%%===================

% ---Read CZI files of whole stack and MIP-----------------------------------------

save_filename='ctrl 1 11-1-2024 s1';

MIP_filename='/Volumes/xiel2lab/Aniket/ctrl 1_SIM_Maximum intensity projection_scaled.czi';

channel1=1;      % TSS CHANNEL, STAYGOLD

scene=1;    % set manually

MIP_out     = ReadImage6D2(MIP_filename, true, scene);
metadata    = MIP_out{2};
MIP_image6d = MIP_out{1};

% % % % % -----------Microscopy parameters ------------------------------------------
start_f = 1;        % start frame  (nucleus registration assumes start_f = 1)

xpixel = 0.0313;    % this is the real pixel value 20x
ypixel = 0.0313;

zs = metadata.SizeT;  % number of frames

time_int = 3;

% --- pipeline switches ---
run_nuclei     = true;     % true -> compute nucleus motion through cellpose segmentation
motion_correct = true;     % true -> seed TS search with host-nucleus motion correction

nuc_method = 'blocks';      % 'blocks' -> nuc_motion_MS2: segment ONCE on the
                            %             time average, then estimate motion by
                            %             registering block averages. 

nuc_block  = 10;            % frames per block
nuc_preproc = 'none';       % 'none' Raw MS2 intensities analyzed


%% ===== 1. PICK THE TS SPOTS =====
% Click each TS on a block average of the first frames. See pick_TS_spots.m.
% To reuse a previous set of positions, pass them with 'manual', C_in_manual.

C_in = pick_TS_spots(MIP_image6d, channel1, ...
        'start_f',  start_f, ...
        'n_frames', 10, ...
        'snap_R',   8, ...
        'invert',   true);

vis = size(C_in,1);   % no of TSS to be tracked


%% ===== 2. NUCLEUS SEGMENTATION, REGISTRATION & SPOT ASSIGNMENT =====
% Nucleus tracking uses track_nuclei_MS2

if run_nuclei
    switch lower(nuc_method)

        case 'blocks'   % registration of block averages (shape-independent)
            nuc = nuc_motion_MS2(MIP_image6d, channel1, C_in, ...
                    'block_size', nuc_block, ...
                    'pad',        100, ...
                    'sigma',      4, ...
                    'max_step',   150, ...
                    'preproc',    nuc_preproc, ...
                    'cell_D_seg', 160);

        case 'masks'    % Cellpose every frame, motion = mask centroid
            nuc = track_nuclei_MS2(MIP_image6d, channel1, C_in, ...
                    'block_size',   nuc_block, ...
                    'min_presence', 0.75, ...
                    'cell_D_seg',   160, ...
                    'max_drift_px', 200, ...
                    'alpha_ema',    0.3);

        otherwise
            error('nuc_method must be ''blocks'' or ''masks''.');
    end

    stable_ids      = nuc.stable_ids;      % nucleus IDs
    linked_masks    = nuc.linked_masks;    % 'blocks': one reference mask
                                           % 'masks' : (m x n x zs) per frame
    nuc_centroids   = nuc.centroids;       % (zs x 2 x n_stable) [y x] per frame
    nuc_shifts      = nuc.shifts;          % (zs x 2 x n_stable) [dy dx] from frame 1
    spot_nuc_gid    = nuc.spot_gid;        % nucleus ID per TS
    spot_nuc_stable = nuc.spot_stable;     % stable index per TS
else
    nuc = [];
    if motion_correct
        error('motion_correct = true requires run_nuclei = true.');
    end
end


% %% ===== 3. TS TRACKING =====
% TS spot tracking uses track_TS_spots.m. 

trk = track_TS_spots(MIP_image6d, channel1, C_in, nuc, motion_correct, ...
        'start_f', start_f, ...
        'Rp',      80, ...
        'Rn',      150, ...
        'm_MS2',   0.3, ...
        'R_c',     5, ...
        'show',    true);

C_cent    = trk.C_cent;      % (zs+1 x 2 x vis) [x y], row k+1 = frame k
I_t       = trk.I_t;         % raw intensity at the fitted spot
sd        = trk.sd;          % per-frame intensity SD
R_spot_n  = trk.R_spot_n;    % search radius used each frame
spot_size = trk.spot_size;   % fitted spot size (kept for compatibility)


% % % =================Intensity correction is histogram SD for HMM fitting=========
I_t_norm = zeros(zs,vis);
sd0 = sd(1,1);

for j=1:vis
for k=1:zs
    if sd(k,j) > 0 && sd0 > 0
I_t_norm(k,j)=I_t(k,j)*sd0/sd(k,j);
    else
I_t_norm(k,j)=I_t(k,j);
    end
end
end
% =================% % % =================% % % =================% % % =============


%% ===== 4. HMM FITTING, ON/OFF DURATIONS, BURST QUANTIFICATION =====
% All of it lives in burst_analysis_MS2.m. Outputs come back in one struct.

bur = burst_analysis_MS2(I_t, I_t_norm, time_int);

fitMS2_2s          = bur.fitMS2_2s;           % 2-state HMM fit
binary_2s          = bur.binary_2s;           % 1/0 ON/OFF states
on_times           = bur.on_times;            % {starts, ends, durations}
off_times          = bur.off_times;
on_cellavg         = bur.on_cellavg;          % per-trajectory mean ON duration
off_cellavg        = bur.off_cellavg;         % per-trajectory mean OFF duration
on_compiled        = bur.on_compiled;         % pooled ON durations
off_compiled       = bur.off_compiled;        % pooled OFF durations
Burst              = bur.Burst;               % per-burst size
Burst_amp          = bur.Burst_amp;           % per-burst amplitude
burst_cellavg      = bur.burst_cellavg;       % per-trajectory mean burst size
burst_amp_cellavg  = bur.burst_amp_cellavg;   % per-trajectory mean burst amplitude
Burst_compiled     = bur.Burst_compiled;      % pooled bursts
Burst_amp_compiled = bur.Burst_amp_compiled;  % pooled burst amplitudes


%%%=======Visualize intensity profile with HMM fitting ==============

length_t = zs;
time = 0:time_int:(length_t-1)*time_int;
figure;

for i = 1:vis
    subplot(vis,1,i);
    plot(time, I_t(:,i), 'Color', [0.2 0.9 0.0],'LineWidth', 3,'DisplayName', 'MS2 intensity');
    hold on;
    plot(time, fitMS2_2s(:,i), 'LineWidth', 3, 'DisplayName', 'HMM 2-state');

    ylabel(sprintf('TS Intensity %d', i));
    if i == 1
        title('MS2 intensity with HMM 2-state fitting');
    end
    if i == vis
        xlabel('Time(min)');
    end
    legend('show');
    grid on;
end


%%%======%%%%======%%%%======%%%%======

%%======SAVE key parameters==========

save(save_filename, 'C_in', 'C_cent' , 'sd' , 'spot_size' , 'I_t' , 'I_t_norm' , 'fitMS2_2s',"binary_2s", ...
    "on_times",'off_times',"on_cellavg","off_cellavg",'on_compiled','off_compiled', ...
    'Burst', 'Burst_amp',"Burst_compiled","Burst_amp_compiled",'burst_cellavg',"burst_amp_cellavg");

%%=====%%=====%%=====%%=====


toc
