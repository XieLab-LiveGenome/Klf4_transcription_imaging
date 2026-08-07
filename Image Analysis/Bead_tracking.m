%% =========Tracking multispectral beads imaged using SIM========

clear
clc
tic

% %============================== PROJECT DIRECTORY SETUP=============================
here = fileparts(mfilename('fullpath'));   % Folder of this main script
addpath(genpath(fullfile(here, 'functions')));   % Folder with the functions
addpath(fullfile(here, 'bioformats'));     % Folder with the bioformats package for image loading
javaaddpath(fullfile(here, 'bioformats', 'bioformats_package.jar'));    % Add the JAVA path for loading bioformats
% %==============================%%==============================%%===================

% % % %=========USER DEFINED PARAMETERS======================

start_f = 1;

save_filename = 'bead tracking Beads 642 1_SIM ROI 4 7-7-2026';
Input_zstack  = '/Volumes/Aniket2/Beads 642 7-7-2026/Beads 642 1_SIM ROI 4.czi';

ch_bead = 1;     % channel containing the beads (EDIT to your CZI channel order)
scene   = 1;     % set manually

R       = 20;    % 2D tracking radius
R_fit   = 8;    % 3D gaussian fitting radius
fit_thr = 0.2;   % minimum R2 for a good 3D gaussian fit
m_bead  = 0.1;  % SNR_inc2 contrast-stretch parameter (tune per dataset)

% % %=================================================

stack_out  = ReadImage6D2(Input_zstack, true, scene);
metadata   = stack_out{2};
full_stack = stack_out{1};   % 6D: [S T Z C Y X]

% ---- generate the MIP directly from the input stack (max over Z, dim 3) ----
MIP_image6d = max(full_stack, [], 3);   % -> [S T 1 C Y X]

xpixel = 0.0313;
ypixel = 0.0313;
zpixel = metadata.ScaleZ;
spacing = [xpixel, ypixel, zpixel];   % um per pixel

z_slice = metadata.SizeZ;   % number of z-slices
total_f = metadata.SizeT;   % total number of frames

stop_f = total_f;           % beads are fixed -> track the whole movie (lower if they bleach)
zs = min(stop_f, total_f);

% % =========PICK BEADS MANUALLY FROM FIRST-FRAME MIP================
detImg = squeeze(MIP_image6d(1, start_f, 1, ch_bead, :, :));
figure; imshow((65536-detImg), [], 'InitialMagnification', 'fit');
title('Click on each bead, then press Enter');
[cx, cy] = ginput;                         % click beads, Enter to finish
close(gcf);

C_in = [cx cy];                            % [x y] = [col row]
vis  = size(C_in, 1);
if vis == 0
    error('No beads selected.');
end
fprintf('Selected %d bead(s) to track.\n', vis);


vis  = size(C_in, 1);

% % ========Declare variables===========

bead_xyz = zeros(zs,3,vis);
C_cent   = zeros(zs+1,3,vis);   % row start_f holds the detected seed
r2       = zeros(zs,1,vis);

for k = 1:vis
    C_cent(start_f,1,k) = C_in(k,1);
    C_cent(start_f,2,k) = C_in(k,2);
    C_cent(start_f,3,k) = round(z_slice/2);
end

% % % %%==================BEAD TRACKING=================%=================
% % % %%==================BEAD TRACKING=================%=================
for j = 1:vis
for k = start_f:zs

    inputImage  = MIP_image6d(1,k,1,ch_bead,:,:);
    inputImage  = squeeze(inputImage);
    inputImage2 = SNR_inc2(inputImage, m_bead);    % denoising / contrast stretch

    stack_B = single_zstack(full_stack,1,k,ch_bead);   % z-stack for this frame (once)

    out_B = track_spot2D(inputImage2, C_cent(k,1,j), C_cent(k,2,j), R, m_bead);

    if out_B(3) >= 1     % at least 1 bead found in MIP
        C_cent(k+1,1,j) = out_B(1);
        C_cent(k+1,2,j) = out_B(2);    % precise 2D position

        guess = [C_cent(k+1,1,j), C_cent(k+1,2,j), C_cent(k,3,j)];
        [x_ref, y_ref, z_ref, params, R2] = fit_Gaussian3D(stack_B, guess, R_fit, spacing);

        if ~isnan(x_ref)   % 3D fit succeeded
            bead_xyz(k,1,j) = x_ref*xpixel;
            bead_xyz(k,2,j) = y_ref*ypixel;
            bead_xyz(k,3,j) = (z_ref-1)*zpixel;
            C_cent(k+1,3,j) = z_ref;
            r2(k,1,j)       = R2;
        else               % fit crashed -> keep 2D, hold z
            bead_xyz(k,1,j) = C_cent(k+1,1,j)*xpixel;
            bead_xyz(k,2,j) = C_cent(k+1,2,j)*ypixel;
            bead_xyz(k,3,j) = (C_cent(k,3,j)-1)*zpixel;
            C_cent(k+1,3,j) = C_cent(k,3,j);
        end
    end

    % --- fallback: 2D tracking failed or bad fit -> wide fit at last position ---
    if k > start_f
        if out_B(3) == 0 || r2(k,1,j) <= fit_thr
            guess = [bead_xyz(k-1,1,j)/xpixel, ...
                     bead_xyz(k-1,2,j)/ypixel, ...
                     bead_xyz(k-1,3,j)/zpixel + 1];
            [x_ref, y_ref, z_ref, params, R2] = fit_Gaussian3D(stack_B, guess, R, spacing);

            bead_xyz(k,1,j) = x_ref*xpixel;
            bead_xyz(k,2,j) = y_ref*ypixel;
            bead_xyz(k,3,j) = (z_ref-1)*zpixel;
            C_cent(k+1,1,j) = x_ref;
            C_cent(k+1,2,j) = y_ref;
            C_cent(k+1,3,j) = z_ref;
            r2(k,1,j)       = R2;
        end
    end

    k

    % %%%%%=======================visualization=======
    label = 'bead';
    position = [C_cent(k+1,1,j) C_cent(k+1,2,j) 5];
    if ~isnan(position(1,1))
        I16bit = uint16(65536 - inputImage2);
        RGB = insertObjectAnnotation(I16bit,"circle",position,label);
        imshow(RGB,[])
    end
    % %%=============

end
end

thr_z = 0.3; %%(300 nm z-movement in ~5 sec is most likely due to bad fitting in the frame)

for j=1:vis
    for k=2:zs-1
        if abs(bead_xyz(k,3,j)-bead_xyz(k-1,3,j)) > thr_z || bead_xyz(k,3,j) < 0
            bead_xyz(k,3,j) = (bead_xyz(k-1,3,j) + bead_xyz(k+1,3,j))/2;
        end
    end
end

% % % =====Localization precision / drift summary (nm)============
bead_stats = struct();
for j = 1:vis
    xyz = bead_xyz(:,:,j);                                      % zs x 3, um
    bead_stats(j).precision_nm = std(xyz,0,1)*1000;            % [sx sy sz]
    bead_stats(j).net_drift_nm = (xyz(end,:)-xyz(1,:))*1000;   % [dx dy dz]
    bead_stats(j).range_nm     = (max(xyz,[],1)-min(xyz,[],1))*1000;
    fprintf(['Bead %d | precision (x,y,z) = %.1f, %.1f, %.1f nm | ' ...
             'net drift = %.1f, %.1f, %.1f nm\n'], j, ...
             bead_stats(j).precision_nm, bead_stats(j).net_drift_nm);
end
% % % ==============================================================

save(save_filename, 'bead_xyz', 'C_cent', 'r2', 'bead_stats', 'C_in', 'spacing');
toc