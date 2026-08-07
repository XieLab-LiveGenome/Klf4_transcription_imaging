% % % % Enhancer-Promoter-MS2-BRD4 hub_TRACKING_PIPELINE.m

% % %  4-channel SR-SIM timecourse:
% % %    Channel 1 = MS2        used for transcription state(bursting/inactive)
% % %    Channel 2 = Enhancer   used for enhancer 3D location
% % %    Channel 3 = Promoter   used for enhancer 3D location
% % %    Channel 4 = Protein hubs     used for searching hubs in the enhancer proximity
% % %
% % %  Pipeline:
% % %    Step 1 : Load CZI (all channels) 
% % %    Step 2 : Enhancer 3D tracking starting from seed co-ordinates
% % %    Step 3 : Promoter 3D tracking starting from enhancer co-ordinates
% % %    Step 4 : MS2 spot tracking starting from enhancer co-ordinates
% % %    Step 5 : HMM fitting of MS2 spot intensity for transcription state
% % %    Step 6 : Proximal condensate/Hub peak detection in enhancer local neighborhood
% % %    Step 7 : Enhancer-Promoter/Enhancer-Hub/Promoter-Hub 3D distance calculations
% % %    Step 8:  Save .mat
% % %
% % %  Dependencies on path:
% % %  ReadImage6D2, SNR_inc2, track_spot2D, single_zstack, fit_Gaussian3D, maskavg, HMM_fit_fun, binary, condensate_search_v9

clear
clc
tic

%============================== PROJECT DIRECTORY SETUP=============================
here = fileparts(mfilename('fullpath'));   % Folder of this main script
addpath(genpath(fullfile(here, 'functions')));   % Folder with the functions
addpath(genpath(fullfile(here, 'HMM fitting')));  % Folder with the HMM fitting package
addpath(fullfile(here, 'bioformats'));     % Folder with the bioformats package for image loading
javaaddpath(fullfile(here, 'bioformats', 'bioformats_package.jar'));    % Add the JAVA path for loading bioformats
%%==============================%%==============================%%===================

% % %=========USER DEFINED PARAMETERS======================
C_in =[262 225];   % Approximate pixel co-ordinates of enhancer in start frame/serves as trajectory seed
% 
start_f = 1;       % START frame
stop_f = 91;       % STOP frame
% % 
save_filename='4 color 4c ctrl 3-24-2026 9_SIM 4000 threshold.mat';
Input_zstack  = '/Volumes/Aniket2/4C g67 3-24-2026/4c ctrl 3-24-2026 9_SIM.czi';
MIP_filename='/Volumes/Aniket2/4C g67 3-24-2026/4c ctrl 3-24-2026 9_SIM_Maximum intensity projection.czi';

channel_MS2=1;      % 405 channel
channel_enh=2;      % 488 channel
channel_prom=3;      % 561 channel
channel_cond=4;      % 642 channel

scene=1;    % set manually

R=40;  %tracking radius for enhancer
R_p = 20; %search radius for promoter
R_m = 20; %search radius for MS2

R_fit = 10;  % 3d gaussian fitting radius
fit_thr = 0.2; % minimum R2 value needed for a good gaussian fit from tracked spot
% % 
m_enh = 0.2;
m_prom = 0.2;
% 

% % % 
% % % %%%=================================================
% % 
stack_out = ReadImage6D2(Input_zstack, true, scene);
metadata = stack_out{2};
full_stack = stack_out{1};

MIP_out = ReadImage6D2(MIP_filename, true, scene);
MIP_image6d = MIP_out{1};

xpixel = 0.0313;
ypixel = 0.0313;
zpixel = metadata.ScaleZ;

spacing = [xpixel, ypixel, zpixel];  % μm per pixel (example from CZI metadata)
z_slice = metadata.SizeZ; % number of z-slices
total_f = metadata.SizeT;  % total number of frames

% % ====================================

zs = min(stop_f,total_f);    % resizing the matrix sizes according to the 

vis=size(C_in,1);  % no of EC pairs to be tracked

% % ========Declare variables===========

% % % JF 646 prom is clear so tracked first then enhancer with JF552

enh_xyz=zeros(zs,3,vis);   
prom_xyz=zeros(zs,3,vis);

C_cent_e=zeros(zs+1,3,vis);   % the 1st position is the user defined initial condition
C_cent_p=zeros(zs+1,3,vis);   % the 1st position is the user defined initial condition
C_cent_TSS = zeros(zs,3,vis);


TSS_xyz=zeros(zs,3,vis);      % TS xyz co-ordinates
MS2_score=zeros(zs,5,vis);    % MS2 intensity information


% =====primary-- this is closest to enhancer ======
Cond_xyz1=zeros(zs,3,vis);  % xyz co-ordinates 
Cond_s1=zeros(zs,vis);   %  SIZE OF CONDENSATE in MIP
Cond_int1=zeros(zs,vis);  % BRIGHTNESS OF CONDENSATE in MIP
Cond_AR1=zeros(zs,vis);   % AR OF CONDENSATE in MIP

% =====secondary -- this is closest to promoter ======
Cond_xyz2=zeros(zs,3,vis);  % xyz co-ordinates 
Cond_s2=zeros(zs,vis);   %  SIZE OF CONDENSATE in MIP
Cond_int2=zeros(zs,vis);  % BRIGHTNESS OF CONDENSATE in MIP
Cond_AR2=zeros(zs,vis);   % AR OF CONDENSATE in MIP

r2=zeros(zs,2,vis);

for k=1:vis
    C_cent_e(start_f,1,k)=C_in(k,1);
    C_cent_e(start_f,2,k)=C_in(k,2);
    C_cent_e(start_f,3,k)= round(z_slice/2);
    C_cent_p(start_f,1,k)=C_in(k,1);
    C_cent_p(start_f,2,k)=C_in(k,2);
    C_cent_p(start_f,3,k)= round(z_slice/2);
end

% % % % % % %%==================ENHANCER%=================%=================%=================
% % % % % % %%==================ENHANCER%=================%=================%=================
% % % 
for j=1:vis

for k=start_f:zs

inputImage_enh = MIP_image6d(1,k,1,channel_enh,:,:);
inputImage_enh = squeeze(inputImage_enh);

inputImage_enh2 = SNR_inc2(inputImage_enh,m_enh);  % denoising through contrast stretching 

out_E = track_spot2D(inputImage_enh2, C_cent_e(k,1,j), C_cent_e(k,2,j),R, m_enh);      % track and 2d gaussian fitting based on guess/previous frame position

% === SKIP FRAME if no puncta detected in MIP ===
if out_E(3) == 0
    enh_xyz(k,1,j) = NaN;
    enh_xyz(k,2,j) = NaN;
    enh_xyz(k,3,j) = NaN;
    C_cent_e(k+1,1,j) = C_cent_e(k,1,j);   % carry forward previous guess
    C_cent_e(k+1,2,j) = C_cent_e(k,2,j);
    C_cent_e(k+1,3,j) = C_cent_e(k,3,j);
    r2(k,1,j) = NaN;
    fprintf('Enhancer | cell %d | frame %d/%d  -- SKIPPED (no puncta)\n', j, k, zs);
    continue
end

% --- At least 1 enhancer spot FOUND in MIP ---
C_cent_e(k+1,1,j) = out_E(1);
C_cent_e(k+1,2,j) = out_E(2);      % precise 2D position of enhancer

stack_e = single_zstack(full_stack,1,k,channel_enh);

guess_e = [C_cent_e(k+1,1,j), C_cent_e(k+1,2,j), C_cent_e(k,3,j)];     % guess latest 2D co-ordinates and 3D co-ordinate from last frame

[x_ref_e, y_ref_e, z_ref_e, params_e,R2_e] = fit_Gaussian3D(stack_e, guess_e, R_fit, spacing);

% Check if 3D gaussian fit was successful or crashes
if ~isnan(x_ref_e) 

enh_xyz(k,1,j) = x_ref_e*xpixel;
enh_xyz(k,2,j) = y_ref_e*ypixel;
enh_xyz(k,3,j) = (z_ref_e-1)*zpixel;
C_cent_e(k+1,3,j) = z_ref_e;
r2(k,1,j) = R2_e;

else

enh_xyz(k,1,j) = C_cent_e(k+1,1,j)*xpixel;
enh_xyz(k,2,j) = C_cent_e(k+1,2,j)*ypixel;
enh_xyz(k,3,j) = (C_cent_e(k,3,j)-1)*zpixel;
C_cent_e(k+1,3,j) = C_cent_e(k,3,j);

end

% --- Fallback: spot found but bad 3D fit (R2 below threshold) ---
if k>start_f && r2(k,1,j) <= fit_thr

guess_e = [enh_xyz(k-1,1,j)/xpixel, enh_xyz(k-1,2,j)/ypixel, enh_xyz(k-1,3,j)/zpixel + 1]; 

[x_ref_e, y_ref_e, z_ref_e, params,R2_e] = fit_Gaussian3D(stack_e, guess_e, R, spacing);     % much wider area for fitting based on last position

enh_xyz(k,1,j) = x_ref_e*xpixel;
enh_xyz(k,2,j) = y_ref_e*ypixel;
enh_xyz(k,3,j) = (z_ref_e-1)*zpixel;

C_cent_e(k+1,1,j) = x_ref_e;
C_cent_e(k+1,2,j) = y_ref_e;
C_cent_e(k+1,3,j) = z_ref_e;

end

fprintf('Enhancer | cell %d | frame %d/%d\n', j, k, zs);

% %%%%%=======================visualization=======
label='enh';
position = [enh_xyz(k,1,j)/xpixel enh_xyz(k,2,j)/ypixel 5];

if ~isnan(position(1,1))

I16bit=uint16(65536-inputImage_enh2);

RGB = insertObjectAnnotation(I16bit,"circle",position,label);

imshow(RGB,[],'InitialMagnification', 800)

end
%%=============%%=============%%=============%%====

end
end

% %%==================PROMOTER%=================%=================%=================
% %%==================PROMOTER%=================%=================%=================
for j=1:vis

for k= start_f:zs

inputImage_prom = MIP_image6d(1,k,1,channel_prom,:,:);
inputImage_prom = squeeze(inputImage_prom);

inputImage_prom2 = SNR_inc2(inputImage_prom,m_prom); % denoising 

out_P = track_spot2D(inputImage_prom2, C_cent_e(k+1,1,j), C_cent_e(k+1,2,j),R_p, m_prom);     % track and 2d gaussian fitting based on guess/previous frame position

% === SKIP FRAME if no puncta detected in MIP ===
if out_P(3) == 0
    prom_xyz(k,1,j) = NaN;
    prom_xyz(k,2,j) = NaN;
    prom_xyz(k,3,j) = NaN;
    C_cent_p(k+1,1,j) = C_cent_p(k,1,j);   % carry forward previous guess
    C_cent_p(k+1,2,j) = C_cent_p(k,2,j);
    C_cent_p(k+1,3,j) = C_cent_p(k,3,j);
    r2(k,2,j) = NaN;
    fprintf('Promoter | cell %d | frame %d/%d  -- SKIPPED (no puncta)\n', j, k, zs);
    continue
end

% --- At least 1 promoter spot FOUND in MIP ---
C_cent_p(k+1,1,j) = out_P(1);
C_cent_p(k+1,2,j) = out_P(2);      % precise 2D position of promoter 

stack_p = single_zstack(full_stack,1,k,channel_prom);

guess_p = [C_cent_p(k+1,1,j), C_cent_p(k+1,2,j), C_cent_e(k+1,3,j)];  % guess latest 2D co-ordinates and 3D co-ordinate from last frame

[x_ref_p, y_ref_p, z_ref_p, params, R2_p] = fit_Gaussian3D(stack_p, guess_p, R_fit, spacing);

% Check if 3D gaussian fit was successful or crashes
if ~isnan(x_ref_p) 

prom_xyz(k,1,j) = x_ref_p*xpixel;
prom_xyz(k,2,j) = y_ref_p*ypixel;
prom_xyz(k,3,j) = (z_ref_p-1)*zpixel;
C_cent_p(k+1,3,j) = z_ref_p;
r2(k,2,j) = R2_p;

else

prom_xyz(k,1,j) = C_cent_p(k+1,1,j)*xpixel;
prom_xyz(k,2,j) = C_cent_p(k+1,2,j)*ypixel;
prom_xyz(k,3,j) = (C_cent_p(k,3,j)-1)*zpixel;
C_cent_p(k+1,3,j) = C_cent_p(k,3,j);

end

% --- Fallback: spot found but bad 3D fit (R2 below threshold) ---
if r2(k,2,j) <= fit_thr

guess_p = [enh_xyz(k,1,j)/xpixel, enh_xyz(k,2,j)/ypixel, enh_xyz(k,3,j)/zpixel + 1]; 

[x_ref_p, y_ref_p, z_ref_p, params,R2_p] = fit_Gaussian3D(stack_p, guess_p, R, spacing);  % much bigger radius to be fitted

prom_xyz(k,1,j) = x_ref_p*xpixel;
prom_xyz(k,2,j) = y_ref_p*ypixel;
prom_xyz(k,3,j) = (z_ref_p-1)*zpixel;

C_cent_p(k+1,1,j) = x_ref_p;
C_cent_p(k+1,2,j) = y_ref_p;
C_cent_p(k+1,3,j) = z_ref_p;

r2(k,2,j) = R2_p;
end

fprintf('Promoter | cell %d | frame %d/%d\n', j, k, zs);

%%%%%=======================visualization=======
label='prom';
position = [C_cent_p(k+1,1,j) C_cent_p(k+1,2,j) 5];

if ~isnan(position(1,1))

I16bit=uint16(65536-inputImage_prom2);

RGB = insertObjectAnnotation(I16bit,"circle",position,label);

imshow(RGB,[],'InitialMagnification', 800)
end
%%=============%%=============%%================

end
end


% % % === Interpolate NaN frames (skipped puncta) then filter odd z-loc ===

thr_z = 0.6; %%(600 nm z-movement in ~60 sec is most likely due to bad fitting in the frame)

% --- Step 1: Linear interpolation to fill NaN frames from skipped puncta ---
for j=1:vis
    for dim=1:3
        enh_col = enh_xyz(start_f:zs, dim, j);
        prom_col = prom_xyz(start_f:zs, dim, j);

        if any(isnan(enh_col)) && ~all(isnan(enh_col))
            enh_xyz(start_f:zs, dim, j) = fillmissing(enh_col, 'linear');
        end
        if any(isnan(prom_col)) && ~all(isnan(prom_col))
            prom_xyz(start_f:zs, dim, j) = fillmissing(prom_col, 'linear');
        end
    end

    n_nan_e = sum(any(isnan(enh_xyz(start_f:zs,:,j)),2));
    n_nan_p = sum(any(isnan(prom_xyz(start_f:zs,:,j)),2));
    if n_nan_e > 0
        fprintf('WARNING: cell %d enhancer has %d frames that could not be interpolated (all NaN neighbors)\n', j, n_nan_e);
    end
    if n_nan_p > 0
        fprintf('WARNING: cell %d promoter has %d frames that could not be interpolated (all NaN neighbors)\n', j, n_nan_p);
    end
end

% --- Step 2: Filter out aberrant z-jumps (>600 nm between consecutive frames) ---
for j=1:vis
    for k=2:zs-1
        if abs(prom_xyz(k,3,j)-prom_xyz(k-1,3,j)) > thr_z || prom_xyz(k,3,j) < 0
            prom_xyz(k,3,j) = (prom_xyz(k-1,3,j) + prom_xyz(k+1,3,j))/2;
        end

        if abs(enh_xyz(k,3,j)-enh_xyz(k-1,3,j)) > thr_z || enh_xyz(k,3,j) < 0
            enh_xyz(k,3,j) = (enh_xyz(k-1,3,j) + enh_xyz(k+1,3,j))/2;
        end
    end
end
% % ==============%%=============%%=============


m_MS2 = 0.1;
% ================MS2 check=============================================
% ================%================%================%================

for j=1:vis
for k= start_f:zs

stack_M = single_zstack(full_stack,1,k,channel_MS2);  % Channel, Timepoint

inputImage_MS2 = MIP_image6d(1,k,1,channel_MS2,:,:);
inputImage_MS2 = squeeze(inputImage_MS2);

out_M = track_spot2D(inputImage_MS2, C_cent_e(k+1,1,j), C_cent_e(k+1,2,j),R_m,m_MS2);

guess_m = [out_M(1), out_M(2), C_cent_e(k+1,3,j)];

[x_ref_m, y_ref_m, z_ref_m, params_m] = fit_Gaussian3D(stack_M, guess_m, R_fit, spacing);

TSS_xyz(k,1,j) = x_ref_m*xpixel;
TSS_xyz(k,2,j) = y_ref_m*ypixel;
TSS_xyz(k,3,j) = (z_ref_m-1)*zpixel;

C_cent_TSS(k,1,j) = round(x_ref_m);
C_cent_TSS(k,2,j) = round(y_ref_m);
C_cent_TSS(k,3,j) = round(z_ref_m);

MS2_score(k,1,j) = 1;
MS2_score(k,2,j) = maskavg(inputImage_MS2,C_cent_TSS(k,1,j),C_cent_TSS(k,2,j),5);

if out_M(3)==0
TSS_xyz(k,1,j) = NaN;
TSS_xyz(k,2,j) = NaN;
TSS_xyz(k,3,j) = NaN;
MS2_score(k,1,j)= -1;
MS2_score(k,2,j) = maskavg(inputImage_MS2,prom_xyz(k,1,j)/xpixel,prom_xyz(k,2,j)/ypixel,5);
end

end
end

% % ================Hidden Markov model (2-state fitting) =================
for j=1:vis
inter = HMM_fit_fun(MS2_score(start_f:zs, 2, j));
MS2_score(start_f:zs,4,j) = inter{1,1};  % 2 state binary model of each trajectory
end

for j=1:vis
MS2_score(start_f:zs,3,j) = binary (MS2_score(start_f:zs,4,j));  % 1 or 0 value assigned to ON/OFF states model of each trajectory
end

for j=1:vis
I_max = max(MS2_score(:,2,j));
MS2_score(:,5,j) = MS2_score(:,2,j)/I_max;
end
% 
% % % % ================% % % ================% % % ================% % % ================
% 
% % % % ==================CONDENSATE%=================%=================%=================
% % % % ==================CONDENSATE%=================%=================%=================
% 
figure 

nuc_ave =4000; % average nuclear background
mult = 3.0;   % hub detection threshold = mult*average nuclear background

for j=1:vis
for k= start_f:zs

inputImage_cond = MIP_image6d(1,k,1,channel_cond,:,:);
inputImage_cond = squeeze(inputImage_cond);


%%=======%%=======%%=======%%=======%%=======%%=======%%=======%%=======

%%%%= for timelapse imaging===== 

inputImage_cond2 = SNR_inc2(inputImage_cond,0.1);     % BACKGROUND REDUCTION AND SIGNAL BOOSTING

sd(k,j)=std(inputImage_cond(:));

t_cond=nuc_ave*sd(k,j)/sd(start_f,j);   % adpative thresolding based on image brightness

BW = imbinarize(inputImage_cond,t_cond*mult);

BW2 = bwareafilt(BW,[10 inf]);   % filtering very small condensates 

BW3 = inputImage_cond().*BW2;

%%%======================

% -------------% Extract condensate stack and get enh/prom position to search -------------% -------------% ---

stack_C = single_zstack(full_stack,1,k,channel_cond);

R_search_c  = 20; %1st search radius is 20 pixel or ~ 600 nm 

enh_pos = [enh_xyz(k,1,j)  enh_xyz(k,2,j)  enh_xyz(k,3,j)];
prom_pos = [prom_xyz(k,1,j)  prom_xyz(k,2,j)  prom_xyz(k,3,j)];

% -------------% -------------% -------------% -------------% -------------% -------------% -------------

%%%%======== condensate closest to enhancer ========

OUT_e = condensate_search_v9(stack_C, BW3, enh_pos, spacing, t_cond, mult, R_search_c);

Cond_xyz1(k,1,j) = OUT_e{1}(1,1);
Cond_xyz1(k,2,j) = OUT_e{1}(1,2);
Cond_xyz1(k,3,j) = OUT_e{1}(1,3);
Cond_s1(k,j) = OUT_e{1}(1,4);
Cond_int1(k,j) = OUT_e{1}(1,5);
Cond_AR1(k,j) = OUT_e{1}(1,6);

%%%%======== condensate closest to promoter ========

if ~isnan(prom_pos(1)) 

OUT_p = condensate_search_v9(stack_C, BW3, prom_pos, spacing, t_cond, mult, R_search_c);

Cond_xyz2(k,1,j) = OUT_p{1}(1,1);
Cond_xyz2(k,2,j) = OUT_p{1}(1,2);
Cond_xyz2(k,3,j) = OUT_p{1}(1,3);
Cond_s2(k,j) = OUT_p{1}(1,4);
Cond_int2(k,j) = OUT_p{1}(1,5);
Cond_AR2(k,j) = OUT_p{1}(1,6);

end


% % %%%%%=======================visualization=======
% First annotation
label1 = 'enh-C';
position1 = [Cond_xyz1(k,1,j)/xpixel Cond_xyz1(k,2,j)/ypixel 5];

% Second annotation
label2 = 'prom-C';
position2 = [Cond_xyz2(k,1,j)/xpixel Cond_xyz2(k,2,j)/ypixel 5];

% Combine all positions and labels
positions = [position1; position2];
labels = {label1, label2};

% Create the annotated image
I16bit2 = uint16(65536-inputImage_cond2);
RGB2 = insertObjectAnnotation(I16bit2, "circle", positions, labels);
imshow(RGB2, [], 'InitialMagnification', 800)

% %%%%%================%%%%%================%================
% 

end
end

% % filter out odd z-loc======

for j=1:vis
    for k=2:zs-1
        if abs(Cond_xyz1(k,3,j)-Cond_xyz1(k-1,3,j)) > thr_z || Cond_xyz1(k,3,j) < 0
            Cond_xyz1(k,3,j) = (Cond_xyz1(k-1,3,j) + Cond_xyz1(k+1,3,j))/2;
        end
    end
end

for j=1:vis
    for k=2:zs-1
        if abs(Cond_xyz2(k,3,j)-Cond_xyz2(k-1,3,j)) > thr_z || Cond_xyz2(k,3,j) < 0
            Cond_xyz2(k,3,j) = (Cond_xyz2(k-1,3,j) + Cond_xyz2(k+1,3,j))/2;
        end
    end
end


% ===================% Enhancer-promoter separation in nm %==============% %===================
E_P_dist_3d=zeros(zs,2,vis);
E_P_dist_2d=zeros(zs,2,vis);

for k=1:zs
for j=1:vis
        E_P_dist_3d(k,1,j) = (((enh_xyz(k,1,j)-prom_xyz(k,1,j))^2 + (enh_xyz(k,2,j)-prom_xyz(k,2,j))^2 + (enh_xyz(k,3,j)-prom_xyz(k,3,j))^2)^0.5)*1000;   % 3D separation ENHANCER -tss in nm
        E_P_dist_3d(k,2,j) = ((  ((C_cent_e(k+1,1,j)-C_cent_p(k+1,1,j))*xpixel)^2 + ((C_cent_e(k+1,2,j)-C_cent_p(k+1,2,j))*ypixel)^2 + (enh_xyz(k,3,j)-prom_xyz(k,3,j))^2)^0.5)*1000;   % 3D alternate measurement
        E_P_dist_2d(k,1,j) = (((enh_xyz(k,1,j)-prom_xyz(k,1,j))^2 + (enh_xyz(k,2,j)-prom_xyz(k,2,j))^2)^0.5)*1000;   % 2D separation ENHANCER -tss in nm
        E_P_dist_2d(k,2,j) = ((  ((C_cent_e(k+1,1,j)-C_cent_p(k+1,1,j))*xpixel)^2 + ((C_cent_e(k+1,2,j)-C_cent_p(k+1,2,j))*ypixel)^2)^0.5)*1000;   % 2D alternate measurement
end
end

% %===================% %===================% %===================% %===================

% %==================% %==================% %==================% %==================% %===================
% %===================% Enhancer/Promoter- matched CONDENSATE 3D separation in nm %===================

E_C_dist1=zeros(k,2,j);
P_C_dist1=zeros(k,2,j);

for j=1:vis
    for k=1:zs
        E_C_dist1(k,1,j) = (((enh_xyz(k,1,j)-Cond_xyz1(k,1,j))^2 + (enh_xyz(k,2,j)-Cond_xyz1(k,2,j))^2 + (enh_xyz(k,3,j)-Cond_xyz1(k,3,j))^2)^0.5)*1000;   % 3D separation ENHANCER - enhancer CONDENSATE in nm
        E_C_dist1(k,2,j) = (((enh_xyz(k,1,j)-Cond_xyz1(k,1,j))^2 + (enh_xyz(k,2,j)-Cond_xyz1(k,2,j))^2 )^0.5)*1000;   % 2D separation ENHANCER -CONDENSATE in n

        P_C_dist1(k,1,j) = (((prom_xyz(k,1,j)-Cond_xyz2(k,1,j))^2 + (prom_xyz(k,2,j)-Cond_xyz2(k,2,j))^2 + (prom_xyz(k,3,j)-Cond_xyz2(k,3,j))^2)^0.5)*1000;   % 3D separation PROMOTER - promoter CONDENSATE in nm
        P_C_dist1(k,2,j) = (((prom_xyz(k,1,j)-Cond_xyz2(k,1,j))^2 + (prom_xyz(k,2,j)-Cond_xyz2(k,2,j))^2 )^0.5)*1000;   % 2D separation PROMOTER -CONDENSATE in nM

    end
end

% % % ================%================%================%=======% %===================% %===================% %===========
% %===================% Enhancer/Promoter- other CONDENSATE 3D separation in nm %===================% %===================

E_C_dist2=zeros(k,2,j);
P_C_dist2=zeros(k,2,j);

for j=1:vis
    for k=1:zs
        E_C_dist2(k,1,j) = (((enh_xyz(k,1,j)-Cond_xyz2(k,1,j))^2 + (enh_xyz(k,2,j)-Cond_xyz2(k,2,j))^2 + (enh_xyz(k,3,j)-Cond_xyz2(k,3,j))^2)^0.5)*1000;   % 3D separation ENHANCER - promoter CONDENSATE in nm
        E_C_dist2(k,2,j) = (((enh_xyz(k,1,j)-Cond_xyz2(k,1,j))^2 + (enh_xyz(k,2,j)-Cond_xyz2(k,2,j))^2 )^0.5)*1000;   % 2D separation ENHANCER -CONDENSATE in n

        P_C_dist2(k,1,j) = (((prom_xyz(k,1,j)-Cond_xyz1(k,1,j))^2 + (prom_xyz(k,2,j)-Cond_xyz1(k,2,j))^2 + (prom_xyz(k,3,j)-Cond_xyz1(k,3,j))^2)^0.5)*1000;   % 3D separation PROMOTER - enhancer CONDENSATE in nm
        P_C_dist2(k,2,j) = (((prom_xyz(k,1,j)-Cond_xyz1(k,1,j))^2 + (prom_xyz(k,2,j)-Cond_xyz1(k,2,j))^2 )^0.5)*1000;   % 2D separation PROMOTER -CONDENSATE in nM

    end
end


save(save_filename, 'enh_xyz', 'prom_xyz','TSS_xyz','C_cent_e','C_cent_p','C_cent_TSS','MS2_score','E_P_dist_3d', 'E_P_dist_2d',"r2", 'Cond_xyz1' ,'Cond_xyz2', 'Cond_s1','Cond_s2', 'Cond_int1','Cond_int2', "E_C_dist1",'E_C_dist2', "P_C_dist1",'P_C_dist2')

toc
