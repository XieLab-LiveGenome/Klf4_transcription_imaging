% % % % Enhancer-Promoter_MS2_TRACKING_PIPELINE.m

% % %  3-channel SR-SIM timecourse (5 second frame-frame time interval):
% % %    Channel 1 = MS2        used for transcription state(bursting/inactive)
% % %    Channel 2 = Enhancer   used for enhancer 3D location
% % %    Channel 3 = Promoter   used for promoter 3D location
% % %
% % %  Pipeline:
% % %    Step 1 : Load CZI (all channels) 
% % %    Step 2 : Enhancer 3D tracking starting from seed co-ordinates
% % %    Step 3 : Promoter 3D tracking using enhancer co-ordinates as search center
% % %    Step 4 : MS2 spot tracking and intensity extraction starting
% % %    Step 5 : HMM fitting of MS2 spot intensity for transcription state
% % %    Step 6 : Enhancer-Promoter 3D distance calculation
% % %    Step 7:  Save .mat
% % %
% % %  Dependencies on path:
% % %  ReadImage6D2, track_spot2D, pkfnd, single_zstack,fit_Gaussian3D, SNR_inc2, HMM_fit_fun, binary

clear
clc
tic

% %============================== PROJECT DIRECTORY SETUP=============================
here = fileparts(mfilename('fullpath'));   % Folder of this main script
addpath(genpath(fullfile(here, 'functions')));   % Folder with the functions
addpath(genpath(fullfile(here, 'HMM fitting')));  % Folder with the HMM fitting package
addpath(fullfile(here, 'bioformats'));     % Folder with the bioformats package for image loading
javaaddpath(fullfile(here, 'bioformats', 'bioformats_package.jar'));    % Add the JAVA path for loading bioformats
% %==============================%%==============================%%===================

% % %=========USER DEFINED PARAMETERS======================
C_in =[193 189];   % Approximate pixel co-ordinates of enhancer in start frame/serves as trajectory seed

start_f = 1;       % Start frame
stop_f = 96;       % Stop frame

save_filename='ctrl 3 8-1-25 cell 1 updated code v9';
Input_zstack  = '/Volumes/My Passport/mNG MS2 EP 8-1-25/ctrl 3 8-1-25 cell 1.czi';
MIP_filename='/Volumes/My Passport/mNG MS2 EP 8-1-25/ctrl 3 8-1-25 cell 1_Maximum intensity projection.czi';


channel_enh=2;      % 642 channel
channel_prom=3;      % RED channel
channel_MS2=1;

scene=1;    % set manually

R=50;  %tracking radius for enhancer
R_p = 25; %search radius for promoter
R_m = 20; %search radius for MS2

R_fit = 10;  % 3d gaussian fitting radius
fit_thr = 0.2; % minimum R2 value needed for a good gaussian fit from tracked spot

m_enh = 0.1;
m_prom = 0.1;

%%=================================================

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

zs = min(stop_f,total_f);

vis=size(C_in,1);  % no of EC pairs to be tracked

% % ========Declare variables===========

% % % JF 646 prom is clear so tracked first then enhancer with JF552

enh_xyz=zeros(zs,3,vis);
prom_xyz=zeros(zs,3,vis);

C_cent_e=zeros(zs+1,3,vis);   % the 1st position is the user defined initial condition
C_cent_p=zeros(zs+1,3,vis);

C_cent_TSS = zeros(zs,3,vis);
TSS_xyz=zeros(zs,3,vis);

MS2_score=zeros(zs,5,vis);    % MS2 intensity information

r2=zeros(zs,2,vis);

for k=1:vis
    C_cent_e(start_f,1,k)=C_in(k,1);
    C_cent_e(start_f,2,k)=C_in(k,2);
    C_cent_e(start_f,3,k)= round(z_slice/2);
    C_cent_p(start_f,1,k)=C_in(k,1);
    C_cent_p(start_f,2,k)=C_in(k,2);
    C_cent_p(start_f,3,k)= round(z_slice/2);
end


% 
% % % % % %%==================ENHANCER%=================%=================%=================
% % % % % %%==================ENHANCER%=================%=================%=================
% % 
for j=1:vis

for k=start_f:zs

inputImage_enh = MIP_image6d(1,k,1,channel_enh,:,:);
inputImage_enh = squeeze(inputImage_enh);

inputImage_enh2 = SNR_inc2(inputImage_enh,m_enh);  % denoising through contrast stretching 

out_E = track_spot2D(inputImage_enh2, C_cent_e(k,1,j), C_cent_e(k,2,j),R, m_enh);      % track and 2D gaussian fitting based on guess/previous frame position

if out_E(3) >= 1     % at least 1 enhancer spot FOUND in MIP

C_cent_e(k+1,1,j) = out_E(1);
C_cent_e(k+1,2,j) = out_E(2);      % precise 2D position of enhancer

stack_e = single_zstack(full_stack,1,k,channel_enh);

guess_e = [C_cent_e(k+1,1,j), C_cent_e(k+1,2,j), C_cent_e(k,3,j)];     % guess latest 2D co-ordinates and Z co-ordinate from last frame

[x_ref_e, y_ref_e, z_ref_e, params_e,R2_e] = fit_Gaussian3D(stack_e, guess_e, R_fit, spacing);   % 3D gaussian fitting for precise XYZ  co-ordinates

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

end

if k>start_f
if out_E(3) == 0 || r2(k,1,j) <= fit_thr    % no enhancer spots found or bad fit from tracked spot

guess_e = [enh_xyz(k-1,1,j)/xpixel, enh_xyz(k-1,2,j)/ypixel, enh_xyz(k-1,3,j)/zpixel + 1]; 

[x_ref_e, y_ref_e, z_ref_e, params,R2_e] = fit_Gaussian3D(stack_e, guess_e, R, spacing);     % much wider area for fitting based on last position

enh_xyz(k,1,j) = x_ref_e*xpixel;
enh_xyz(k,2,j) = y_ref_e*ypixel;
enh_xyz(k,3,j) = (z_ref_e-1)*zpixel;

C_cent_e(k+1,1,j) = x_ref_e;
C_cent_e(k+1,2,j) = y_ref_e;
C_cent_e(k+1,3,j) = z_ref_e;

end
end

fprintf('Analyzing frame %d\n', k);    % Frame counter

% %%%%%=======================visualization=======
label='enh';
position = [C_cent_e(k+1,1,j) C_cent_e(k+1,2,j) 5];

if ~isnan(position(1,1))

I16bit=uint16(65536-inputImage_enh2);

RGB = insertObjectAnnotation(I16bit,"circle",position,label);

imshow(RGB,[],'InitialMagnification', 800)

end
%%=============

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

if out_P(3) >= 1       % at least 1 promoter SPOT FOUND in MIP

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

end


if out_P(3) == 0 || r2(k,2,j) <= fit_thr       %%% this is when even 2D tracking fails or very bad gaussian fit

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


fprintf('Analyzing frame %d\n', k);    % Frame counter

%%%%%=======================visualization=======
label='prom';
position = [C_cent_p(k+1,1,j) C_cent_p(k+1,2,j) 5];

if ~isnan(position(1,1))

I16bit=uint16(65536-inputImage_prom2);

RGB = insertObjectAnnotation(I16bit,"circle",position,label);

imshow(RGB,[],'InitialMagnification', 800)
end
%%=============

end
end


% % % filter out odd z-spike movements======

thr_z = 0.5; %%(500 nm z-movement in ~20 sec is most likely due to bad fitting in the frame)

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
% % % ==============


m_MS2 = 0.2;
% % % ================MS2 check=============================================
% % % ================%================%================%================
% 
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

% % % ================Hidden Markov model (2-state fitting) =================
for j=1:vis
inter = HMM_fit_fun (MS2_score(:,2,j));
MS2_score(:,4,j) = inter{1,1};  % 2 state binary model of each trajectory
end

for j=1:vis
MS2_score(:,3,j) = binary (MS2_score(:,4,j));  % 1 or 0 value assigned to ON/OFF states model of each trajectory
end

for j=1:vis
I_max = max(MS2_score(:,2,j));
MS2_score(:,5,j) = MS2_score(:,2,j)/I_max;
end

% % ===================% Enhancer-promoter separation in nm %=====
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

% ================%================%================%========

save(save_filename, 'prom_xyz', 'enh_xyz','TSS_xyz','C_cent_e','C_cent_p','C_cent_TSS','MS2_score','E_P_dist_3d', 'E_P_dist_2d',"r2");

toc

