% % % % Enhancer-MS2-BRD4 hub_TRACKING_PIPELINE.m

% % %  3-channel SR-SIM timecourse:
% % %    Channel 1 = MS2        used for transcription state(bursting/inactive)
% % %    Channel 2 = Enhancer   used for enhancer 3D location
% % %    Channel 3 = Protein hubs     used for searching hubs in the enhancer proximity
% % %
% % %  Pipeline:
% % %    Step 1 : Load CZI (all channels) 
% % %    Step 2 : Enhancer 3D tracking starting from seed co-ordinates
% % %    Step 3 : MS2 spot tracking starting from enhancer co-ordinates
% % %    Step 4 : HMM fitting of MS2 spot intensity for transcription state
% % %    Step 5 : Condensate/Hub peak detection in enhancer local neighborhood
% % %    Step 6 : Primary/secondary/tertiary hub 3D co-ordinates
% % %    Step 7 : Enhancer-Primary/secondary/tertiary hub 3D distance calculation
% % %    Step 8:  Save .mat
% % %
% % %  Dependencies on path:
% % %  ReadImage6D2, SNR_inc2, track_spot2D, single_zstack, fit_Gaussian3D, maskavg, HMM_fit_fun, binary, condensate_search_v6


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


save_filename='Condensate parameters Brd4 ctrl 1 1-8-2026 cell 1 3x thr 1500 v7 code';

C_in=[99 96];     % Approximate pixel co-ordinates of enhancer in start frame/serves as trajectory seed

start_f = 1;      % START frame
stop_f = 50;      % STOP frame

% % % ---Read CZI files of whole stack and MIP---------------------------------------------------------------------------------
% % -------------------------------------------------------------------------------------------------------------------------
%%===========Video read==========
Input_zstack='/Volumes/My Passport/Brd4 MS2 1-8-2026/off-on/ctrl 1 1-8-2026 cell 1.czi';
MIP_filename='/Volumes/My Passport/Brd4 MS2 1-8-2026/off-on/ctrl 1 1-8-2026 cell 1_Maximum intensity projection.czi';
scene=1;    % set manually

MIP_out = ReadImage6D2(MIP_filename, true, scene);
stack_out = ReadImage6D2(Input_zstack, true, scene);
metadata = stack_out{2};
full_stack = stack_out{1};
MIP_image6d = MIP_out{1};

%%===========DNA/RNA/PROTEIN labeling channels==========

channel_MS2=1;  % MS2, mStaygold
channel_enh=2;  % Enhancer channel, SNAPTag JF552
channel_cond=3;  % Condensate/Hub channel, HaloTag 642

% % % -----------------------Imaging parameters load from stack----------------------------------------------------------------
% % -------------------------------------------------------------------------------------------------------------------------

xpixel = 0.0313;    % pixel scaling, um/pixel
ypixel = 0.0313;   
zpixel = metadata.ScaleZ;

spacing = [xpixel ypixel zpixel];

z_slice = metadata.SizeZ; % number of z-slices
total_f = metadata.SizeT;  % total number of frames

%%====================================

zs = min(stop_f,total_f); 

vis=size(C_in,1);  % no of EC pairs to be tracked

sd = zeros(zs,vis);

%========Declare variables==================
%%==========================================

C_cent_e=zeros(zs+1,3,vis);         % integer co-ordinates in voxels
C_cent_cond=zeros(zs,3,vis);
C_cent_TSS=zeros(zs,3,vis);

enh_xyz=zeros(zs,3,vis);          % actual co-ordinates in microns
TSS_xyz=zeros(zs,3,vis);

MS2_score=zeros(zs,5,vis);    % MS2 intensity information

r2=zeros(zs,2,vis);

% =====primary(1st)======
Cond_xyz1=zeros(zs,3,vis);  % xyz co-ordinates 
Cond_s1=zeros(zs,vis);   %  SIZE OF CONDENSATE in MIP
Cond_int1=zeros(zs,vis);  % BRIGHTNESS OF CONDENSATE in MIP
Cond_AR1=zeros(zs,vis);   % AR OF CONDENSATE in MIP

% =====secondary(2nd)======
Cond_xyz2=zeros(zs,3,vis);  % xyz co-ordinates 
Cond_s2=zeros(zs,vis);   %  SIZE OF CONDENSATE in MIP
Cond_int2=zeros(zs,vis);  % BRIGHTNESS OF CONDENSATE in MIP
Cond_AR2=zeros(zs,vis);   % AR OF CONDENSATE in MIP

% =====tertiary(3rd)======
Cond_xyz3=zeros(zs,3,vis);  % xyz co-ordinates 
Cond_s3=zeros(zs,vis);   %  SIZE OF CONDENSATE in MIP
Cond_int3=zeros(zs,vis);  % BRIGHTNESS OF CONDENSATE in MIP
Cond_AR3=zeros(zs,vis);   % AR OF CONDENSATE in MIP

%%====================================
%%====================================

R=40;       % SEARCH RADIUS FOR NEXT TIMEPOINT
R_fit = 10;  % 3d gaussian fitting radius
fit_thr = 0.2; % minimum R2 value needed for a good gaussian fit from tracked spot

m_enh = 0.1;

% % % % % ========END OF variables===========

% %%===========Enahncer image extract and process============
% 

for j=1:vis
    C_cent_e(start_f,1,j)=C_in(j,1);
    C_cent_e(start_f,2,j)=C_in(j,2);
    C_cent_e(start_f,3,j)= round(z_slice/2);
end

figure 
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

% %%%%%=======================visualization===============
label='enh';
position = [C_cent_e(k+1,1,j) C_cent_e(k+1,2,j) 5];

if ~isnan(position(1,1))

I16bit=uint16(65536-inputImage_enh2);

RGB = insertObjectAnnotation(I16bit,"circle",position,label);

imshow(RGB,[],'InitialMagnification', 800)

end
%%=============%%=============%%=============%%=============

end
end


% % filter out odd z-spikes======
thr_z = 0.5;

for j=1:vis
    for k=2:zs-1
        if abs(enh_xyz(k,3,j)-enh_xyz(k-1,3,j)) > thr_z || enh_xyz(k,3,j) < 0
            enh_xyz(k,3,j) = (enh_xyz(k-1,3,j) + enh_xyz(k+1,3,j))/2;
        end
    end
end

R_m = 30;
m_MS2 = 0.1;
% % % % ================MS2 check=============================================
% % % % ================%================%================%================
% % 
for j=1:vis
for k= start_f:zs

inputImage_MS2 = MIP_image6d(1,k,1,channel_MS2,:,:); 
inputImage_MS2 = squeeze(inputImage_MS2);  

out_M = track_spot2D(inputImage_MS2, C_cent_e(k+1,1,j), C_cent_e(k+1,2,j),R_m,m_MS2);

guess_m = [out_M(1), out_M(2), C_cent_e(k+1,3,j)];

stack_M = single_zstack(full_stack,1,k,channel_MS2);  % Channel, Timepoint

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
MS2_score(k,2,j) = maskavg(inputImage_MS2,C_cent_e(k+1,1,j),C_cent_e(k+1,2,j),5);
end

end
end

% % ================Hidden Markov model (2-state fitting) =================
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


% % % % % % % % % ================% % ================
% % % % % ==================CONDENSATE/HUB positions=================%=================%=================
% % % % % ==================CONDENSATE/HUB positions=================%=================%=================
% % 
figure 

nuc_ave = 1500; % average nuclear background
mult = 3.0;  % hub detection threshold = mult*average nuclear background

for j=1:vis
for k=start_f:zs 

inputImage_cond = MIP_image6d(1,k,1,channel_cond,:,:);
inputImage_cond = squeeze(inputImage_cond);

inputImage_cond2 = SNR_inc2(inputImage_cond,0.1);     % BACKGROUND REDUCTION AND SIGNAL BOOSTING

sd(k,j)=std(inputImage_cond(:));

t_cond=nuc_ave*sd(k,j)/sd(start_f,j);   % adpative thresolding based on whole image brightness

BW = imbinarize(inputImage_cond,t_cond*mult);

BW2 = bwareafilt(BW,[10 inf]);   % filtering very small condensates 

BW3 = inputImage_cond().*BW2;

% -------------% -------------% -------------% -------------% -------------% -------------% -------------

stack_C = single_zstack(full_stack,1,k,channel_cond);

R_search_c  = 10; 

enh_pos = [enh_xyz(k,1,j)  enh_xyz(k,2,j)  enh_xyz(k,3,j)];

% -------------% -------------% -------------% -------------% -------------% -------------% -------------

OUT = condensate_search_v6(stack_C, BW3, enh_pos, spacing, t_cond, mult, R_search_c);


%%========Primary(1st) hub==============
Cond_xyz1(k,1,j) = OUT{1}(1,1);
Cond_xyz1(k,2,j) = OUT{1}(1,2);
Cond_xyz1(k,3,j) = OUT{1}(1,3);
Cond_s1(k,j) = OUT{1}(1,4);
Cond_int1(k,j) = OUT{1}(1,5);
Cond_AR1(k,j) = OUT{1}(1,6);

%%========Primary(2nd) hub==============
Cond_xyz2(k,1,j) = OUT{2}(1,1);
Cond_xyz2(k,2,j) = OUT{2}(1,2);
Cond_xyz2(k,3,j) = OUT{2}(1,3);
Cond_s2(k,j) = OUT{2}(1,4);
Cond_int2(k,j) = OUT{2}(1,5);
Cond_AR2(k,j) = OUT{2}(1,6);

%%========Primary(3rd) hub==============
Cond_xyz3(k,1,j) = OUT{3}(1,1);
Cond_xyz3(k,2,j) = OUT{3}(1,2);
Cond_xyz3(k,3,j) = OUT{3}(1,3);
Cond_s3(k,j) = OUT{3}(1,4);
Cond_int3(k,j) = OUT{3}(1,5);
Cond_AR3(k,j) = OUT{3}(1,6);

% 
% % %%%%%=======================visualization=======
% First annotation
label1 = '1st';
position1 = [Cond_xyz1(k,1,j)/xpixel Cond_xyz1(k,2,j)/ypixel 5];

% Second annotation
label2 = '2nd';
position2 = [Cond_xyz2(k,1,j)/xpixel Cond_xyz2(k,2,j)/ypixel 5];

% Third annotation
label3 = '3rd';
position3 = [Cond_xyz3(k,1,j)/xpixel Cond_xyz3(k,2,j)/ypixel 5];

% Combine all positions and labels
positions = [position1; position2; position3];
labels = {label1, label2, label3};

% Create the annotated image
I16bit2 = uint16(65536-inputImage_cond2);
RGB2 = insertObjectAnnotation(I16bit2, "circle", positions, labels);
imshow(RGB2, [], 'InitialMagnification', 800)
% %%%%%================%%%%%================
% 

end
end

% %===================% Enhancer-TSS separation in nm %===================% %===================

E_TSS_dist=zeros(k,2,j);

for j=1:vis
    for k=1:zs
    if MS2_score(k,1,j)==1
        E_TSS_dist(k,1,j) = (((enh_xyz(k,1,j)-TSS_xyz(k,1,j))^2 + (enh_xyz(k,2,j)-TSS_xyz(k,2,j))^2 + (enh_xyz(k,3,j)-TSS_xyz(k,3,j))^2)^0.5)*1000;   % 3D separation ENHANCER -tss in nm
        E_TSS_dist(k,2,j) = (((enh_xyz(k,1,j)-TSS_xyz(k,1,j))^2 + (enh_xyz(k,2,j)-TSS_xyz(k,2,j))^2 )^0.5)*1000;   % 2D separation ENHANCER -tss in nm
    end
    end
end

% %==================% %==================% %==================% %==================
% %===================% Enhancer-CONDENSATE 3d separation in nm %===================% %===================

E_C_dist1=zeros(k,2,j);
E_C_dist2=zeros(k,2,j);
E_C_dist3=zeros(k,2,j);

for j=1:vis
    for k=1:zs
        E_C_dist1(k,1,j) = (((enh_xyz(k,1,j)-Cond_xyz1(k,1,j))^2 + (enh_xyz(k,2,j)-Cond_xyz1(k,2,j))^2 + (enh_xyz(k,3,j)-Cond_xyz1(k,3,j))^2)^0.5)*1000;   % 3D separation ENHANCER -CONDENSATE in nm
        E_C_dist1(k,2,j) = (((enh_xyz(k,1,j)-Cond_xyz1(k,1,j))^2 + (enh_xyz(k,2,j)-Cond_xyz1(k,2,j))^2 )^0.5)*1000;   % 2D separation ENHANCER -CONDENSATE in n

        E_C_dist2(k,1,j) = (((enh_xyz(k,1,j)-Cond_xyz2(k,1,j))^2 + (enh_xyz(k,2,j)-Cond_xyz2(k,2,j))^2 + (enh_xyz(k,3,j)-Cond_xyz2(k,3,j))^2)^0.5)*1000;   % 3D separation ENHANCER -CONDENSATE in nm
        E_C_dist2(k,2,j) = (((enh_xyz(k,1,j)-Cond_xyz2(k,1,j))^2 + (enh_xyz(k,2,j)-Cond_xyz2(k,2,j))^2 )^0.5)*1000;   % 2D separation ENHANCER -CONDENSATE in n

        E_C_dist3(k,1,j) = (((enh_xyz(k,1,j)-Cond_xyz3(k,1,j))^2 + (enh_xyz(k,2,j)-Cond_xyz3(k,2,j))^2 + (enh_xyz(k,3,j)-Cond_xyz3(k,3,j))^2)^0.5)*1000;   % 3D separation ENHANCER -CONDENSATE in nm
        E_C_dist3(k,2,j) = (((enh_xyz(k,1,j)-Cond_xyz3(k,1,j))^2 + (enh_xyz(k,2,j)-Cond_xyz3(k,2,j))^2 )^0.5)*1000;   % 2D separation ENHANCER -CONDENSATE in n
    end
end

% %===================% TSS-condensate separation in nm %===================% %===================

TSS_cond_dist=zeros(k,2,j);

for j=1:vis
    for k=1:zs
    if MS2_score(k,1,j)==1
        TSS_cond_dist(k,1,j) = (((Cond_xyz1(k,1,j)-TSS_xyz(k,1,j))^2 + (Cond_xyz1(k,2,j)-TSS_xyz(k,2,j))^2 + (Cond_xyz1(k,3,j)-TSS_xyz(k,3,j))^2)^0.5)*1000;   % 3D separation ENHANCER -tss in nm
        TSS_cond_dist(k,2,j) = (((Cond_xyz1(k,1,j)-TSS_xyz(k,1,j))^2 + (Cond_xyz1(k,2,j)-TSS_xyz(k,2,j))^2)^0.5)*1000; 
    end
    end
end


save(save_filename, 'MS2_score', 'enh_xyz', 'TSS_xyz', 'Cond_xyz1' ,'Cond_xyz2','Cond_xyz3', 'Cond_s1','Cond_s2', 'Cond_s3', 'Cond_int1','Cond_int2','Cond_int3', "E_TSS_dist", "E_C_dist1",'E_C_dist2','E_C_dist3', "TSS_cond_dist")

toc

