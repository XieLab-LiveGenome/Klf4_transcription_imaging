% % % % Enhancer-MS2-BRD4 hub_configuration_PIPELINE.m

% % %  3-channel SR-SIM timecourse:
% % %    Channel 1 = MS2        used for transcription state(bursting/inactive)
% % %    Channel 2 = Enhancer   used for enhancer 3D location
% % %    Channel 3 = Protein hubs     used for searching hubs in the enhancer proximity
% % %
% % %  Pipeline:
% % %    Step 1 : Load CZI (all channels) 
% % %    Step 2 : Refine enhancer 3D position via 3D Gaussian fit
% % %    Step 3 : MS2 spot refined 3D position starting from enhancer co-ordinates and spot intensity extraction
% % %    Step 4 : Nucleus segmentation using Cellpose
% % %    Step 5 : Condensate/Hub peak detection in enhancer local neighborhood in each nucleus
% % %    Step 5 : Primary hub 3D co-ordinates 
% % %    Step 7 : Enhancer-hub 3D distance calculation
% % %    Step 8:  Save .mat
% % %
% % %  Dependencies on path:
% % %  ReadImage6D2, SNR_inc2, track_spot2D, single_zstack, fit_Gaussian3D, maskavg, HMM_fit_fun, binary, condensate_search1 , point_location2D

clear
clc
tic

% %============================== PROJECT DIRECTORY SETUP=============================
here = fileparts(mfilename('fullpath'));   % folder this script lives in
addpath(fullfile(here, 'functions'));
addpath(fullfile(here, 'bioformats'));
javaaddpath(fullfile(here, 'bioformats', 'bioformats_package.jar'));
% %==============================%%==============================%%===================

% % ---Read CZI files of whole stack and MIP---------------------------------------------------------------------------------
% % % -------------------------------------------------------------------------------------------------------------------------
% % ===========Video read==========
save_filename='/Users/janaa/Desktop/MS2 transcription/Compiled data/Condensate E MS2 imaging -3 color/BRD4/Condensate parameters brg1 10 um 2.5 hr 1 BRD4 100 ms_SIM.mat';
Input_zstack='/Volumes/Aniket2/4c MED14 dtag 4-22-2026/M14 enh BRD4/m14 dtag 1.5 hr 1_SIM.czi';
MIP_filename='/Volumes/Aniket2/4c MED14 dtag 4-22-2026/M14 enh BRD4/m14 dtag 1.5 hr 1_SIM_Maximum intensity projection.czi';

scene=1;    % set manually

MIP_out = ReadImage6D2(MIP_filename, true, scene);
stack_out = ReadImage6D2(Input_zstack, true, scene);
metadata = stack_out{2};
full_stack = stack_out{1};
MIP_image6d = MIP_out{1};

% % ===========Video read END==========

channel_MS2=1;  % MS2 CHANNEL, GF
channel_enh=2;  % ENH CHANNEL, RED snap 552
channel_cond=3;  % condensate channel, 642

C_in = [1533.13099	1837.188498
1422.619808	422.5559105
1460.638978	2223.70607
1668.658147	2170.702875
531.5654952	1300.638978
2489.744409	1502.651757
2222.715655	439.5527157
1153.610224	408.5303514];    % Approximate pixel co-ordinates of enhancer (Verify the cell has all 3 channel labels in ImageJ)

% % % % % % -----------------------Imaging parameters load from stack----------------------------------------------------------------
% % % % % % % -------------------------------------------------------------------------------------------------------------------------
% 
xpixel = 0.0313;    % this is the real pixel value
ypixel = 0.0313;   
zpixel = metadata.ScaleZ;

spacing = [xpixel ypixel zpixel];

z_slice = metadata.SizeZ; % number of z-slices
zs = metadata.SizeT;  % number of frames
vis=size(C_in,1);  % no of EC pairs to be tracked

% 
% %========Declare variables===========
C_cent_e=zeros(zs,3,vis);         % integer co-ordinates in voxels
C_cent_cond=zeros(zs,3,vis);
C_cent_TSS=zeros(zs,3,vis);

enh_xyz=zeros(zs,3,vis);          % actual co-ordinates in microns
Cond_xyz=zeros(zs,3,vis);
TSS_xyz=zeros(vis,3);


MS2_score=zeros(zs,2,vis);    % MS2 intensity information

Cond_s=zeros(zs,vis);   %  SIZE OF CONDENSATE in MIP
Cond_int=zeros(zs,vis);  % BRIGHTNESS OF CONDENSATE in MIP
Cond_AR=zeros(zs,vis);   %  SIZE OF CONDENSATE in MIP
% 
% 
R=12;      % enhancer tracking radius
% 
start_f = 1;
% % 
% % % % % % ========END OF variables===========
% % 
% % % % %%===========Enhancer image extract and process============
% % % 
% figure 
for j=1:vis

guess_e = [C_in(j,1), C_in(j,2), round(z_slice/2)];

for k=start_f:zs

inputImage_enh = MIP_image6d(1,k,1,channel_enh,:,:);
inputImage_enh = squeeze(inputImage_enh);

m=size(inputImage_enh,1); % y dimension of image
n=size(inputImage_enh,2); % x dimension of image

stack_e = single_zstack(full_stack,1,k,channel_enh);

[x_ref_e, y_ref_e, z_ref_e, params_e] = fit_Gaussian3D(stack_e, guess_e, R, spacing);

enh_xyz(k,1,j) = x_ref_e*xpixel;
enh_xyz(k,2,j) = y_ref_e*ypixel;
enh_xyz(k,3,j) = (z_ref_e-1)*zpixel;

guess_e = [x_ref_e, y_ref_e, z_ref_e];

C_cent_e(k,1,j) = round(x_ref_e) ;
C_cent_e(k,2,j) = round(y_ref_e) ;
C_cent_e(k,3,j) = round(z_ref_e) ;

%%%%%=======================visualization=======
label='enh';
position = [C_cent_e(k,1,j) C_cent_e(k,2,j) 5];

I16bit=uint16(inputImage_enh);

RGB = insertObjectAnnotation(I16bit,"circle",position,label);

imshow(RGB,[])
%%=============


end
end

% % ================MS2 check=============================================
% % ================%================%================%================

for j=1:vis
for k= 1:zs

inputImage_MS2 = MIP_image6d(1,k,1,channel_MS2,:,:); 
inputImage_MS2 = squeeze(inputImage_MS2); 
t_MS2=0.05*max(inputImage_MS2(:));

masked_img_MS2 = imagemask(inputImage_MS2,C_cent_e(k,1,j),C_cent_e(k,2,j),30);

pk_MS2 = pkfnd(masked_img_MS2,t_MS2,5); % -------------set this-------------------
s=size(pk_MS2,1);

if s>=1

bright_spots_y_n = pk_MS2(:,2);
bright_spots_x_n = pk_MS2(:,1);
num_spots_n=size(pk_MS2,1);
I_spots_n=zeros(num_spots_n,1);

for i=1:num_spots_n
    index1=bright_spots_x_n(i,1);
    index2=bright_spots_y_n(i,1);
I_spots_n(i,1)=inputImage_MS2(index2,index1);
end

b_max_n=max(I_spots_n);

for i=1:num_spots_n
if I_spots_n(i,1)==b_max_n
    tss_n=i;
end
end

TSS_y_n=pk_MS2(tss_n,2);    % approximate spot location 
TSS_x_n=pk_MS2(tss_n,1);

z_m = enh_xyz(k,3,j)/zpixel +1 ;

guess_m = [TSS_x_n, TSS_y_n, z_m];

stack_M = single_zstack(full_stack,1,k,channel_MS2);

[x_ref_m, y_ref_m, z_ref_m, params_m] = fit_Gaussian3D(stack_M, guess_m, R, spacing);

TSS_xyz(k,1,j) = x_ref_m*xpixel;
TSS_xyz(k,2,j) = y_ref_m*ypixel;
TSS_xyz(k,3,j) = (z_ref_m-1)*zpixel;

C_cent_TSS(k,1,j) = round(x_ref_m);
C_cent_TSS(k,2,j) = round(y_ref_m);
C_cent_TSS(k,3,j) = round(z_ref_m);

MS2_score(k,1,j)=1;
MS2_score(k,2,j) = maskavg(inputImage_MS2,C_cent_TSS(k,1,j),C_cent_TSS(k,2,j),5);

elseif s==0
TSS_xyz(k,1,j) = NaN;
TSS_xyz(k,2,j) = NaN;
TSS_xyz(k,3,j) = NaN;
MS2_score(k,1,j)= -1;
MS2_score(k,2,j) = maskavg(inputImage_MS2,C_cent_e(k,1,j),C_cent_e(k,2,j),5);
end

end
end

% % % ==================CONDENSATE%=================%=================%=================
% % % ==================CONDENSATE%=================%=================%=================

mult = 3; 

figure 

for k = 1:zs

inputImage_cond = MIP_image6d(1,k,1,channel_cond,:,:);
inputImage_cond = squeeze(inputImage_cond);
stack_cond = single_zstack(full_stack,1,k,channel_cond);

seg_OUT = cellpose_seg(inputImage_cond,150);  % roughly segment the nuclei
nuc_ave = mean(seg_OUT{1,2}(:));   % average nuclear intensity in image

if isnan(nuc_ave)
    nuc_ave = mean (inputImage_cond(:));
end

for j=1:vis
    x_loc = round (C_cent_e(k+1,1,j));
    y_loc = round (C_cent_e(k+1,2,j));

    id = point_location2D(seg_OUT{1,3}, x_loc, y_loc);   % extracting nucleus id corresponding to enhancer location
    if ~isnan(id)

    t_cond = seg_OUT{1,2}(id,1);  % get intensity of nucleus corresponding to enhancer

    inputImage_cond2 = inputImage_cond().*seg_OUT{1,1}{id,1};   % use that specific nuclear mask

    B = imbinarize(inputImage_cond2,mult*t_cond/2);    % binarize to get rough estimate of condensates

    B2 = bwareafilt(B,[10 inf]);   % filtering very small condensates from this rough binary image

    B3 = inputImage_cond2().*B2;

    else 
        B = imbinarize(inputImage_cond,mult/2*nuc_ave);

        B2 = bwareafilt(B,[10 inf]);   % filtering very small condensates 

        B3 = inputImage_cond().*B2;
        t_cond = nuc_ave;

    end


    R_search_c  = 15; 

enh_pos = [enh_xyz(k,1,j)  enh_xyz(k,2,j)  enh_xyz(k,3,j)];

% -------------% -------------% -------------% -------------% -------------% -------------% -------------

OUT = condensate_search1(stack_cond, B3, enh_pos, spacing, t_cond, mult, R_search_c);

for n=2:20
if isnan(OUT{1})
    OUT = condensate_search1(stack_cond, B3, enh_pos, spacing, t_cond, mult, n*R_search_c);
end
end

Cond_xyz(k,1,j) = OUT{1}(1,1);
Cond_xyz(k,2,j) = OUT{1}(1,2);
Cond_xyz(k,3,j) = OUT{1}(1,3);
Cond_s(k,j) = OUT{1}(1,4)*xpixel*1000;
Cond_int(k,j) = OUT{1}(1,5);
Cond_AR(k,j) = OUT{1}(1,6);

% %%%%%=======================visualization=======
label='condensate';
position = [Cond_xyz(k,1,j)/xpixel Cond_xyz(k,2,j)/ypixel 5];

I16bit2=uint16(65536-B3);

RGB2 = insertObjectAnnotation(I16bit2,"circle",position,label);

imshow(I16bit2,[])
% %%%%%================%%%%%================

end
end

% % %===================% Enhancer-CONDENSATE 3d separation in nm %===================% %===================
% % 
E_C_dist=zeros(zs,2,vis);

for j=1:vis
    for k=1:zs
        E_C_dist(k,1,j) = (((enh_xyz(k,1,j)-Cond_xyz(k,1,j))^2 + (enh_xyz(k,2,j)-Cond_xyz(k,2,j))^2 + (enh_xyz(k,3,j)-Cond_xyz(k,3,j))^2)^0.5)*1000;   % 3D separation ENHANCER -CONDENSATE in nm
        E_C_dist(k,2,j) = (((enh_xyz(k,1,j)-Cond_xyz(k,1,j))^2 + (enh_xyz(k,2,j)-Cond_xyz(k,2,j))^2 )^0.5)*1000;   % 2D separation ENHANCER -CONDENSATE in nm
    end
end
% 

if zs == 1 

enh_xyz = (reshape(enh_xyz,3,[]))';
Cond_xyz = (reshape(Cond_xyz,3,[]))';
TSS_xyz = (reshape(TSS_xyz,3,[]))';
MS2_score = (reshape(MS2_score,2,[]))';
E_C_dist = (reshape(E_C_dist,2,[]))';

end

save(save_filename,'MS2_score','TSS_xyz', 'enh_xyz', 'Cond_xyz', 'Cond_s','Cond_int', 'Cond_AR', "E_C_dist")

toc