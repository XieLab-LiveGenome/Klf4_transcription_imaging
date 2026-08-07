function bur = burst_analysis_MS2(I_t, I_t_norm, time_int)
%BURST_ANALYSIS_MS2  Two-state HMM fitting, ON/OFF durations and burst metrics.
%
%   bur = BURST_ANALYSIS_MS2(I_t, I_t_norm, time_int)
%
%   Inputs:
%       I_t       (zs x vis) raw TS intensity, used for the burst integrals
%       I_t_norm  (zs x vis) SD-corrected intensity, used for the HMM fit
%       time_int  frame interval (durations come out in these units)
%
%   Output struct bur:
%       .fitMS2_2s .binary_2s
%       .on_times .off_times .on_cellavg .off_cellavg .on_compiled .off_compiled
%       .Burst .Burst_amp .burst_cellavg .burst_amp_cellavg
%       .Burst_compiled .Burst_amp_compiled .base_v
%
%   The logic is unchanged from the original inline block; only the variable
%   set-up and the output packing are new.

[zs, vis] = size(I_t);

fitMS2_2s = zeros(zs,vis);  % HMM 2-state fitting
binary_2s = zeros(zs,vis);  % binary 1/0 converted states
on_times  = cell(vis,3);    % duration of ON state (1 burst)
off_times = cell(vis,3);    % duration of OFF state
on_cellavg  = zeros(vis,1);       % Average ON duration for 1 trajectory (cell)
off_cellavg = zeros(vis,1);       % Average OFF duration for 1 trajectory (cell)
burst_cellavg = zeros(vis,1);     % Average BURST SIZE for 1 trajectory (cell)
burst_amp_cellavg = zeros(vis,1); % Average BURST AMPLITUDE for 1 trajectory (cell)


% ================Hidden Markov model (2-state fitting) =================
for j=1:vis
inter = HMM_fit_fun (I_t_norm(:,j));
% inter = HMM_fit_fun (I_t(:,j));
fitMS2_2s(:,j) = inter{1,1};  % 2 state binary model of each trajectory
end

for j=1:vis
binary_2s(:,j) = binary (fitMS2_2s(:,j));  % 1 or 0 value assigned to ON/OFF states model of each trajectory
end

% % % ===============% % % ===============% % % =======================

% % % =============% % % ===ON & OFF duration calculations=============


for j = 1:vis

v=(binary_2s(:,j))';

[starts1, ends1, lengths1, starts0, ends0, lengths0] = binaryv2(v);

% if size(ends1,2) >= 2
%
% if starts1(1)==1
%     starts1(1)=[];
%     ends1(1)=[];
%     lengths1(1)=[];  %% filter out incomplete ON at start
% end
% if ends1(end)==zs
%      starts1(end)=[];
%      ends1(end)=[];
%      lengths1(end)=[];  %% filter out incomplete ON at end
% end
%
% end
%
%
% if starts0(1)==1
%     starts0(1)=[];
%     ends0(1)=[];
%     lengths0(1)=[];     %% filter out incomplete OFF at start
% end
% if ends0(end)==zs
%      starts0 = starts0(1:end-1);
%      ends0 = ends0(1:end-1);
%      lengths0 = lengths0(1:end-1);    %% filter out incomplete OFF at end
% end

on_times{j,1} = starts1;
on_times{j,2} = ends1;
on_times{j,3} = (lengths1)*time_int;

off_times{j,1} = starts0;
off_times{j,2} = ends0;
off_times{j,3} = (lengths0)*time_int;
end

for j = 1:vis
on_cellavg(j,1) = mean(on_times{j,3});
off_cellavg(j,1) = mean(off_times{j,3});
end

vec1 = cellfun(@(x) x(:), on_times(:,3), 'UniformOutput', false); % Flatten each element
on_compiled = vertcat(vec1{:});  % compilation of ON durations

vec0 = cellfun(@(x) x(:), off_times(:,3), 'UniformOutput', false); % Flatten each element
off_compiled = vertcat(vec0{:});  % compilation of Off durations


% % % ===================================================================
% % % ===============Burst amplitude and size calculations===============

num_t=vis;
length_t=zs;


base=zeros(length_t,num_t);
base_v=zeros(1,num_t);

size_B=50;

B=zeros(size_B,num_t);
Burst=zeros(size_B,num_t);      % burst size
Burst_amp=zeros(size_B,num_t);  % burst amplitude

%%%=======Actual Intensity at OFF==================
for j=1:num_t       % cycle through all trajectories

for i=1:length_t
if binary_2s(i,j)==0
    base(i,j)= I_t (i,j);
end
end

base_v(1,j) = mean(nonzeros(base(:,j)));

end

%%%====================Burst size & amplitude calculations================

for j=1:num_t

on_s1 = on_times{j,1};
on_s2 = on_times{j,2};


if ~isnan(on_s1)

m=size(nonzeros(on_s1),1);

for i=1:m
    s1=on_s1(1,i);
    s2=on_s2(1,i);
    for p=s1:s2
        B(i,j)=B(i,j)+I_t(p,j);
    end
    Burst(i,j)=(B(i,j)-(s2-s1+1)*base_v(1,j))*time_int;  %%  ith burst of jth trajectory

    Burst_amp(i,j)=(Burst(i,j)/(s2-s1+1))/time_int;

end

end

end

for j=1:vis
burst_cellavg(j,1) = mean(nonzeros(Burst(:,j)));
burst_amp_cellavg(j,1) = mean(nonzeros(Burst_amp(:,j)));
end

Burst_compiled=nonzeros(Burst);          % all bursts from all trajectories
Burst_amp_compiled=nonzeros(Burst_amp);  % all burst amplitudes from all trajectories

% % % ===================================================================
% % % =========================PACK OUTPUT===============================

bur.fitMS2_2s          = fitMS2_2s;
bur.binary_2s          = binary_2s;
bur.on_times           = on_times;
bur.off_times          = off_times;
bur.on_cellavg         = on_cellavg;
bur.off_cellavg        = off_cellavg;
bur.on_compiled        = on_compiled;
bur.off_compiled       = off_compiled;
bur.Burst              = Burst;
bur.Burst_amp          = Burst_amp;
bur.burst_cellavg      = burst_cellavg;
bur.burst_amp_cellavg  = burst_amp_cellavg;
bur.Burst_compiled     = Burst_compiled;
bur.Burst_amp_compiled = Burst_amp_compiled;
bur.base_v             = base_v;

end
