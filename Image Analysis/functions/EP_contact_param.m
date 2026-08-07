function out = EP_contact_param(EP_distances, thr, time_int)

out = cell(4,1);

multiple_trajectories = EP_distances';

[num_trajectories, num_time] = size(multiple_trajectories);

vis = num_trajectories;

on_times=cell(vis,3);       % duration of ON state (EP looping)
off_times=cell(vis,3);      % duration of OFF state (EP unlooped)
on_cellavg = zeros(vis,1);  % Average ON duration for 1 trajectory (cell)
off_cellavg = zeros(vis,1);       % Average OFF duration for 1 trajectory (cell)

for j=1:vis
traj = (multiple_trajectories(j,:))';

traj_n = traj_clean(traj);  % NaN values removed from end and gaps taken care of

binary_2s(traj_n>thr) = 0;
binary_2s(traj_n<=thr) = 1;
v = binary_2s;

% % % % =============% % % ===ON & OFF duration calculations=============

[starts1, ends1, lengths1, starts0, ends0, lengths0] = binaryv2(v);

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
% 
vec1 = cellfun(@(x) x(:), on_times(:,3), 'UniformOutput', false); % Flatten each element
on_compiled = vertcat(vec1{:});  % compilation of ON durations

vec0 = cellfun(@(x) x(:), off_times(:,3), 'UniformOutput', false); % Flatten each element
off_compiled = vertcat(vec0{:});  % compilation of OFF durations

out{1} = on_compiled;
out{2} = off_compiled;
out{3} = on_cellavg;
out{4} = off_cellavg;

end