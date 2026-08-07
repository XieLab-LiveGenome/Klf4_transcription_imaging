clear
clc

file_path = '/Users/janaa/Desktop/MS2 transcription/EP datasheets/Fast EP /EP/enhancer_promoter fast dynamics_8-18-25.xlsx';   %Excel sheet fielpath here
sheet_name = 'Distance compiled_ctrl_3c';

% Read data
ec_dist = readmatrix(file_path, 'Sheet', sheet_name, 'Range', 'B2:AD252');
multiple_trajectories = ec_dist';

[num_trajectories, num_time] = size(multiple_trajectories);

% Define bin edges
edges = [0, 120, 360, 600, Inf];  % 4 bins

% Assign bin indices: 1 to 4; NaN stays NaN
bin_idx = discretize(multiple_trajectories, edges);  % result: NaNs preserved

% Map bins to 1–4, NaNs to 5 (for color index)
color_data = bin_idx;
color_data(isnan(bin_idx)) = 5;

% Define colormap with 5 colors: 4 for bins, 1 for NaN (gray)
cmap = [1 0 0;         % Red
        0.4 0.4 0.9;   % Light blue
        0.2 0.2 1;     % Blue
        0 0 1;         % Dark blue
        0.6 0.6 0.6];  % Gray for NaNs



% Create scaled x-axis values: 1 point = 5 units
x_units = linspace(0, (num_time - 1) * 5, num_time);

% Plot image with indexed colors
figure
h = imagesc(x_units, 1:num_trajectories, color_data);   % Data values are 1–5
colormap(cmap);
colorbar;

% Manually set colorbar ticks to match bins + NaN
cb = colorbar;
cb.Ticks = [1 2 3 4 5];
cb.TickLabels = {'<120', '120–360', '360–600', '>600', 'NaN'};

% Label axes
xlabel('Time (seconds)', 'FontWeight', 'bold', 'FontSize', 20);
ylabel('Trajectory index', 'FontWeight', 'bold', 'FontSize', 20);
cb.Label.String = '3D EP separation (nm)';
cb.Label.FontWeight = 'bold';
cb.Label.FontSize = 20;

% set(gca, 'FontSize', 20, 'FontWeight', 'bold');

% Ensure axis and color scaling are tight and consistent
caxis([1 5]);  % map 1–5 directly to colormap rows
