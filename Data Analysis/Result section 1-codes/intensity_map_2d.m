clear
clc

file_path = '/Users/janaa/Desktop/MS2 transcription/MS2 data sheets/Master sheets/Latest/Unedited master sheet latest 5-17-25.xlsx';

CTRL_TRAJ = readmatrix(file_path,'Sheet','PARENTAL CTRL TRAJ');
n= 121;  % max timepoints

CTRL_TRAJ = CTRL_TRAJ(1:n, :);  % if some trajectories have more timepoints

CTRL_TRAJ(:, all(isnan(CTRL_TRAJ))) = [];   % delete empty columns

multiple_trajectories = CTRL_TRAJ';

num_trajectories = size(multiple_trajectories, 1);

% Define the x-axis coordinates for imagesc
% This maps your data's columns (e.g., 121 points) to the desired time range (1 to 360 min)
num_time_points_in_data = size(multiple_trajectories, 2);
x_coordinates = linspace(0, 360, num_time_points_in_data); % Generates evenly spaced points from 1 to 360

% Set predefined limits for the color scale (intensity)
predefined_min_intensity = 100;
predefined_max_intensity = 600;

% Create the figure and display the image using the custom x-coordinates
figure;
imagesc(x_coordinates, 1:num_trajectories, multiple_trajectories);

% figure;
% imagesc(multiple_trajectories); % imagesc automatically scales colors based on data range
% --- Apply the Predefined Color Scale ---
clim([predefined_min_intensity predefined_max_intensity]);
colormap('gray'); % Choose a suitable colormap (e.g., 'hot', 'jet', 'parula', 'gray')
colorbar; % Show the color scale
% title('Parental ctrl', 'FontSize', 30, 'FontWeight', 'bold');
xlabel('Time (minutes)','FontSize', 20, 'FontWeight', 'bold');
ylabel('Trajectory Index', 'FontSize', 20, 'FontWeight', 'bold');
% Set x-axis ticks at the desired intervals (1, 4, 7, ..., 360)
% These ticks will be placed on the 'x_coordinates' scale
xticks(0:60:360);
xlim([0 360]);

ax.XAxis.FontSize = 16; 
ax.XAxis.FontName = 'Arial';

ax.YAxis.FontSize = 16; 
ax.YAxis.FontName = 'Arial';
% Ensure x-axis limits span the full desired range (1 to 360)


% % Export as high-quality TIFF (300 DPI)
% exportgraphics(gcf, 'Klf4 ctrl trajectories.tiff', 'Resolution', 600);








