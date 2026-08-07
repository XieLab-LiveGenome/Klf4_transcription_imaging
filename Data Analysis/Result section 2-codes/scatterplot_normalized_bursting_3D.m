clear 
clc

colors = [
    1.0, 0.0, 0.0;      
    0.0, 0.5, 0.0;      
    0.2, 1.0, 0.0;      
    1.0, 0.0, 1; 
    1.0, 0.5, 0.8;    
    0.0, 1.0, 1.0; 
    0.0, 0.5, 1.0;  
    0.75 0.6, 0.1;     
    1.0, 0.5, 0.0;      
    0.8 0.8 0.8
];


data = [0.5616	0.4275	0.4939
1.0278	0.327	0.7487
1.3472	0.6015	0.4814
1.0575	0.3467	0.9483
1.2106	1.1033	1.3506
0.7456	0.5341	0.5546
0.7629	0.8133	1.1751
0.7305	0.6616	0.3988
1.2268	0.4839	0.6226
0.4751	0.5804	0.5808];

x = data(:, 1);
y = data(:, 2);
z = data(:, 3);

regulator_names = {'BRD4', 'MED14', 'MED12', 'RAD21', ...
                   'CTCF', 'SOX2', 'TAF2', ...
                   'P300', 'BRG1', 'del OCT4/SOX2'};   % del OCT4/SOX2 - same as SOX2;

figure('Color','w','Position',[100 100 1200 900]);
hold on;

% --- axis limits ---
xl = [0.2 max(x)*1.05];
yl = [0.2 max(y)*1.05];
zl = [0.2 max(z)*1.05];

% --- light shaded floor and back walls ---
[Fx, Fy] = meshgrid(xl, yl);
surf(Fx, Fy, zl(1)*ones(size(Fx)), ...
    'FaceColor', [0.95 0.95 0.95], 'FaceAlpha', 0.5, ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');

[Wx, Wz] = meshgrid(xl, zl);
surf(Wx, yl(2)*ones(size(Wx)), Wz, ...
    'FaceColor', [0.97 0.97 0.97], 'FaceAlpha', 0.4, ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');

[Wy, Wz2] = meshgrid(yl, zl);
surf(xl(1)*ones(size(Wy)), Wy, Wz2, ...
    'FaceColor', [0.97 0.97 0.97], 'FaceAlpha', 0.4, ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');

% --- depth-scaled spheres (subtle) ---
az = 37.5; el = 30;
cam_dist = 5;
cam_pos = cam_dist * [cosd(el)*cosd(az), cosd(el)*sind(az), sind(el)];
dists = sqrt((x - cam_pos(1)).^2 + (y - cam_pos(2)).^2 + (z - cam_pos(3)).^2);
min_d = min(dists);
max_d = max(dists);

base_radius = 0.04 * max([max(x); max(y); max(z)]);
min_scale = 0.75;   % was 0.55 — subtler now
max_scale = 1.05;   % was 1.15

for i = 1:10
    t = (dists(i) - min_d) / (max_d - min_d + eps);
    scale = max_scale - t * (max_scale - min_scale);
    radius = base_radius * scale;

    [sx, sy, sz] = sphere(40);

    surf(sx*radius + x(i), sy*radius + y(i), sz*radius + z(i), ...
        'FaceColor', colors(i,:), 'EdgeColor', 'none', ...
        'FaceLighting', 'gouraud', ...
        'AmbientStrength', 0.45, ...    % raised: more base color visible
        'DiffuseStrength', 0.75, ...
        'SpecularStrength', 0.7, ...    % toned down: less white washout
        'SpecularExponent', 20, ...
        'BackFaceLighting', 'lit', ...
        'DisplayName', regulator_names{i});
end
hold off;

% --- lighting ---
light('Position', [1 1 1.5], 'Style', 'infinite');
light('Position', [-1 -0.5 0.5], 'Style', 'infinite', 'Color', [0.4 0.4 0.5]);
camlight('headlight');

grid off;
xlabel('Burst duration', 'FontWeight', 'bold', 'FontSize', 30);
ylabel('Burst frequency', 'FontWeight', 'bold', 'FontSize', 30);
zlabel('Burst amplitude', 'FontWeight', 'bold', 'FontSize', 30);

xlim(xl);
ylim(yl);
zlim(zl);

view(37.5, 30);
camproj('perspective');

lgd = legend('Location', 'bestoutside', 'FontWeight', 'bold', 'TextColor', 'K');
lgd.Color = 'W';
lgd.EdgeColor = 'K';

ax = gca;
ax.LineWidth = 2;
ax.FontWeight = 'bold';
ax.FontSize = 20;
ax.Color = 'W';
ax.XColor = 'K';
ax.YColor = 'K';
ax.ZColor = 'K';
ax.BoxStyle = 'back';
box on;
set(gcf, 'Color', 'W');

% --- fit figure to PDF page ---
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [14 10]);
set(gcf, 'PaperPosition', [0 0 14 10]);

rotate3d on;