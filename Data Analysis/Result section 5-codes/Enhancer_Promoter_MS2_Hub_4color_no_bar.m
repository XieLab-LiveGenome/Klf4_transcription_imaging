% =========================================================
%  4-panel distance + MS2 intensity plot  (no error bars)
%  data_matrix : N×4  [MS2 intensity | E-C | E-P | P-C]
%  x           : N×1  time vector (seconds)
% =========================================================

x = [-10:1:12]';

data_matrix = [5356.135802	995.6345484	228.770435	1069.289021
4880.45679	577.6713864	407.5684144	569.9042131
6229.382716	937.9033896	366.3879017	959.6636296
3776.432099	821.229528	455.4496842	467.7229548
3838.315789	1075.257319	372.4096813	745.9282439
2599.173333	1265.49553	223.1423814	1205.854141
4466.617284	933.7898057	363.9713194	1189.978739
4710.641975	826.7192318	485.7850957	848.3383831
4362.7875	711.1520852	406.0201322	405.7362399
6185.469136	458.8982634	709.1913802	509.4381905
10229.85185	488.2542073	815.8440064	345.5465844
14484.54321	549.8829665	683.8969072	168.3107584
9275.098765	588.853271	572.8418184	76.32043753
26569.03704	452.1243838	397.1481749	309.5250113
20341.40741	506.2471346	525.4420286	73.99880258
7635.506173	581.6971299	616.3063591	359.933
8906.111111	961.5665061	830.6542367	147.8839111
13074.95062	748.1036511	683.0224504	84.05697987
8568.82716	879.6432409	828.0180807	276.172352
4687.391892	557.7597566	502.1551757	431.5353371
6167.185185	606.4677842	642.0540926	563.5504469
6285.54321	685.516234	581.0657416	950.4004208
4288.090909	718.8800176	695.665001	768.9265927];

% --- Unpack columns ---
ms2 = data_matrix(:,1);
ec  = data_matrix(:,2);
ep  = data_matrix(:,3);
pc  = data_matrix(:,4);

% --- Colors ---
col_ms2 = [0.12 0.47 0.71];   % blue
col_ec  = [0.07 0.45 0.13];   % deep green
col_ep  = [0.85 0.65 0.00];   % gold
col_pc  = [0.80 0.07 0.07];   % red

lw = 5;    % line width
fs = 18;   % font size

dist_labels = {'E–P distance (nm)', 'E–C distance (nm)', 'P–C distance (nm)'};
ydata       = {ep,  ec,  pc };
dist_cols   = {col_ep, col_ec, col_pc};

% --- Figure ---
fig = figure;
set(fig, 'Color', 'w');
tl = tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Panels 1–3 (individual distances) ---
for k = 1:3
    nexttile;

    yyaxis left;
    hold on;
    plot(x, ydata{k}, 'Color', dist_cols{k}, 'LineWidth', lw);
    ylabel(dist_labels{k}, 'Color', dist_cols{k});
    ylim([0 1500]);

    yyaxis right;
    hold on;
    plot(x, ms2, 'Color', col_ms2, 'LineWidth', lw);
    ylabel('MS2 intensity', 'Color', col_ms2);
    ylim([0 30000]);

    xline(0, '--k', 'LineWidth', 1.5, 'Alpha', 0.4);

    ax = gca;
    ax.Color          = 'w';
    ax.XColor         = 'k';
    ax.YAxis(1).Color = dist_cols{k};
    ax.YAxis(2).Color = col_ms2;
    ax.TickDir        = 'out';
    ax.LineWidth      = 2;
    ax.XTickLabel     = {};
    set(ax, 'FontSize', fs, 'FontWeight', 'bold');
    xlim([x(1) x(end)]);
    hold off;
end

% --- Panel 4 (all three distances overlaid) ---
nexttile;

yyaxis left;
hold on;
for k = 1:3
    plot(x, ydata{k}, 'Color', dist_cols{k}, 'LineWidth', lw, ...
        'DisplayName', strrep(dist_labels{k}, ' (nm)', ''));
end
ylabel('Distance (nm)', 'Color', 'k');
ylim([0 1500]);

yyaxis right;
hold on;
plot(x, ms2, 'Color', col_ms2, 'LineWidth', lw, 'DisplayName', 'MS2 intensity');
ylabel('MS2 intensity', 'Color', col_ms2);
ylim([0 30000]);

xline(0, '--k', 'LineWidth', 1.5, 'Alpha', 0.4);
xlabel('Time (seconds)');
xlim([x(1) x(end)]);

ax4 = gca;
ax4.Color          = 'w';
ax4.XColor         = 'k';
ax4.YAxis(1).Color = 'k';
ax4.YAxis(2).Color = col_ms2;
ax4.TickDir        = 'out';
ax4.LineWidth      = 2;
set(ax4, 'FontSize', fs, 'FontWeight', 'bold');

hold off;