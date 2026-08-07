% =========================================================
%  4-panel distance + MS2 intensity plot
%  data_matrix : N×4  [MS2 intensity | E-C | E-P | P-C]
%  err_matrix  : N×4  [MS2 intensity | E-C | E-P | P-C]
%  x           : N×1  time vector (seconds)
% =========================================================

% --- Paste your data here --------------------------------
x = [-16:1:14]';          % replace with your time vector

data_matrix = [0.65351189	337.3531062	159.6951601	211.8831207
0.741958421	397.0107934	276.8082637	403.4935769
0.763969259	389.2390524	301.8217935	341.9117003
0.631475799	488.8121611	340.9126393	477.6355325
0.655411934	403.6901088	294.0906847	456.9304549
0.63718145	389.9989834	345.8011515	332.9187526
0.65295703	451.2098625	277.4356116	417.9518325
0.707270889	486.4367206	267.2583464	372.2869682
0.65181299	664.8807072	449.6478483	494.9298801
0.655232286	628.680787	393.6823934	360.6838159
0.642794797	582.2639876	280.7076809	502.8193733
0.647943891	612.3310753	387.9340372	563.2925773
0.640129852	722.9615634	540.9165361	545.9248338
0.566155541	638.6920685	506.818125	442.0492692
0.652849511	608.5828601	533.1403485	425.3731558
0.68436438	597.8231443	505.1387672	435.8728222
0.616547105	565.6381494	387.9503253	501.538625
0.279808895	706.5889321	430.432423	540.7620327
0.291368789	485.7802443	397.7273115	484.4055829
0.316323258	796.041617	475.9062892	794.3099307
0.296129342	795.4598602	477.8196147	669.7707582
0.278387784	719.5935489	485.5740362	584.3080623
0.288239265	611.4308065	363.7162252	709.5328408
0.308466681	742.7230886	329.7965233	711.4575761
0.33909837	834.0817554	350.9369816	734.6529263
0.389597682	829.4031556	315.9053009	774.6523493
0.39304127	998.7451101	431.2300302	1007.473495
0.416819598	843.2035973	502.0938599	686.7644446
0.32692075	1185.287162	482.0377672	1162.69595
0.381288484	911.5570296	513.8958585	842.6456686
0.382050599	949.7474984	361.8938898	979.3124265          % N rows × 4 cols
];

err_matrix = [0.125329371	50.11603383	46.31606724	43.91215281
0.080068746	96.92139682	103.5858584	107.6752947
0.070744838	91.17143917	65.4935598	68.35349506
0.085963112	13.64991442	56.98609559	35.20494557
0.061754554	47.20894029	46.62286927	88.92022384
0.064950405	63.64060397	55.68478892	73.98425183
0.043981399	68.7822117	34.31309007	63.27852969
0.061901047	66.7085836	50.3391005	51.59522875
0.064981675	69.7349078	68.35401997	69.69938378
0.058614944	82.52089639	50.71682808	65.85775924
0.047855139	66.01299012	39.95641886	70.32073515
0.057739234	87.97176418	69.10719631	82.05446985
0.040199634	89.63654672	85.11771911	93.02719465
0.046098211	80.81745043	81.42506551	71.97616397
0.040474716	92.76719673	90.39433177	56.81008307
0.045604098	66.31737031	80.40489478	64.09429017
0.043116845	67.14290319	47.98942877	70.80845265
0.034929961	133.0561414	71.56935383	97.59127886
0.03231423	61.42889051	91.85077298	71.89623813
0.033539118	143.5005916	76.21930661	133.6523976
0.039595896	181.0851364	96.66717913	148.7275479
0.042431501	118.8468123	104.5727628	113.6741242
0.05937967	114.0234414	67.78466497	212.7921205
0.057603672	107.5846611	45.2155377	170.0483745
0.062940076	209.2881108	63.22946138	205.8557908
0.064797045	141.4167013	48.0988328	153.906666
0.071712887	472.6029744	106.1517386	457.4757971
0.103076071	224.08568	94.96490527	190.3003407
0.075724296	365.1883992	78.95445106	361.2769584
0.08395765	176.0069273	46.73228882	93.87429327
0.107793478	231.1678544	79.43494055	318.2854091        % N rows × 4 cols
];
% ---------------------------------------------------------

% Unpack columns
ms2  = data_matrix(:,1);   ec  = data_matrix(:,2);   ep  = data_matrix(:,3);   pc  = data_matrix(:,4);
ems2 = err_matrix(:,1);    eec = err_matrix(:,2);     eep = err_matrix(:,3);     epc = err_matrix(:,4);

% Colors
col_ms2 = [0.12 0.47 0.71];   % blue      (MS2, secondary axis)
col_ec  = [0.07 0.45 0.13];   % deep green (E–C)
col_ep  = [0.85 0.65 0.00];   % gold       (E–P)
col_pc  = [0.80 0.07 0.07];   % red        (P–C)

lw  = 5;      % line width
fa  = 0.20;   % fill alpha
fs  = 18;     % font size

dist_labels = {'E–P distance (nm)', 'E–C distance (nm)', 'P–C distance (nm)'};
ydata       = {ep,  ec,  pc };
yerr        = {eep, eec, epc};
dist_cols   = {col_ep, col_ec, col_pc};

fig = figure;
set(fig, 'Color', 'w');

tl = tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% ---- helper: draw MS2 band + line on current yyaxis right ----
function plot_ms2(x, ms2, ems2, col_ms2, lw, fa)
    ub = ms2 + ems2;
    lb = ms2 - ems2;
    fill([x; flipud(x)], [lb; flipud(ub)], col_ms2, ...
        'FaceAlpha', fa, 'EdgeColor', 'none');
    plot(x, ms2, 'Color', col_ms2, 'LineWidth', lw);
end

% --- Panels 1–3 (individual) ----------------------------
for k = 1:3
    nexttile;

    y   = ydata{k};
    e   = yerr{k};
    col = dist_cols{k};
    ub  = y + e;
    lb  = y - e;

    % Left axis – distance
    yyaxis left;
    hold on;
    fill([x; flipud(x)], [lb; flipud(ub)], col, ...
        'FaceAlpha', fa, 'EdgeColor', 'none');
    plot(x, y, 'Color', col, 'LineWidth', lw);
    ylabel(dist_labels{k}, 'Color', col);

    % Right axis – MS2
    yyaxis right;
    hold on;
    plot_ms2(x, ms2, ems2, col_ms2, lw, fa);
    ylabel('MS2 intensity', 'Color', col_ms2);

    xline(0, '--k', 'LineWidth', 1.5, 'Alpha', 0.4);

    ax = gca;
    ax.Color          = 'w';
    ax.XColor         = 'k';
    ax.YAxis(1).Color = col;
    ax.YAxis(2).Color = col_ms2;
    ax.TickDir        = 'out';
    ax.LineWidth      = 2;
    ax.XTickLabel     = {};
    set(ax, 'FontSize', fs, 'FontWeight', 'bold');
    xlim([x(1) x(end)]);
    hold off;
end

% Left axis – distance
    yyaxis left;
    hold on;
    fill([x; flipud(x)], [lb; flipud(ub)], col, ...
        'FaceAlpha', fa, 'EdgeColor', 'none');
    plot(x, y, 'Color', col, 'LineWidth', lw);
    ylabel(dist_labels{k}, 'Color', col);
    % ylim([200 1000]);                          % <-- add this

    % Right axis – MS2
    yyaxis right;
    hold on;
    plot_ms2(x, ms2, ems2, col_ms2, lw, fa);
    ylabel('MS2 intensity', 'Color', col_ms2);
    ylim([0 1]);                               % <-- add this

% --- Panel 4 (merged) ------------------------------------
nexttile;

% Left axis – all three distances
yyaxis left;
hold on;

for k = 1:3
    y   = ydata{k};
    e   = yerr{k};
    col = dist_cols{k};
    ub  = y + e;
    lb  = y - e;

    fill([x; flipud(x)], [lb; flipud(ub)], col, ...
        'FaceAlpha', fa, 'EdgeColor', 'none');
    plot(x, y, 'Color', col, 'LineWidth', lw, ...
        'DisplayName', strrep(dist_labels{k}, ' (nm)', ''));
end

ylabel('3D Distance (nm)', 'Color', 'k');

% Right axis – MS2
yyaxis right;

% Left axis – all three distances
yyaxis left;
hold on;
% ... fill and plot calls ...
ylabel('Distance (nm)', 'Color', 'k');
ylim([0 1500]);                              % <-- add this

% Right axis – MS2
yyaxis right;
hold on;
% ... fill and plot calls ...
ylabel('MS2 intensity', 'Color', col_ms2);
ylim([0 1]);                                   % <-- add this
hold on;
fill([x; flipud(x)], [ms2-ems2; flipud(ms2+ems2)], col_ms2, ...
    'FaceAlpha', fa, 'EdgeColor', 'none');
plot(x, ms2, 'Color', col_ms2, 'LineWidth', lw, 'DisplayName', 'MS2 intensity');
ylabel('MS2 intensity', 'Color', col_ms2);

xline(0, '--k', 'LineWidth', 1.5, 'Alpha', 0.4);

xlabel('Time (seconds)');
xlim([0.1 1.0]);

ax4 = gca;
ax4.Color          = 'w';
ax4.XColor         = 'k';
ax4.YAxis(1).Color = 'k';
ax4.YAxis(2).Color = col_ms2;
ax4.TickDir        = 'out';
ax4.LineWidth      = 2;
set(ax4, 'FontSize', fs, 'FontWeight', 'bold');
xlim([x(1) x(end)]);

% lg = legend('Location', 'northwest', 'Box', 'off');
% lg.FontSize   = fs - 4;
% lg.FontWeight = 'bold';

hold off;