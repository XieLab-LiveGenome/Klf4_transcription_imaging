
clear; clc;

%% ---------------- USER SETTINGS -------------------------------------
xlsxFile = '/Users/janaa/Desktop/MS2 transcription/figures ai format/Latest version/TTseq vs Bursting.xlsx';   % <-- put your workbook path here

% Sheet names in the workbook. Edit these to match yours exactly.
sheetDur   = 'Burst duration';
sheetFreq  = 'Burst frequency';
sheetAmp   = 'Burst amplitude';
sheetTTseq = 'TTsequencing';

% Range for the 9 factors. B2:C10 = Mean_norm, SEM_norm for 9 rows.
dataRange = 'c2:d10';
nameRange = 'b2:b10';   % factor names for a sanity check

%% ---------------- Expected factor order & colors --------------------
factors = {'BRD4','MED14','MED12','RAD21','CTCF','SOX2','TAF2','P300','BRG1'};

% ---- Per-panel axis limits (tune freely) ----------------------------
xlimDur   = [0.0 1.8];   % Norm. Burst duration
xlimFreq  = [0.0 1.4];   % Norm. Burst frequency
xlimAmp   = [0.0 1.6];   % Norm. Burst amplitude
ylimAll   = [0.0 1.4];   % shared y-axis (Norm. Klf4 expression)


% RGB (0-1). BRD4 = red. Tweak any row to fine-tune shades.
colors = [ ...
    0.84 0.15 0.16 ;   % BRD4  - red
    0.62 0.71 0.22 ;   % MED14 - olive / yellow-green
    0.17 0.63 0.35 ;   % MED12 - green
    0.83 0.24 0.55 ;   % SOX2  - magenta
    0.20 0.37 0.73 ;   % RAD21 - blue
    0.45 0.76 0.96 ;   % CTCF  - light blue
    0.92 0.68 0.85 ;   % TAF2  - light pink
    0.47 0.47 0.28 ;   % P300  - dark olive
    0.90 0.55 0.22 ];  % BRG1  - orange

%% ---------------- Read the 4 sheets ---------------------------------
[dur_mean,   dur_sem]   = readSheet(xlsxFile, sheetDur,   dataRange, nameRange, factors);
[freq_mean,  freq_sem]  = readSheet(xlsxFile, sheetFreq,  dataRange, nameRange, factors);
[amp_mean,   amp_sem]   = readSheet(xlsxFile, sheetAmp,   dataRange, nameRange, factors);
[ttseq_mean, ttseq_sem] = readSheet(xlsxFile, sheetTTseq, dataRange, nameRange, factors);

%% ---------------- Build the figure ----------------------------------
fig = figure('Color','w','Units','pixels','Position',[100 100 1380 440]);
tl  = tiledlayout(fig,1,3,'Padding','compact','TileSpacing','compact');

ax1 = nexttile;
plotCorrPanel(ax1, dur_mean,  dur_sem,  ttseq_mean, ttseq_sem, colors, 'Norm. Burst duration', xlimDur, ylimAll);

ax2 = nexttile;
plotCorrPanel(ax2, freq_mean, freq_sem, ttseq_mean, ttseq_sem, colors, 'Norm. Burst frequency', xlimFreq, ylimAll);

ax3 = nexttile;
plotCorrPanel(ax3, amp_mean, amp_sem, ttseq_mean, ttseq_sem, colors, 'Norm. Burst amplitude', xlimAmp, ylimAll);

% ---- Build the legend from CLEAN white-filled circle handles --------
% (errorbar handles carry the error-bar line into the legend, which we
% don't want. Draw invisible plot markers off-axes purely for legend.)
hLeg = gobjects(numel(factors),1);
hold(ax3,'on');
for i = 1:numel(factors)
    hLeg(i) = plot(ax3, NaN, NaN, 'o', ...
        'MarkerEdgeColor', colors(i,:), 'MarkerFaceColor', 'w', ...
        'MarkerSize', 9, 'LineWidth', 1.7, 'LineStyle', 'none');
end
lgd = legend(ax3, hLeg, factors, 'Location','eastoutside','Box','off','FontSize',11);
lgd.ItemTokenSize = [14 14];

%% ---------------- Export --------------------------------------------
exportgraphics(fig,'Klf4_burst_correlation.png','Resolution',300);
exportgraphics(fig,'Klf4_burst_correlation.pdf','ContentType','vector');

%% ===================================================================
function [mu, sem] = readSheet(file, sheet, dataRange, nameRange, expectedOrder)
% Read Mean_norm / SEM_norm from one sheet and verify factor order.
    M   = readmatrix(file, 'Sheet', sheet, 'Range', dataRange);
    mu  = M(:,1);
    sem = M(:,2);

    % Sanity check: names in col A should match the expected order.
    try
        names = readcell(file, 'Sheet', sheet, 'Range', nameRange);
        names = string(names(:));
        expected = string(expectedOrder(:));
        if numel(names) == numel(expected) && ~all(strcmpi(strtrim(names), strtrim(expected)))
            warning('Factor order in sheet "%s" does not match expected order.\nFound:    %s\nExpected: %s', ...
                sheet, strjoin(names, ', '), strjoin(expected, ', '));
        end
    catch
        % Silently skip the check if col A is missing or unreadable.
    end
end

function plotCorrPanel(ax, x, xsem, y, ysem, colors, xlab, xl, yl)
% One correlation panel with per-point white-filled colored markers +
% x/y error bars, dotted least-squares trend line, and Pearson r annotation.
    x = x(:); xsem = xsem(:); y = y(:); ysem = ysem(:);
    hold(ax,'on');
    n = numel(x);

    % --- least-squares trend line (drawn first so markers sit on top) ---
    p   = polyfit(x, y, 1);
    pad = 0.05*(max(x)-min(x));
    xf  = linspace(min(x)-pad, max(x)+pad, 50);
    plot(ax, xf, polyval(p, xf), ':', 'Color', [0.35 0.45 0.85], 'LineWidth', 1.5);

    % --- per-point error bars + white-filled colored circles ---
    for i = 1:n
        errorbar(ax, x(i), y(i), ysem(i), ysem(i), xsem(i), xsem(i), 'o', ...
            'Color', colors(i,:), 'MarkerEdgeColor', colors(i,:), ...
            'MarkerFaceColor', 'w', 'MarkerSize', 10, ...
            'LineWidth', 1.7, 'CapSize', 3);
    end

    % Pearson correlation (base MATLAB; no toolbox)
    [R, P] = corrcoef(x, y);
    r = R(1,2);  pval = P(1,2);

    % --- cosmetics ---
    xlim(ax, xl); ylim(ax, yl);
    axis(ax,'square'); box(ax,'off');
    set(ax,'FontSize',12,'LineWidth',1.2,'TickDir','out');
    xlabel(ax, xlab, 'FontWeight','bold','FontSize',13);
    ylabel(ax, 'Norm. Klf4 expression', 'FontWeight','bold','FontSize',13);

    % r/p annotation positioned relative to the panel's own limits
    xt = xl(1) + 0.05*(xl(2)-xl(1));
    yt = yl(1) + 0.95*(yl(2)-yl(1));
    text(ax, xt, yt, sprintf('r = %.3f\np = %.2g', r, pval), ...
        'FontSize', 12, 'FontWeight','bold', 'VerticalAlignment','top');
end