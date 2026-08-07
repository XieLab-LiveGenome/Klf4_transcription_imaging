clear
clc

% --- Data Loading ---
file_path = '/Users/janaa/Desktop/MS2 transcription/MS2 bursting compiled all_5-20-2026.xlsx';

WT_ON_cellavg  = readmatrix(file_path,'Sheet','WT ALL', 'Range','B3:B600');
WT_OFF_cellavg = readmatrix(file_path,'Sheet','WT ALL', 'Range','C3:C600');
WT_amp_cellavg = readmatrix(file_path,'Sheet','WT ALL', 'Range','F3:F600');
WT_ON_pooled   = readmatrix(file_path,'Sheet','WT ALL', 'Range','D3:D600');
WT_OFF_pooled  = readmatrix(file_path,'Sheet','WT ALL', 'Range','E3:E600');
WT_amp_pooled  = readmatrix(file_path,'Sheet','WT ALL', 'Range','H3:H600');

G67_ON_cellavg  = readmatrix(file_path,'Sheet','G67 ALL', 'Range','B3:B600');
G67_OFF_cellavg = readmatrix(file_path,'Sheet','G67 ALL', 'Range','C3:C600');
G67_amp_cellavg = readmatrix(file_path,'Sheet','G67 ALL', 'Range','F3:F600');
G67_ON_pooled   = readmatrix(file_path,'Sheet','G67 ALL', 'Range','D3:D600');
G67_OFF_pooled  = readmatrix(file_path,'Sheet','G67 ALL', 'Range','E3:E600');
G67_amp_pooled  = readmatrix(file_path,'Sheet','G67 ALL', 'Range','H3:H600');

%% ========================== FIGURE BUILD ==========================
clean = @(x) x(~isnan(x) & isfinite(x));

cWT  = [0.3 0.3 0.3];             % WT = black
cDel = [0.85 0.2 0.2];      % ΔOCT4/SOX2 = blue

ylabs = {'ON duration (min)', 'OFF duration (min)', 'Burst amplitude'};

% Cell-average arrays: {ON, OFF, AMP}
avg_WT  = {WT_ON_cellavg,  WT_OFF_cellavg,  WT_amp_cellavg};
avg_G67 = {G67_ON_cellavg, G67_OFF_cellavg, G67_amp_cellavg};

fig = figure('Units','inches','Position',[1 2 10 3.5],'Color','w');
tiledlayout(1, 3, 'TileSpacing','compact', 'Padding','compact');

for r = 1:3
    nexttile;
    a = clean(avg_WT{r});
    b = clean(avg_G67{r});
    m1 = mean(a);  sem1 = std(a)/sqrt(numel(a));
    m2 = mean(b);  sem2 = std(b)/sqrt(numel(b));

    % --- Filled bars ---
    bar(1, m1, 0.55, 'FaceColor', cWT,  'EdgeColor', cWT,  'LineWidth', 1.2); hold on;
    bar(2, m2, 0.55, 'FaceColor', cDel, 'EdgeColor', cDel, 'LineWidth', 1.2);

    % --- SEM error bars ---
    errorbar(1, m1, sem1, 'Color', 'k', 'LineStyle','none', 'LineWidth', 1.3, 'CapSize', 10);
    errorbar(2, m2, sem2, 'Color', 'k', 'LineStyle','none', 'LineWidth', 1.3, 'CapSize', 10);

    % --- t-test (Welch's) ---
    if numel(a) > 1 && numel(b) > 1
        [~, pT] = ttest2(a, b, 'Vartype','unequal');
    else
        pT = NaN;
    end

    % --- Significance bar ---
    yMax = max([m1+sem1, m2+sem2]);
    if ~isfinite(yMax) || yMax <= 0, yMax = 1; end
    yBar = yMax * 1.18;
    yTxt = yMax * 1.32;
    plot([1 2], [yBar yBar], 'k-', 'LineWidth', 1.1);
    plot([1 1], [yBar*0.97 yBar], 'k-', 'LineWidth', 1.1);
    plot([2 2], [yBar*0.97 yBar], 'k-', 'LineWidth', 1.1);
    text(1.5, yTxt, formatP(pT), 'HorizontalAlignment','center', ...
         'FontSize',11, 'Color',[0.25 0.25 0.25], 'FontWeight','bold');

    % --- Cosmetics ---
    xlim([0.4 2.6]);
    ylim([0, yTxt * 1.15]);
    xticks([1 2]);
    xticklabels({'WT','\DeltaOCT4/SOX2'});
    ylabel(ylabs{r}, 'FontWeight','bold', 'FontSize',12);
    set(gca,'FontSize',11,'LineWidth',1,'TickDir','out','Box','off');
    hold off;
end

%% ---------- Save ----------
exportgraphics(fig, 'WT_vs_delOCT4SOX2_bursting.pdf', 'ContentType','vector');
exportgraphics(fig, 'WT_vs_delOCT4SOX2_bursting.png', 'Resolution', 400);

%% ---------- Helper ----------
function s = formatP(p)
    if isnan(p)
        s = 'p=NaN';
    elseif p >= 1e-3
        s = sprintf('p=%.4f', p);
    else
        e = floor(log10(p));
        m = p / 10^e;
        s = sprintf('p=%.2f*10^{%d}', m, e);
    end
end