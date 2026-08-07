%% =====================================================================
%  E-P distance vs BRD4-hub distance: per-condition summary & correlation
%  Conditions: Control | R21 depl | del OCT4/SOX2 (G67) | G67+R21 | M14 dTAG
%
%  For each condition computes mean/SEM and median/[25-75 pctile] of:
%     - Enhancer-Promoter (E-P) distance        (column G)
%     - Enhancer-BRD4 hub distance              (column E)
%     - Promoter-BRD4 hub distance              (column K)
%  Then makes two scatter versions (mean+/-SEM and median+IQR), each with
%  5 points (one per condition), E-P distance on the x-axis, and reports
%  the Pearson correlation across conditions.
% =====================================================================
clear; clc;

%% ---------------- Inputs ----------------
files = { ...
    '/Users/janaa/Desktop/MS2 transcription/4 color compiled/4 COLOR_6-29-2026.xlsx', ...
    '/Users/janaa/Desktop/MS2 transcription/4 color compiled/4 COLOR g67.xlsx', ...
    '/Users/janaa/Desktop/MS2 transcription/4 color compiled/4 COLOR R21 depl.xlsx', ...
    '/Users/janaa/Desktop/MS2 transcription/4 color compiled/4 COLOR g67 +R21 dep.xlsx', ...
    '/Users/janaa/Desktop/MS2 transcription/4 color compiled/4 COLOR m14 dTAG 4-23-26.xlsx'};

conditions = {'Control','del OCT4/SOX2','RAD21(-)','del OCT4/SOX2 + RAD21(-)','MED14(-)'};

sheet     = 'ms intensity vs distance';
unitLabel = 'nm';          % <-- change to 'pixels' or '\mum' if needed
nCond     = numel(files);

thr   = 0.2*65536;   % MS2 intensity threshold (kept for downstream ON/OFF work)
thr_d = 1000;        % distance threshold (not used in this script)

%% ---------------- Load data ----------------
EP  = cell(1,nCond);   % enhancer-promoter distance      (col G)
dE  = cell(1,nCond);   % enhancer - BRD4 hub distance     (col E)
dP  = cell(1,nCond);   % promoter - BRD4 hub distance     (col K)
MS2 = cell(1,nCond);   % MS2 intensity                    (col D)

for i = 1:nCond
    dE{i}  = readmatrix(files{i},'Sheet',sheet,'Range','E2:E300');
    dP{i}  = readmatrix(files{i},'Sheet',sheet,'Range','K2:K300');
    EP{i}  = readmatrix(files{i},'Sheet',sheet,'Range','G2:G300');
    MS2{i} = readmatrix(files{i},'Sheet',sheet,'Range','D2:D300');

    % NaN removal per column (matches original behaviour)
    dE{i}  = dE{i}(~isnan(dE{i}));
    dP{i}  = dP{i}(~isnan(dP{i}));
    EP{i}  = EP{i}(~isnan(EP{i}));
    MS2{i} = MS2{i}(~isnan(MS2{i}));
end

%% ---------------- Per-condition summary statistics ----------------
sem = @(x) std(x)/sqrt(numel(x));

S = struct();
for i = 1:nCond
    S(i).cond = conditions{i};

    % E-P distance
    S(i).EP_n   = numel(EP{i});
    S(i).EP_mean= mean(EP{i});        S(i).EP_sem = sem(EP{i});
    S(i).EP_med = median(EP{i});
    S(i).EP_p25 = prctile(EP{i},25);  S(i).EP_p75 = prctile(EP{i},75);

    % Enhancer - BRD4 hub
    S(i).dE_n   = numel(dE{i});
    S(i).dE_mean= mean(dE{i});        S(i).dE_sem = sem(dE{i});
    S(i).dE_med = median(dE{i});
    S(i).dE_p25 = prctile(dE{i},25);  S(i).dE_p75 = prctile(dE{i},75);

    % Promoter - BRD4 hub
    S(i).dP_n   = numel(dP{i});
    S(i).dP_mean= mean(dP{i});        S(i).dP_sem = sem(dP{i});
    S(i).dP_med = median(dP{i});
    S(i).dP_p25 = prctile(dP{i},25);  S(i).dP_p75 = prctile(dP{i},75);
end

%% ---------------- Summary tables ----------------
SummaryTable = table( conditions(:), ...
    [S.EP_n]', [S.EP_mean]', [S.EP_sem]', [S.EP_med]', [S.EP_p25]', [S.EP_p75]', ...
    [S.dE_n]', [S.dE_mean]', [S.dE_sem]', [S.dE_med]', [S.dE_p25]', [S.dE_p75]', ...
    [S.dP_n]', [S.dP_mean]', [S.dP_sem]', [S.dP_med]', [S.dP_p25]', [S.dP_p75]', ...
    'VariableNames', {'Condition', ...
    'N_EP','EP_mean','EP_SEM','EP_median','EP_p25','EP_p75', ...
    'N_Ehub','Ehub_mean','Ehub_SEM','Ehub_median','Ehub_p25','Ehub_p75', ...
    'N_Phub','Phub_mean','Phub_SEM','Phub_median','Phub_p25','Phub_p75'});

disp('=================== SUMMARY STATISTICS ===================');
disp(SummaryTable);
% writetable(SummaryTable,'EP_hub_summary.xlsx');   % <-- uncomment to export

%% ---------------- Assemble plotting vectors ----------------
% Mean +/- SEM
x_mean  = [S.EP_mean]';   x_sem  = [S.EP_sem]';
yE_mean = [S.dE_mean]';   yE_sem = [S.dE_sem]';
yP_mean = [S.dP_mean]';   yP_sem = [S.dP_sem]';

% Median + [25-75] (asymmetric error lengths)
x_med  = [S.EP_med]';   x_lo  = x_med  - [S.EP_p25]';   x_hi  = [S.EP_p75]' - x_med;
yE_med = [S.dE_med]';   yE_lo = yE_med - [S.dE_p25]';   yE_hi = [S.dE_p75]' - yE_med;
yP_med = [S.dP_med]';   yP_lo = yP_med - [S.dP_p25]';   yP_hi = [S.dP_p75]' - yP_med;

%% ---------------- Pearson correlation across the 5 conditions ----------------
[rE_mean, pE_mean] = corr(x_mean, yE_mean);   % E-P vs E-hub  (means)
[rP_mean, pP_mean] = corr(x_mean, yP_mean);   % E-P vs P-hub  (means)
[rE_med,  pE_med ] = corr(x_med,  yE_med );   % E-P vs E-hub  (medians)
[rP_med,  pP_med ] = corr(x_med,  yP_med );   % E-P vs P-hub  (medians)

fprintf('\n================= PEARSON CORRELATION (n = %d conditions) =================\n', nCond);
fprintf('  Mean-based:   E-P vs E-hub:  r = %+.3f  (p = %.3f)\n', rE_mean, pE_mean);
fprintf('                E-P vs P-hub:  r = %+.3f  (p = %.3f)\n', rP_mean, pP_mean);
fprintf('  Median-based: E-P vs E-hub:  r = %+.3f  (p = %.3f)\n', rE_med,  pE_med);
fprintf('                E-P vs P-hub:  r = %+.3f  (p = %.3f)\n', rP_med,  pP_med);
fprintf('  (n = 5 -> low power; interpret p-values as descriptive.)\n');

%% ---------------- Plot ----------------
colors = [0.45 0.45 0.45;   % Control              – gray
          0.93 0.7 0.12;   % del OCT4/SOX2        – golden yellow
          0.90 0.2 0.2;   % R21 depl             – burnt orange
          0.50 0.10 0.50;   % G67 + R21            – purple
          0.45 0.55 0.00];  % M14 dTAG             – olive green

%% ---------------- Figure 1: Mean +/- SEM ----------------
% figure('Color','w','Position',[80 120 1100 470]);
% 
% subplot(1,2,1);
% drawPanel(x_mean,yE_mean,x_sem,x_sem,yE_sem,yE_sem,conditions,colors, ...
%     sprintf('E\\itP\\rm distance (%s)',unitLabel), ...
%     sprintf('Enhancer\\itBRD4\\rm hub distance (%s)',unitLabel), ...
%     sprintf('Mean \\pm SEM   |   Pearson r = %.2f, p = %.3f', rE_mean, pE_mean));
% 
% subplot(1,2,2);
% drawPanel(x_mean,yP_mean,x_sem,x_sem,yP_sem,yP_sem,conditions,colors, ...
%     sprintf('E\\itP\\rm distance (%s)',unitLabel), ...
%     sprintf('Promoter\\itBRD4\\rm hub distance (%s)',unitLabel), ...
%     sprintf('Mean \\pm SEM   |   Pearson r = %.2f, p = %.3f', rP_mean, pP_mean));
% 
% sgtitle('E\itP\rm distance vs BRD4-hub distance  (mean \pm SEM)');
% 
% %% ---------------- Figure 2: Median + [25-75 percentile] ----------------
figure('Color','w','Position',[120 160 1100 470]);

subplot(1,2,1);
drawPanel(x_med,yE_med,x_lo,x_hi,yE_lo,yE_hi,conditions,colors, ...
    sprintf('E\\itP\\rm distance (%s)',unitLabel), ...
    sprintf('Enhancer\\itBRD4\\rm hub distance (%s)',unitLabel), ...
    sprintf('Median + IQR   |   Pearson r = %.2f, p = %.3f', rE_med, pE_med));
ylim([0 2500]);

subplot(1,2,2);
drawPanel(x_med,yP_med,x_lo,x_hi,yP_lo,yP_hi,conditions,colors, ...
    sprintf('E\\itP\\rm distance (%s)',unitLabel), ...
    sprintf('Promoter\\itBRD4\\rm hub distance (%s)',unitLabel), ...
    sprintf('Median + IQR   |   Pearson r = %.2f, p = %.3f', rP_med, pP_med));
ylim([0 2500]);
sgtitle('E\itP\rm distance vs BRD4-hub distance  (median + 25\rm-\rm75 percentile)');

%% ---------------- Figure 3: E-P distance box plot ----------------
% Combine all EP data with integer group index
allEP  = vertcat(EP{:});
grpIdx = [];
for i = 1:nCond
    grpIdx = [grpIdx; i*ones(numel(EP{i}),1)]; %#ok<AGROW>
end

figure('Color','w','Position',[200 200 550 520]);
bp = boxplot(allEP, grpIdx, 'Widths',0.55, 'Symbol','');  % no outlier markers
hold on;

% --- Save median line coordinates before patches cover them ---
hMed = findobj(gca,'Tag','Median');
medX = cell(length(hMed),1);  medY = cell(length(hMed),1);
for j = 1:length(hMed)
    medX{j} = get(hMed(j),'XData');
    medY{j} = get(hMed(j),'YData');
end

% --- Style boxes: white fill, BLACK edges ---
hBox = findobj(gca,'Tag','Box');
for j = 1:length(hBox)
    xd = get(hBox(j),'XData');  yd = get(hBox(j),'YData');
    patch(xd, yd, 'w', 'EdgeColor','k', 'LineWidth',2.5, 'FaceColor','w');
end

% --- Redraw median lines ON TOP of patches, colored per condition ---
for j = 1:length(hMed)
    idx = nCond - j + 1;
    plot(medX{j}, medY{j}, '-', 'Color',colors(idx,:), 'LineWidth',5);
end

% --- Whiskers & caps: gray dashed ---
set(findobj(gca,'Tag','Whisker'),              'Color',[0.55 0.55 0.55],'LineStyle','--','LineWidth',1);
set(findobj(gca,'Tag','Upper Adjacent Value'), 'Color',[0.55 0.55 0.55],'LineWidth',1);
set(findobj(gca,'Tag','Lower Adjacent Value'), 'Color',[0.55 0.55 0.55],'LineWidth',1);

ylim([0 1500]);
% --- Labels ---
set(gca, 'XTickLabel', conditions, 'FontSize',11);
xtickangle(30);
ylabel(sprintf('Enhancer-Promoter 3D distance (%s)',unitLabel));
box on;

%% ---------------- Local function (must be at end of script) ----------------
function drawPanel(x,y,xneg,xpos,yneg,ypos,conds,cmap,xlab,ylab,ttl)
    hold on; box on;
    for k = 1:numel(x)
        % Error bars only (no marker)
        errorbar(x(k), y(k), yneg(k), ypos(k), xneg(k), xpos(k), ...
            'Color',cmap(k,:),'LineStyle','none','LineWidth',1.2, ...
            'CapSize',6,'Marker','none','HandleVisibility','off');
        % Marker only (thicker edge, white fill)
        plot(x(k), y(k), 'o', ...
            'Color',cmap(k,:),'MarkerFaceColor','w', ...
            'MarkerSize',13,'LineWidth',2.4,'DisplayName',conds{k});
    end
    pf = polyfit(x,y,1);                       % least-squares trend line
    xx = linspace(min(x),max(x),100);
    plot(xx, polyval(pf,xx),'k--','LineWidth',1,'HandleVisibility','off');
    xlabel(xlab); ylabel(ylab); title(ttl);
    legend('Location','best'); set(gca,'FontSize',11);
end