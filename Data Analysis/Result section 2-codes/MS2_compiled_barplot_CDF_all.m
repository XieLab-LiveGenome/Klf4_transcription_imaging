%% =====================================================================
%  MS2 transcriptional bursting — compiled barplot and CDF figure
%
%    Rows 1-3 : bar charts of per-cell averages (ON, OFF, amplitude)
%    Rows 4-6 : pooled-burst ksdensity CDFs (ON, OFF, amplitude)
%    Columns  : BRD4 | MED14 | MED12 | RAD21 | CTCF | SOX2 | TAF2 | P300i | BRG1i
%
% =====================================================================
clear; clc; close all;

%% ---------------------- 1. INPUT / FACTOR TABLE ----------------------
file_path = '/Users/janaa/Desktop/MS2 transcription/MS2 data sheets/MS2 bursting compiled all_7-30-2026.xlsx';
out_stem  = 'ms2_bursting_9factor_183x240';

% One entry per FIGURE column. Add/remove entries here and the panel grid,
% stats and summary tables all follow automatically.
cfg = struct( ...
    'name',       {'BRD4','MED14','MED12','RAD21','CTCF','SOX2','TAF2','P300i','BRG1i'}, ...
    'ctrlSheet',  {'BRD4 CTRL ALL','MED14 CTRL ALL','MED12 CTRL ALL','R21 CTRL ALL', ...
                   'CTCF CTRL ALL','SOX2 CTRL ALL','TAF2 CTRL ALL', ...
                   'unedited ctrl all','unedited ctrl all'}, ...
    'treatSheet', {'BRD4 DTAG ALL','MED14 DTAG ALL','MED12 DTAG ALL','R21 DTAG ALL', ...
                   'CTCF DTAG ALL','SOX2 DTAG ALL','TAF2 DTAG ALL', ...
                   'P300_inh ALL','BRG1_inh ALL'}, ...
    'ctrlLab',    {'Uninduced','Uninduced','Uninduced','Uninduced','Uninduced','Uninduced','Uninduced','WT','WT'}, ...
    'treatLab',   {'dTAG(+)','dTAG(+)','dTAG(+)','dTAG(+)','dTAG(+)','dTAG(+)','dTAG(+)','P300i','BRG1i'} );

% Extra factors: quantified and normalised, but NOT drawn in the figure.
cfgExtra = struct( ...
    'name',       {'del OCT4/SOX2'}, ...
    'ctrlSheet',  {'WT ALL'}, ...
    'treatSheet', {'G67 ALL'}, ...
    'ctrlLab',    {'WT'}, ...
    'treatLab',   {'del OCT4/SOX2'} );

nF      = numel(cfg);            % figure columns
factors = {cfg.name};

cfgAll     = [cfg, cfgExtra];    % everything that gets quantified
nA         = numel(cfgAll);
factorsAll = {cfgAll.name};

fprintf('Figure factors: %d  |  quantified factors: %d (extra: %s)\n', ...
        nF, nA, strjoin({cfgExtra.name}, ', '));

% Sheet layout: read B:H once per sheet, then slice.
%   B = per-cell ON   C = per-cell OFF   F = per-cell amplitude
%   D = pooled  ON    E = pooled  OFF    H = pooled  amplitude
dataRange = 'B2:H700';
iAvg  = [1 2 5];        % B, C, F  within B..H
iPool = [3 4 7];        % D, E, H  within B..H

params    = {'ON_duration','OFF_duration','Burst_amplitude'};
paramLab  = {'ON duration (min)','OFF duration (min)','Burst amplitude'};
nP        = numel(params);

%% ---------------------- 2. LAYOUT (all mm) ---------------------------
L.pageW    = 200;    L.pageH   = 220;
L.mLeft    = 12.5;   L.mRight  = 2.0;
L.mTop     = 5.2;    L.mBot    = 2.0;
L.colGap   = 5.2;    % 6.2 in the 7-column spec; tightened to fit 9 columns
L.rowGap   = 2.2;    % between rows within a block
L.blockGap = 3.0;    % between the bar block and the CDF block
L.labBar   = 9.6;    % label strip under each bar row  (rotated tick labels)
L.labCDF   = 9.4;    % label strip under each CDF row  (ticks + x title)

panelW = (L.pageW - L.mLeft - L.mRight - (nF-1)*L.colGap) / nF;
panelH = (L.pageH - L.mTop - L.mBot - nP*L.labBar - nP*L.labCDF ...
          - 4*L.rowGap - L.blockGap) / 6;

fprintf('Grid: %d cols  |  panel = %.2f x %.2f mm  |  col gap = %.2f mm\n', ...
        nF, panelW, panelH, L.colGap);
if panelW < 12
    warning('Panel width %.1f mm is very tight — reduce L.colGap or L.mLeft.', panelW);
end

% Column left edges and row bottom edges (mm, from page bottom-left)
xL = L.mLeft + (0:nF-1)*(panelW + L.colGap);
yB = zeros(1,6);
yCursor = L.pageH - L.mTop;
for r = 1:6
    yB(r) = yCursor - panelH;
    if r <= nP, strip = L.labBar; else, strip = L.labCDF; end
    if r == nP,      gap = L.blockGap;
    elseif r == 6,   gap = 0;
    else,            gap = L.rowGap;
    end
    yCursor = yCursor - panelH - strip - gap;
end
% Normalised position of panel (row r, column i)
posn = @(r,i) [xL(i)/L.pageW, yB(r)/L.pageH, panelW/L.pageW, panelH/L.pageH];

%% ---------------------- 3. STYLE -------------------------------------
S.font      = 'Arial';
S.fs        = 7;        % all text
S.cC        = [0 0 0];          % control  = black
S.cD        = [0.851 0 0];      % treated  = red
S.axGray    = [0.149 0.149 0.149];
S.lwAxis    = 0.5;      % spines & ticks
S.lwCDF     = 0.7;      % CDF curves
S.lwBarEdge = 0.4;      % bar edges
S.lwErr     = 0.6;      % error bars
S.capSize   = 1.5;
S.barWidth  = 0.6;
S.tickLen   = [0.035 0.035];
S.xTitleCol = ceil(nF/2);   % which column carries the x-axis title (0 = all)

% Axis settings per parameter (ON, OFF, AMP)
xlims   = {[0 240],      [0 480],      [0 1200]};
xticksv = {0:120:240,    0:240:480,    0:600:1200};   % 3 ticks: panels are narrow
bws     = [0.4, 0.4, 2];        % ksdensity CDF bandwidths

clean = @(x) x(~isnan(x) & isfinite(x));

%% ---------------------- 4. LOAD DATA ---------------------------------
sheets = string(sheetnames(file_path));
cache  = containers.Map('KeyType','char','ValueType','any');

for i = 1:nA
    for k = {'ctrlSheet','treatSheet'}
        s = cfgAll(i).(k{1});
        if ~isKey(cache, s)
            cache(s) = readBlock(file_path, s, dataRange, sheets, 7);
        end
    end
end

avg_C  = cell(nP,nA);  avg_D  = cell(nP,nA);
pool_C = cell(nP,nA);  pool_D = cell(nP,nA);
for i = 1:nA
    A = cache(cfgAll(i).ctrlSheet);
    B = cache(cfgAll(i).treatSheet);
    for r = 1:nP
        avg_C{r,i}  = A(:, iAvg(r));    avg_D{r,i}  = B(:, iAvg(r));
        pool_C{r,i} = A(:, iPool(r));   pool_D{r,i} = B(:, iPool(r));
    end
end

% Flag shared controls (P300i / BRG1i both use 'unedited ctrl all')
[uSheets, ~, ic] = unique({cfgAll.ctrlSheet});
for u = 1:numel(uSheets)
    hits = find(ic == u);
    if numel(hits) > 1
        fprintf('Note: %s share the control sheet "%s" — their black bars/curves are identical.\n', ...
                strjoin(factorsAll(hits), ', '), uSheets{u});
    end
end

%% ---------------------- 5. FIGURE (9 columns only) -------------------
fig = figure('Units','centimeters', ...
             'Position',[1 1 L.pageW/10 L.pageH/10], ...
             'Color','w', 'InvertHardcopy','off', ...
             'PaperUnits','centimeters', ...
             'PaperSize',[L.pageW L.pageH]/10, ...
             'PaperPosition',[0 0 L.pageW L.pageH]/10, ...
             'PaperPositionMode','manual');

% ============ ROWS 1-3 : per-cell averages, bars + Welch t-test ============
for r = 1:nP
    for i = 1:nF
        ax = axes('Parent',fig,'Units','normalized','Position',posn(r,i)); 
        hold(ax,'on');

        a = clean(avg_C{r,i});
        b = clean(avg_D{r,i});

        if isempty(a) || isempty(b)
            text(ax,0.5,0.5,'no data','Units','normalized', ...
                 'HorizontalAlignment','center','FontName',S.font,'FontSize',S.fs);
            styleAxis(ax,S); axis(ax,'off');
            if r == 1, title(ax,factors{i},'FontName',S.font,'FontSize',S.fs,'FontWeight','bold'); end
            continue
        end

        m1 = mean(a);  sem1 = std(a)/sqrt(numel(a));
        m2 = mean(b);  sem2 = std(b)/sqrt(numel(b));

        bar(ax, 1, m1, S.barWidth, 'FaceColor',S.cC, 'EdgeColor',S.cC, 'LineWidth',S.lwBarEdge);
        bar(ax, 2, m2, S.barWidth, 'FaceColor',S.cD, 'EdgeColor',S.cD, 'LineWidth',S.lwBarEdge);
        errorbar(ax, 1, m1, sem1, 'Color',S.axGray, 'LineStyle','none', ...
                 'LineWidth',S.lwErr, 'CapSize',S.capSize);
        errorbar(ax, 2, m2, sem2, 'Color',S.axGray, 'LineStyle','none', ...
                 'LineWidth',S.lwErr, 'CapSize',S.capSize);

        if numel(a) > 1 && numel(b) > 1
            [~, pT] = ttest2(a, b, 'Vartype','unequal');
        else
            pT = NaN;
        end

        % --- y-axis autoscaled to the data; top tick tracks the bars, ---
        % --- with a fixed band reserved above it for the p-value bracket ---
        dataTop = max([m1+sem1, m2+sem2]);
        if ~isfinite(dataTop) || dataTop <= 0, dataTop = 1; end

        yTick   = niceCeil(dataTop);     % highest y-tick sits just above the taller bar
        yLimTop = yTick * 1.40;          % top ~28% of the panel holds the bracket + p label

        yBar = dataTop + 0.08*yTick;     % bracket just clears the taller bar/error cap
        tk   = 0.035*yTick;
        plot(ax,[1 2],[yBar yBar],'-','Color',S.axGray,'LineWidth',S.lwAxis);
        plot(ax,[1 1],[yBar-tk yBar],'-','Color',S.axGray,'LineWidth',S.lwAxis);
        plot(ax,[2 2],[yBar-tk yBar],'-','Color',S.axGray,'LineWidth',S.lwAxis);
        text(ax, 1.5, yBar, formatP(pT), 'HorizontalAlignment','center', ...
             'VerticalAlignment','bottom', 'FontName',S.font, 'FontSize',S.fs, ...
             'FontWeight','bold', 'Color',S.axGray, 'Interpreter','tex');

        xlim(ax,[0.35 2.65]);
        ylim(ax,[0 yLimTop]);
        xticks(ax,[1 2]);
        xticklabels(ax,{cfg(i).ctrlLab, cfg(i).treatLab});
        xtickangle(ax,45);
        yticks(ax,[0 yTick/2 yTick]);
        if i == 1
            ylabel(ax, paramLab{r}, 'FontName',S.font,'FontSize',S.fs,'FontWeight','bold');
        end
        if r == 1
            title(ax, factors{i}, 'FontName',S.font,'FontSize',S.fs,'FontWeight','bold');
        end
        styleAxis(ax,S);
        hold(ax,'off');
    end
end

% ============ ROWS 4-6 : pooled ksdensity CDFs + KS test ============
for r = 1:nP
    for i = 1:nF
        ax = axes('Parent',fig,'Units','normalized','Position',posn(nP+r,i)); 
        hold(ax,'on');

        x1 = clean(pool_C{r,i});
        x2 = clean(pool_D{r,i});
        x_query = linspace(xlims{r}(1), xlims{r}(2), 2000);

        if numel(x1) > 2
            [f1, xi1] = ksdensity(x1,'function','cdf','Bandwidth',bws(r));
            f1q = min(max(interp1(xi1,f1,x_query,'linear','extrap'),0),1);
            plot(ax, x_query, f1q, 'LineWidth',S.lwCDF, 'Color',S.cC);
        end
        if numel(x2) > 2
            [f2, xi2] = ksdensity(x2,'function','cdf','Bandwidth',bws(r));
            f2q = min(max(interp1(xi2,f2,x_query,'linear','extrap'),0),1);
            plot(ax, x_query, f2q, 'LineWidth',S.lwCDF, 'Color',S.cD);
        end

        if numel(x1) > 1 && numel(x2) > 1
            [~, pKS] = kstest2(x1, x2);
        else
            pKS = NaN;
        end

        xlim(ax, xlims{r});  ylim(ax,[0 1.02]);
        xticks(ax, xticksv{r});  yticks(ax,[0 0.5 1]);
        xtickangle(ax,45);
        if i == 1
            ylabel(ax,'Cumulative freq.','FontName',S.font,'FontSize',S.fs,'FontWeight','bold');
        else
            yticklabels(ax,[]);      % identical 0/0.5/1 scale — label once per row
        end
        if S.xTitleCol == 0 || i == S.xTitleCol
            xlabel(ax, paramLab{r}, 'FontName',S.font,'FontSize',S.fs,'FontWeight','bold');
        end

        text(ax, 0.96, 0.10, formatP(pKS), 'Units','normalized', ...
             'HorizontalAlignment','right', 'FontName',S.font, 'FontSize',S.fs, ...
             'FontWeight','bold', 'Color',S.axGray, 'Interpreter','tex');

        styleAxis(ax,S);
        hold(ax,'off');
    end
end

%% ---------------------- 6. EXPORT ------------------------------------
% print (not exportgraphics) so the page stays exactly at the spec size
print(fig, [out_stem '.pdf'], '-dpdf',  '-painters');
print(fig, [out_stem '.png'], '-dpng',  '-r600');
fprintf('Figure written: %s.pdf / %s.png\n', out_stem, out_stem);

%% ---------------------- 7. SUMMARY TABLE -----------------------------
% Includes del OCT4/SOX2 even though it is not drawn in the figure.
ensMean_C = nan(nP,nA); ensMean_D = nan(nP,nA);
ensSEM_C  = nan(nP,nA); ensSEM_D  = nan(nP,nA);
poolMed_C = nan(nP,nA); poolMed_D = nan(nP,nA);
pTT       = nan(nP,nA); pKSm      = nan(nP,nA);
nCell_C   = zeros(nP,nA); nCell_D = zeros(nP,nA);
nBurst_C  = zeros(nP,nA); nBurst_D= zeros(nP,nA);

for r = 1:nP
    for i = 1:nA
        a  = clean(avg_C{r,i});   b  = clean(avg_D{r,i});
        x1 = clean(pool_C{r,i});  x2 = clean(pool_D{r,i});

        nCell_C(r,i)  = numel(a);   nCell_D(r,i)  = numel(b);
        nBurst_C(r,i) = numel(x1);  nBurst_D(r,i) = numel(x2);

        if ~isempty(a), ensMean_C(r,i) = mean(a); ensSEM_C(r,i) = std(a)/sqrt(numel(a)); end
        if ~isempty(b), ensMean_D(r,i) = mean(b); ensSEM_D(r,i) = std(b)/sqrt(numel(b)); end
        if ~isempty(x1), poolMed_C(r,i) = median(x1); end
        if ~isempty(x2), poolMed_D(r,i) = median(x2); end

        if numel(a) > 1 && numel(b) > 1,   [~, pTT(r,i)]  = ttest2(a,b,'Vartype','unequal'); end
        if numel(x1) > 1 && numel(x2) > 1, [~, pKSm(r,i)] = kstest2(x1,x2); end
    end
end

fid = fopen('bursting_summary_table.txt','w');
fprintf(fid, '=== FACTOR / SHEET MAP ===\n');
fprintf(fid, '(%s is quantified below but is NOT a column in the figure)\n', cfgExtra(1).name);
fprintf(fid, 'Factor\tControl sheet\tTreated sheet\tControl label\tTreated label\tIn figure\n');
for i = 1:nA
    if i <= nF, inFig = 'yes'; else, inFig = 'no'; end
    fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\n', factorsAll{i}, cfgAll(i).ctrlSheet, ...
            cfgAll(i).treatSheet, cfgAll(i).ctrlLab, cfgAll(i).treatLab, inFig);
end

fprintf(fid, '\n=== ENSEMBLE CELL-AVERAGE MEAN +/- SEM (n = cells) ===\n');
fprintf(fid, 'Parameter\tCondition');
fprintf(fid, '\t%s', factorsAll{:});  fprintf(fid, '\n');
for r = 1:nP
    fprintf(fid, '%s\tcontrol', params{r});
    for i = 1:nA
        fprintf(fid, '\t%.2f +/- %.2f (n=%d)', ensMean_C(r,i), ensSEM_C(r,i), nCell_C(r,i));
    end
    fprintf(fid, '\n%s\ttreated', params{r});
    for i = 1:nA
        fprintf(fid, '\t%.2f +/- %.2f (n=%d)', ensMean_D(r,i), ensSEM_D(r,i), nCell_D(r,i));
    end
    fprintf(fid, '\n%s\tWelch t-test p', params{r});
    for i = 1:nA, fprintf(fid, '\t%s', formatPplain(pTT(r,i))); end
    fprintf(fid, '\n');
end

fprintf(fid, '\n=== POOLED MEDIAN (n = bursts) ===\n');
fprintf(fid, 'Parameter\tCondition');
fprintf(fid, '\t%s', factorsAll{:});  fprintf(fid, '\n');
for r = 1:nP
    fprintf(fid, '%s\tcontrol', params{r});
    for i = 1:nA, fprintf(fid, '\t%.2f (n=%d)', poolMed_C(r,i), nBurst_C(r,i)); end
    fprintf(fid, '\n%s\ttreated', params{r});
    for i = 1:nA, fprintf(fid, '\t%.2f (n=%d)', poolMed_D(r,i), nBurst_D(r,i)); end
    fprintf(fid, '\n%s\tKS-test p', params{r});
    for i = 1:nA, fprintf(fid, '\t%s', formatPplain(pKSm(r,i))); end
    fprintf(fid, '\n');
end
fclose(fid);
fprintf('Summary table saved to: bursting_summary_table.txt (tab-delimited)\n');

%% ------- 8. CONTROL-NORMALISED PER-CELL TABLES (CSV) -----------------
% For every factor x parameter, take the MEAN of that factor's own
% control/uninduced per-cell values as the reference (ref).
%   invertNorm(r) == false  -> normalised value = per_cell_value / ref
%                              ("value_over_control_mean"; control mean = 1,
%                               treated = fold change).
%   invertNorm(r) == true   -> normalised value = ref / per_cell_value
%                              ("burst_frequency"; used for OFF duration,
%                               since burst frequency ~ 1/OFF). Each cell
%                               gets ref / (its own OFF), rescaled so the
%                               control column is centred at exactly 1 the
%                               same way ON and amplitude are; treated < 1
%                               means the cells burst less often.
% The maths lives in normPair() so sections 8 and 9 can never drift apart.
invertNorm = [false, true, false];   % ON, OFF, AMP
paramNorm  = {'ON_duration','burst_frequency','Burst_amplitude'};  % field name in outputs
export_percell_values = true;

Fac = {}; Par = {}; Cond = {}; Nn = []; Mn = []; Sdn = []; Semn = []; Ref = []; Nrm = {};
pFac = {}; pPar = {}; pCtrl = []; pTreat = []; pFold = []; pPct = [];
pSDpct = []; pSEMpct = []; pN_C = []; pN_D = []; pP = []; pNrm = {};
valKeys = {}; valCols = {};

for i = 1:nA
    for r = 1:nP
        a = clean(avg_C{r,i});
        b = clean(avg_D{r,i});

        [an, bn, ref, ok] = normPair(a, b, invertNorm(r));
        if ~ok
            warning('%s / %s: control mean unusable — normalisation skipped.', ...
                    factorsAll{i}, params{r});
            continue
        end
        if invertNorm(r)
            normLabel = 'burst_frequency';
        else
            normLabel = 'value_over_control_mean';
        end
        parName = paramNorm{r};

        % --- normalised summary rows (control, then treated) ---
        Fac(end+1)  = {factorsAll{i}};   Par(end+1) = {parName};   %#ok<*SAGROW>
        Cond(end+1) = {cfgAll(i).ctrlLab};
        Nn(end+1)   = numel(an);
        Mn(end+1)   = mean(an);
        Sdn(end+1)  = std(an);
        Semn(end+1) = std(an)/sqrt(numel(an));
        Ref(end+1)  = ref;
        Nrm(end+1)  = {normLabel};

        Fac(end+1)  = {factorsAll{i}};   Par(end+1) = {parName};
        Cond(end+1) = {cfgAll(i).treatLab};
        Nn(end+1)   = numel(bn);
        Mn(end+1)   = mean(bn);
        Sdn(end+1)  = std(bn);
        Semn(end+1) = std(bn)/sqrt(numel(bn));
        Ref(end+1)  = ref;
        Nrm(end+1)  = {normLabel};

        % --- percent change (treated vs its own control mean) ---
        % p on the same cleaned arrays the figure uses.
        if numel(a) > 1 && numel(b) > 1
            [~, pv] = ttest2(a, b, 'Vartype','unequal');
        else
            pv = NaN;
        end
        pFac(end+1)   = {factorsAll{i}};   pPar(end+1) = {parName};
        pCtrl(end+1)  = ref;               pTreat(end+1) = mean(b);
        pFold(end+1)  = mean(bn);
        pPct(end+1)   = (mean(bn) - 1)*100;
        pSDpct(end+1) = std(bn)*100;
        pSEMpct(end+1)= std(bn)/sqrt(numel(bn))*100;
        pN_C(end+1)   = numel(a);          pN_D(end+1) = numel(b);
        pP(end+1)     = pv;
        pNrm(end+1)   = {normLabel};

        if export_percell_values
            valKeys(end+1) = {sprintf('%s_%s_%s', factorsAll{i}, parName, cfgAll(i).ctrlLab)};
            valCols(end+1) = {an(:)};
            valKeys(end+1) = {sprintf('%s_%s_%s', factorsAll{i}, parName, cfgAll(i).treatLab)};
            valCols(end+1) = {bn(:)};
        end
    end
end

Tnorm = table(Fac(:), Par(:), Cond(:), Nn(:), Mn(:), Sdn(:), Semn(:), Ref(:), Nrm(:), ...
    'VariableNames', {'Factor','Parameter','Condition','n_cells', ...
                      'Mean_norm','SD_norm','SEM_norm','Control_mean_raw', ...
                      'Normalization'});
writetable(Tnorm, 'bursting_normalized_summary.csv');

Tpct = table(pFac(:), pPar(:), pCtrl(:), pTreat(:), pFold(:), pPct(:), ...
             pSDpct(:), pSEMpct(:), pN_C(:), pN_D(:), pP(:), pNrm(:), ...
    'VariableNames', {'Factor','Parameter','Control_mean_raw','Treated_mean_raw', ...
                      'Fold_change','Percent_change','SD_percent_points', ...
                      'SEM_percent_points','n_control_cells','n_treated_cells', ...
                      'p_Welch_ttest','Normalization'});
writetable(Tpct, 'bursting_percent_change.csv');

fprintf('CSV written: bursting_normalized_summary.csv  (%d rows)\n', height(Tnorm));
fprintf('CSV written: bursting_percent_change.csv      (%d rows)\n', height(Tpct));

if export_percell_values && ~isempty(valCols)
    maxLen = max(cellfun(@numel, valCols));
    Mvals  = nan(maxLen, numel(valCols));
    for k = 1:numel(valCols)
        Mvals(1:numel(valCols{k}), k) = valCols{k};
    end
    Tvals = array2table(Mvals, 'VariableNames', ...
                        matlab.lang.makeValidName(valKeys, 'Prefix','v_'));
    writetable(Tvals, 'bursting_normalized_percell_values.csv');
    fprintf('CSV written: bursting_normalized_percell_values.csv  (%d x %d)\n', ...
            maxLen, numel(valCols));
end

%% ------- 9. GRID CSVs: perturbation x parameter ----------------------
% Row order = FIGURE column order, with the extra factors appended last.
% Columns   = Burst duration | Burst frequency | Burst Amplitude
%             (ON duration, 1/OFF duration, amplitude — all normalised to
%              each factor's OWN control, whose mean is 1 by construction).
gridRows = factorsAll;
gridCols = {'Burst duration','Burst frequency','Burst Amplitude'};

gMean_cell = nan(nA,nP); gSEM_cell = nan(nA,nP); gSD_cell = nan(nA,nP);
gP_cell    = nan(nA,nP);
gMean_pool = nan(nA,nP); gSEM_pool = nan(nA,nP); gMed_pool = nan(nA,nP);
gP_pool    = nan(nA,nP);
gN_cellC   = zeros(nA,nP); gN_cellD  = zeros(nA,nP);
gN_burstC  = zeros(nA,nP); gN_burstD = zeros(nA,nP);

for i = 1:nA
    for r = 1:nP
        % ---- per-cell averages (t-test file) ----
        a = clean(avg_C{r,i});   b = clean(avg_D{r,i});
        gN_cellC(i,r) = numel(a);  gN_cellD(i,r) = numel(b);

        [~, bn, ~, ok] = normPair(a, b, invertNorm(r));
        if ok
            gMean_cell(i,r) = mean(bn);
            gSD_cell(i,r)   = std(bn);
            gSEM_cell(i,r)  = std(bn)/sqrt(numel(bn));
        end
        if numel(a) > 1 && numel(b) > 1
            [~, gP_cell(i,r)] = ttest2(a, b, 'Vartype','unequal');
        end

        % ---- pooled bursts (KS file) ----
        x1 = clean(pool_C{r,i});  x2 = clean(pool_D{r,i});
        gN_burstC(i,r) = numel(x1);  gN_burstD(i,r) = numel(x2);

        [~, y2, ~, ok2] = normPair(x1, x2, invertNorm(r));
        if ok2
            gMean_pool(i,r) = mean(y2);
            gSEM_pool(i,r)  = std(y2)/sqrt(numel(y2));
            gMed_pool(i,r)  = median(y2);
        end
        if numel(x1) > 1 && numel(x2) > 1
            [~, gP_pool(i,r)] = kstest2(x1, x2);   % KS on the raw pooled values
        end
    end
end

% ---- file 1: per-cell normalised means + Welch t-test ----
b1 = gridBlock('Normalised mean - treated (each factor vs its own control mean = 1)', gMean_cell, 'num');
b1(2) = gridBlock('SEM (normalised units)',                       gSEM_cell,  'num');
b1(3) = gridBlock('SD (normalised units)',                        gSD_cell,   'num');
b1(4) = gridBlock('p-value - Welch t-test on per-cell averages',  gP_cell,    'p');
b1(5) = gridBlock('n cells - control',                            gN_cellC,   'int');
b1(6) = gridBlock('n cells - treated',                            gN_cellD,   'int');
writeGrid('bursting_grid_percell_ttest.csv', gridRows, gridCols, b1, ...
    'Per-cell averages. Burst duration = ON; Burst frequency = 1/OFF; values normalised to each factor own control mean (control = 1). Statistic: two-sample Welch t-test across cells.');

% ---- file 2: pooled-burst normalised means + KS test ----
b2 = gridBlock('Normalised mean - treated (pooled bursts vs own control mean = 1)', gMean_pool, 'num');
b2(2) = gridBlock('SEM (normalised units)',                        gSEM_pool,  'num');
b2(3) = gridBlock('Normalised median - treated',                   gMed_pool,  'num');
b2(4) = gridBlock('p-value - two-sample KS test on pooled bursts', gP_pool,    'p');
b2(5) = gridBlock('n bursts - control',                            gN_burstC,  'int');
b2(6) = gridBlock('n bursts - treated',                            gN_burstD,  'int');
writeGrid('bursting_grid_pooled_kstest.csv', gridRows, gridCols, b2, ...
    'Pooled bursts (all cells combined). Burst duration = ON; Burst frequency = 1/OFF; values normalised to each factor own control mean (control = 1). Statistic: two-sample Kolmogorov-Smirnov test on the raw pooled distributions - identical to the p-values printed on the CDF panels.');

fprintf('CSV written: bursting_grid_percell_ttest.csv   (%d factors x %d parameters)\n', nA, nP);
fprintf('CSV written: bursting_grid_pooled_kstest.csv   (%d factors x %d parameters)\n', nA, nP);

%% ---------------------- 10. LOCAL FUNCTIONS --------------------------
function M = readBlock(file_path, sheetName, dataRange, sheets, nCol)
% Read one sheet block once; pad/trim to nCol columns; NaN if sheet missing.
    hit = find(strcmpi(sheets, sheetName), 1);
    if isempty(hit)
        warning('Sheet "%s" not found in the workbook — filling with NaN.', sheetName);
        M = nan(1, nCol);  return
    end
    M = readmatrix(file_path, 'Sheet', char(sheets(hit)), 'Range', dataRange);
    if isempty(M), M = nan(1, nCol); return; end
    if size(M,2) < nCol, M(:, size(M,2)+1:nCol) = NaN; end
    M = M(:, 1:nCol);
end

function [an, bn, ref, ok] = normPair(a, b, invert)
% Normalise control (a) and treated (b) to the control mean.
%   invert == false : an = a/ref,          bn = b/ref
%   invert == true  : an = (ref./a)/s,     bn = (ref./b)/s   with
%                     s = mean(ref./a), so mean(an) == 1 exactly and the
%                     treated values sit on the SAME control scale.
% Zero values are dropped (and ref recomputed) only in the inverse case.
    an = []; bn = []; ref = NaN; ok = false;
    a = a(:); b = b(:);
    if isempty(a) || isempty(b), return; end
    if invert
        if any(a == 0) || any(b == 0)
            a = a(a ~= 0);  b = b(b ~= 0);
            if isempty(a) || isempty(b), return; end
        end
        ref = mean(a);
        if ~isfinite(ref) || ref == 0, return; end
        s  = mean(ref ./ a);
        if ~isfinite(s) || s == 0, return; end
        an = (ref ./ a) / s;
        bn = (ref ./ b) / s;
    else
        ref = mean(a);
        if ~isfinite(ref) || ref == 0, return; end
        an = a / ref;
        bn = b / ref;
    end
    ok = true;
end

function blk = gridBlock(titleStr, data, kind)
% One labelled block of a grid CSV. kind = 'num' | 'p' | 'int'.
    blk = struct('title', titleStr, 'data', data, 'kind', kind);
end

function writeGrid(fname, rowNames, colNames, blk, headerLine)
% Write stacked labelled grids: rows = perturbations, cols = parameters.
    fid = fopen(fname, 'w');
    if fid < 0
        warning('Could not open %s for writing.', fname);  return
    end
    if nargin >= 5 && ~isempty(headerLine)
        fprintf(fid, '"%s"\n\n', headerLine);
    end
    for k = 1:numel(blk)
        fprintf(fid, '"%s"\n', blk(k).title);
        fprintf(fid, '"Perturbation"');
        fprintf(fid, ',"%s"', colNames{:});
        fprintf(fid, '\n');
        D = blk(k).data;
        for i = 1:numel(rowNames)
            fprintf(fid, '"%s"', rowNames{i});
            for c = 1:numel(colNames)
                v = D(i, c);
                switch blk(k).kind
                    case 'p'
                        fprintf(fid, ',%s', formatPplain(v));
                    case 'int'
                        fprintf(fid, ',%d', round(v));
                    otherwise
                        if isfinite(v), fprintf(fid, ',%.4f', v);
                        else,           fprintf(fid, ',NaN');
                        end
                end
            end
            fprintf(fid, '\n');
        end
        if k < numel(blk), fprintf(fid, '\n'); end
    end
    fclose(fid);
end

function styleAxis(ax, S)
    set(ax, 'FontName',S.font, 'FontSize',S.fs, 'LineWidth',S.lwAxis, ...
            'TickDir','out', 'Box','off', 'TickLength',S.tickLen, ...
            'XColor',S.axGray, 'YColor',S.axGray, 'Layer','top');
end

function v = niceCeil(x)
% Smallest "nice" number >= x, for 3-tick y-axes.
    if ~isfinite(x) || x <= 0, v = 1; return; end
    e  = floor(log10(x));
    m  = x/10^e;
    st = [1 1.5 2 2.5 3 4 5 6 8 10];
    v  = st(find(st >= m - 1e-9, 1, 'first')) * 10^e;
end

function s = formatP(p)
% Compact TeX p-value that fits a ~14 mm panel at 7 pt.
    if isnan(p)
        s = 'p=n.d.';
    elseif p >= 0.001
        s = sprintf('p=%.3f', p);
    else
        e = floor(log10(p));
        m = p/10^e;
        if m < 1.05
            s = sprintf('p=10^{%d}', e);
        else
            s = sprintf('p=%.0f\\times10^{%d}', m, e);
        end
    end
end

function s = formatPplain(p)
    if isnan(p)
        s = 'NaN';
    elseif p >= 0.001
        s = sprintf('%.4f', p);
    else
        e = floor(log10(p));
        s = sprintf('%.2fE%d', p/10^e, e);
    end
end