clear; clc;

%% ---- Point to PARENT folder containing the 4 category subfolders -------
parentFolder = '/Users/janaa/Desktop/MS2 transcription/Compiled data/Condensate BRD4-MED14 dual';
categories   = {'Ctrl','iBRG1','iP300','JQ1'};   % subfolder names / bar order
nCat         = numel(categories);

% ---- Per-cell table schema: {source field in results, output column} ----
pc_map = {
    'BRD4_count_per_nuc'        'BRD4_count'
    'MED14_count_per_nuc'       'MED14_count'
    'BRD4_coloc_frac_per_nuc'   'BRD4_coloc_frac'
    'BRD4_isol_frac_per_nuc'    'BRD4_isol_frac'
    'MED14_coloc_frac_per_nuc'  'MED14_coloc_frac'
    'MED14_isol_frac_per_nuc'   'MED14_isol_frac'
    'M2_coloc_BRD4_per_nuc'     'M2_coloc_BRD4'
    'M2_isol_BRD4_per_nuc'      'M2_isol_BRD4'
    'M1_coloc_MED14_per_nuc'    'M1_coloc_MED14'
    'M1_isol_MED14_per_nuc'     'M1_isol_MED14'
};

%% ---- Loop categories: extract per-cell table, compute mean +/- SEM -----
nCells      = zeros(nCat,1);
med_coloc_m = nan(nCat,1);  med_coloc_s = nan(nCat,1);
med_isol_m  = nan(nCat,1);  med_isol_s  = nan(nCat,1);
brd_coloc_m = nan(nCat,1);  brd_coloc_s = nan(nCat,1);
brd_isol_m  = nan(nCat,1);  brd_isol_s  = nan(nCat,1);
allTables   = {};

for c = 1:nCat
    folder = fullfile(parentFolder, categories{c});
    if ~isfolder(folder)
        warning('Category folder not found: %s', folder);
        continue;
    end
    T = extractCellFractions(folder, pc_map);
    nCells(c) = height(T);
    fprintf('%-8s : %d cells from %d file(s)\n', categories{c}, height(T), ...
            numel(dir(fullfile(folder,'*.mat'))));

    if height(T) == 0
        warning('No per-cell data in category %s.', categories{c});
        continue;
    end

    T.Category = repmat(string(categories{c}), height(T), 1);
    allTables{end+1} = T; %#ok<SAGROW>

    [med_coloc_m(c), med_coloc_s(c)] = meanSEM(T.MED14_coloc_frac);
    [med_isol_m(c),  med_isol_s(c)]  = meanSEM(T.MED14_isol_frac);
    [brd_coloc_m(c), brd_coloc_s(c)] = meanSEM(T.BRD4_coloc_frac);
    [brd_isol_m(c),  brd_isol_s(c)]  = meanSEM(T.BRD4_isol_frac);
end

%% ---- Save combined per-cell table + category summary ------------------
if ~isempty(allTables)
    cellFractions = vertcat(allTables{:});
    cellFractions.Category = categorical(cellFractions.Category, categories);  % keep order
    cellFractions = movevars(cellFractions, 'Category', 'Before', 'FileName');
    save(fullfile(parentFolder,'per_cell_fractions_all_categories.mat'), 'cellFractions');
    writetable(cellFractions, fullfile(parentFolder,'per_cell_fractions_all_categories.csv'));
end

catSummary = table(categories(:), nCells, ...
    100*med_coloc_m, 100*med_coloc_s, 100*med_isol_m, 100*med_isol_s, ...
    100*brd_coloc_m, 100*brd_coloc_s, 100*brd_isol_m, 100*brd_isol_s, ...
    'VariableNames', {'Category','nCells', ...
    'MED14_coloc_pct','MED14_coloc_SEM','MED14_isol_pct','MED14_isol_SEM', ...
    'BRD4_coloc_pct','BRD4_coloc_SEM','BRD4_isol_pct','BRD4_isol_SEM'});
writetable(catSummary, fullfile(parentFolder,'category_coloc_isol_summary.csv'));

%% ---- Console summary ---------------------------------------------------
fprintf('\n==========================================================================\n');
fprintf('  category   n     MED14 coloc%%      MED14 isol%%       BRD4 coloc%%       BRD4 isol%%\n');
fprintf('--------------------------------------------------------------------------\n');
for c = 1:nCat
    fprintf('  %-8s %4d   %5.1f +/- %-4.1f   %5.1f +/- %-4.1f   %5.1f +/- %-4.1f   %5.1f +/- %-4.1f\n', ...
        categories{c}, nCells(c), ...
        100*med_coloc_m(c), 100*med_coloc_s(c), 100*med_isol_m(c), 100*med_isol_s(c), ...
        100*brd_coloc_m(c), 100*brd_coloc_s(c), 100*brd_isol_m(c), 100*brd_isol_s(c));
end
fprintf('==========================================================================\n');

%% ========================================================================
%  STATISTICS: colocalized fraction vs Control (two-sample t-test)
%  ========================================================================
if exist('cellFractions','var') && ~isempty(cellFractions)
    ctrlName = 'Ctrl';
    testCats = categories(~strcmp(categories, ctrlName));   % non-control groups

    statVars = {'MED14_coloc_frac','BRD4_coloc_frac'};
    statLbl  = {'MED14 coloc','BRD4 coloc'};

    ctrlMask = cellFractions.Category == ctrlName;
    statRows = {};

    fprintf('\n==========================================================================\n');
    fprintf('  Welch two-sample t-test vs %s  (per-cell colocalized fraction)\n', ctrlName);
    fprintf('--------------------------------------------------------------------------\n');
    fprintf('  variable     group     Ctrl%%    grp%%    dPct     p-value    sig\n');
    fprintf('--------------------------------------------------------------------------\n');

    for v = 1:numel(statVars)
        xc = cellFractions.(statVars{v})(ctrlMask);  xc = xc(~isnan(xc));
        for c = 1:numel(testCats)
            m  = cellFractions.Category == testCats{c};
            xt = cellFractions.(statVars{v})(m);  xt = xt(~isnan(xt));
            if isempty(xc) || isempty(xt)
                p = NaN;  tstat = NaN;
            else
                [~, p, ~, st] = ttest2(xc, xt, 'Vartype','unequal');  % Welch's t
                tstat = st.tstat;
            end
            dPct = 100*(mean(xt) - mean(xc));
            fprintf('  %-11s  %-8s  %5.1f   %5.1f   %+6.1f   %9.2g   %s\n', ...
                statLbl{v}, testCats{c}, 100*mean(xc), 100*mean(xt), dPct, p, sigStars(p));
            statRows(end+1,:) = {statLbl{v}, testCats{c}, numel(xc), numel(xt), ...
                                 100*mean(xc), 100*mean(xt), dPct, tstat, p, sigStars(p)}; %#ok<SAGROW>
        end
    end
    fprintf('==========================================================================\n');

    statsTable = cell2table(statRows, 'VariableNames', ...
        {'Variable','Group','n_ctrl','n_grp','Ctrl_pct','Grp_pct', ...
         'Delta_pct','tstat','p_value','sig'});
    writetable(statsTable, fullfile(parentFolder,'coloc_ttest_vs_control.csv'));
    fprintf('Saved t-test results to coloc_ttest_vs_control.csv\n');
end

%% ========================================================================
%  STACKED BAR PLOT  (row 1 = MED14, row 2 = BRD4)
%  ========================================================================
x = 1:nCat;
xlabels = cell(nCat,1);
for c = 1:nCat, xlabels{c} = sprintf('%s\\newline(n=%d)', categories{c}, nCells(c)); end

% colors: colocalized = saturated (bottom), isolated = light tint (top)
med_coloc_col = [0.20 0.60 0.20];   med_isol_col = [0.75 0.90 0.75];   % greens (MED14)
brd_coloc_col = [0.70 0.15 0.70];   brd_isol_col = [0.90 0.75 0.90];   % magentas (BRD4)

figure('Position',[400 80 640 740]);
tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% ---- MED14 (row 1) ----
nexttile;
Ymed = 100*[med_coloc_m, med_isol_m];                 % nCat x 2  (coloc | isol)
b = bar(x, Ymed, 'stacked', 'BarWidth',0.6, 'EdgeColor',[0.2 0.2 0.2], 'LineWidth',0.5);
b(1).FaceColor = med_coloc_col;   b(2).FaceColor = med_isol_col;
hold on;
% SEM whisker at the coloc/isol boundary (= mean colocalized %)
errorbar(x, 100*med_coloc_m, 100*med_coloc_s, 'k', 'LineStyle','none', ...
         'LineWidth',1.2, 'CapSize',10);
ylim([0 100]); ylabel('MED14 (% of condensates)');
set(gca,'XTick',x,'XTickLabel',xlabels);
title('MED14: colocalized vs isolated');
legend({'Colocalized','Isolated'}, 'Location','eastoutside', 'Box','off');

% ---- BRD4 (row 2) ----
nexttile;
Ybrd = 100*[brd_coloc_m, brd_isol_m];
b2 = bar(x, Ybrd, 'stacked', 'BarWidth',0.6, 'EdgeColor',[0.2 0.2 0.2], 'LineWidth',0.5);
b2(1).FaceColor = brd_coloc_col;  b2(2).FaceColor = brd_isol_col;
hold on;
errorbar(x, 100*brd_coloc_m, 100*brd_coloc_s, 'k', 'LineStyle','none', ...
         'LineWidth',1.2, 'CapSize',10);
ylim([0 100]); ylabel('BRD4 (% of condensates)');
set(gca,'XTick',x,'XTickLabel',xlabels);
title('BRD4: colocalized vs isolated');
legend({'Colocalized','Isolated'}, 'Location','eastoutside', 'Box','off');

title(tl, sprintf('Colocalized vs isolated (%d categories, mean \\pm SEM)', nCat), ...
      'FontWeight','bold');
set(findall(gcf,'-property','FontSize'),'FontSize',13,'FontName','Arial');

% ---- Save figure --------------------------------------------------------
savefig(gcf, fullfile(parentFolder,'coloc_isol_stacked_bars.fig'));
exportgraphics(gcf, fullfile(parentFolder,'coloc_isol_stacked_bars.png'), 'Resolution',300);
fprintf('\nSaved figure + summary tables to:\n  %s\n', parentFolder);


%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%  ========================================================================
function T = extractCellFractions(folder, pc_map)
% Build a per-cell table (one row per nucleus) for every .mat file in a
% folder, using the robust guarded extraction. Returns an empty table if
% nothing usable is found.
    matFiles = dir(fullfile(folder, '*.mat'));
    perCell  = {};
    for f = 1:numel(matFiles)
        r = loadResultsStruct(fullfile(folder, matFiles(f).name));
        if isempty(r), continue; end
        n = inferNumNuclei(r, pc_map);
        if n == 0, continue; end
        Tf = table();
        Tf.FileName = repmat(string(matFiles(f).name), n, 1);
        Tf.CellID   = (1:n)';
        for k = 1:size(pc_map,1)
            src = pc_map{k,1};  col = pc_map{k,2};
            v = getcol(r, src);
            if numel(v) == n
                Tf.(col) = v;
            else
                Tf.(col) = nan(n,1);
            end
        end
        perCell{end+1} = Tf; %#ok<AGROW>
    end
    if isempty(perCell), T = table(); else, T = vertcat(perCell{:}); end
end

function [m, s] = meanSEM(x)
% Mean and standard error of the mean, ignoring NaNs.
    x = x(~isnan(x));
    if isempty(x), m = NaN; s = NaN; return; end
    m = mean(x);
    s = std(x) / sqrt(numel(x));
end

function r = loadResultsStruct(matPath)
% Load the analysis struct from a .mat file. Prefers a variable named
% 'results'; otherwise returns the first struct-valued variable found.
    S = load(matPath);
    if isfield(S, 'results') && isstruct(S.results)
        r = S.results;
        return;
    end
    r = [];
    fn = fieldnames(S);
    for i = 1:numel(fn)
        if isstruct(S.(fn{i}))
            r = S.(fn{i});
            return;
        end
    end
end

function v = getcol(r, name)
% Return results.(name) as a column vector, or [] if the field is
% absent/empty. Guards every field access so version drift never crashes.
    if isstruct(r) && isfield(r, name) && ~isempty(r.(name))
        v = double(r.(name));
        v = v(:);
    else
        v = [];
    end
end

function n = inferNumNuclei(r, pc_map)
% Determine nuclei count robustly: trust results.numNuclei when present and
% consistent with the per-nucleus vector lengths; otherwise take the most
% common per-nucleus vector length.
    lens = [];
    for k = 1:size(pc_map,1)
        v = getcol(r, pc_map{k,1});
        if ~isempty(v), lens(end+1) = numel(v); end %#ok<AGROW>
    end
    nDeclared = [];
    if isfield(r, 'numNuclei') && ~isempty(r.numNuclei)
        nDeclared = double(r.numNuclei);
    end
    if isempty(lens)
        if ~isempty(nDeclared), n = nDeclared; else, n = 0; end
        return;
    end
    nData = mode(lens);
    if ~isempty(nDeclared) && nDeclared == nData
        n = nDeclared;
    else
        n = nData;
    end
end

function s = sigStars(p)
% Significance stars for a p-value.
    if     isnan(p),    s = 'n/a';
    elseif p < 0.001,   s = '***';
    elseif p < 0.01,    s = '**';
    elseif p < 0.05,    s = '*';
    else,               s = 'ns';
    end
end