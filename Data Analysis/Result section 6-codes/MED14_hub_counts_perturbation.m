%% BRD4 Condensate Size — xxBox Plot + t-tests vs Ctrl

% each treatment/depletion against Ctrl.
%
% Requires: Statistics and Machine Learning Toolbox (ttest2, ranksum, adtest)
% ============================================================
% STEP 1: FILE & RANGE SETUP  (edit these)
% ============================================================
filename  = '/Users/janaa/Desktop/MS2 transcription/Med14_condensate counts.xlsx';  
sheetName = 'Sheet1';                             

% Column range for each condition (rows 5-500)
ranges = {'B2:B500','r2:r500','m2:m500', 'g2:g500'};
labels = {'Ctrl', 'iBRG1', 'iP300', 'JQ1' };

% ---  groups to include ( Ctrl is index 1, first) ---
pick   = [1 2 3 4];          
ranges = ranges(pick);
labels = labels(pick);

nGroups = numel(labels);   


% Which p-value drives the significance stars: p_ttest (raw) or p_ttest_bonf
useBonferroniStars = false;   % set true to annotate with Bonferroni-corrected p

% ============================================================
% STEP 2: READ & CLEAN DATA
% ============================================================
data = cell(1, nGroups);
for i = 1:nGroups
    col = readmatrix(filename, 'Sheet', sheetName, 'Range', ranges{i});
    % --- older MATLAB fallback (pre-R2019a): ---
    % col = xlsread(filename, sheetName, ranges{i});
    col = col(:);                 % force column vector
    col = col(~isnan(col));       % drop blank/non-numeric cells
    data{i} = col;
    fprintf('%-15s n = %d\n', labels{i}, numel(col));
end

% ============================================================
% STEP 3: DESCRIPTIVE STATS + PAIRWISE TESTS vs CTRL
% ============================================================
ctrl    = data{1};
n_ctrl  = numel(ctrl);
mean_ctrl = mean(ctrl);
sd_ctrl   = std(ctrl);
nComp   = nGroups - 1;            % 6 comparisons against Ctrl

Condition  = labels(2:end)';
N_group    = zeros(nComp,1);
Mean_group = zeros(nComp,1);
SD_group   = zeros(nComp,1);
p_normal   = nan(nComp,1);
p_ttest    = nan(nComp,1);
p_ranksum  = nan(nComp,1);
cohens_d   = nan(nComp,1);

for k = 1:nComp
    g  = data{k+1};
    ng = numel(g);
    N_group(k)    = ng;
    if ng < 2,  continue;  end
    Mean_group(k) = mean(g);
    SD_group(k)   = std(g);

    if ng >= 4
        [~, p_normal(k)] = adtest(g);          % Anderson-Darling normality
    end

    % Two-sample t-test (classic, equal-variance — like your original).
    % For unequal variances use: ttest2(ctrl, g, 'Vartype','unequal')  (Welch)
    [~, p_ttest(k)] = ttest2(ctrl, g);

    % Non-parametric alternative (robust to non-normal distributions)
    p_ranksum(k) = ranksum(ctrl, g);

    % Cohen's d (pooled SD)
    sp = sqrt(((n_ctrl-1)*sd_ctrl^2 + (ng-1)*SD_group(k)^2) / (n_ctrl+ng-2));
    cohens_d(k) = (mean_ctrl - Mean_group(k)) / sp;
end

% Bonferroni correction across the 6 comparisons
p_ttest_bonf   = min(p_ttest   * nComp, 1);
p_ranksum_bonf = min(p_ranksum * nComp, 1);

% ============================================================
% STEP 4: STYLED BOX PLOT WITH SCATTER OVERLAY
% ============================================================
% Per-group colors (Ctrl = gray)
colors = [0.50 0.50 0.50;    % Ctrl
          0.90 0.40 0.40;    % JQ1
          0.95 0.61 0.33;    % iP300
          0.85 0.72 0.25;    % iBRG1
          0.45 0.70 0.45;    % MED14 depl.
          0.40 0.55 0.80;    % RAD21 depl.
          0.60 0.45 0.75];   % del OCT4/SOX2

% Stack data for boxplot
all_data = []; grp = [];
for i = 1:nGroups
    all_data = [all_data; data{i}];
    grp      = [grp; i*ones(numel(data{i}),1)];
end

figure('Position', [100 100 1150 650], 'Color', 'w');
h = boxplot(all_data, grp, 'Labels', labels, 'Widths',0.6, 'Symbol', '');
hold on;

set(h, 'LineWidth', 1.8);                                  % whiskers/caps
set(findobj(gca,'Tag','Median'), 'Color','k', 'LineWidth', 5.0);  % median

% Color each box outline by group (sort by x-position = version-robust)
boxes = findobj(gca, 'Tag', 'Box');
xpos  = arrayfun(@(b) mean(get(b,'XData')), boxes);
[~, ord] = sort(xpos);
boxes = boxes(ord);
for i = 1:nGroups
    set(boxes(i), 'Color', colors(i,:)*0.7, 'LineWidth', 2.5);
end

% Hollow jittered scatter
jitter = 0.25;
for i = 1:nGroups
    x = i + jitter*(rand(numel(data{i}),1) - 0.5);
    scatter(x, data{i}, 60, 'o', ...
        'MarkerEdgeColor', colors(i,:), ...
        'MarkerFaceColor', 'none', ...
        'MarkerEdgeAlpha', 0.95, ...
        'LineWidth', 1.2);
end

% Axis formatting
ylabel('# of BRD4 hubs/nucleus', 'FontSize', 22, 'FontWeight', 'bold');
ymax = max(all_data);
% ylim([0, ymax*1.18]);
ylim([0 120]);
set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'LineWidth', 1.5, ...
         'XTickLabelRotation', 30);
box on;

% Significance stars (each treatment vs Ctrl)
annotP = p_ttest;
if useBonferroniStars, annotP = p_ttest_bonf; end
yl = ylim;
for k = 1:nComp
    pj = annotP(k);
    if     pj < 0.001, s = '***';
    elseif pj < 0.01,  s = '**';
    elseif pj < 0.05,  s = '*';
    else,              s = 'ns';
    end
    text(k+1, yl(2)*0.95, s, 'HorizontalAlignment','center', ...
        'FontSize', 17, 'FontWeight','bold');
end
hold off;

% Optional: high-res export for figures
% exportgraphics(gcf, 'BRD4_condensate_boxplot.png', 'Resolution', 600);

% ============================================================
% STEP 5: PRINT RESULTS
% ============================================================
fprintf('\n=== Ctrl: Mean = %.2f, SD = %.2f, n = %d ===\n', mean_ctrl, sd_ctrl, n_ctrl);

Results = table(Condition, N_group, Mean_group, SD_group, ...
                p_normal, p_ttest, p_ttest_bonf, p_ranksum, p_ranksum_bonf, cohens_d, ...
    'VariableNames', {'Condition','N','Mean','SD','p_normal', ...
                      'p_ttest','p_ttest_Bonf','p_ranksum','p_ranksum_Bonf','Cohens_d'});
disp(Results)