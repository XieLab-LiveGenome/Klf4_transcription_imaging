% % % Boxplot code for visualing SNR data of Klf4 enhancer and MS2 transcription site (TS)

clear
clc

% Enhancer SNR data
enh = [28.47473177
24.14148474
19.92515919
58.51268032
15.90131044
33.37410196
19.52229413
11.06126342
21.83811311
9.617318817
14.55019054
9.304649528
13.39731893
14.55247875
16.48833109
35.42811351
34.64165938
5.766403355
19.16682669
15.66296122
15.77314475
22.34660336
32.3066978
55.80742869
42.76368349
21.1655097
16.48470011
10.92591019
34.4048816
64.86932999
30.76358517
44.81610924
20.16073221
17.43842868
10.64025669
31.61580188
8.160241237

15.33497432
44.36694114
43.23891948
21.09090697
52.66504805
23.49580526
25.30506052
51.45406321
30.36576468
17.23605889
38.16065617
6.352122156
23.73627156
43.82353973
22.0133622
14.88705571
34.33267463
4.603709493
17.37662995
16.65219098
34.20235844
14.68097553
26.68827716
26.55329422
27.47713191
21.07332431
12.38262264

45.28345837
43.42630714
51.41785753
44.64677538
54.04167716
30.09175482
36.36387122
13.44810186
13.47106463
15.37328861
7.836886839
23.55421984
25.48414318
31.16101616
37.04138742
46.51202008
40.47168593
36.36227661
27.65455452
12.34349915
23.92229097
22.81042353
15.87273382
17.7736921
8.474433895
13.38170818
7.466980773
38.8397312
8.184497587
30.84075784


6.026590403
12.26911382
39.82159701
19.41488489
18.08192892
14.49203678
34.330536
13.34843506
23.96079709
42.87916447
25.19569299
33.41067886
15.93901177
34.38732082
15.2426349
33.72028821
25.56963881
22.53256809
33.19794898
28.31470872
24.63387279
19.27365431
22.10827933
36.91360629
22.51327849
17.03554428
19.88879625
19.30511633
18.70751449
39.0816768
26.73910886
20.65516712
33.76057381
23.18999944
12.22515736];

% MS2 spot SNR data 
MS2 = [23.36626646
28.30176659
39.07269486
43.29840636

19.45175847
55.19706784
39.96384111
33.17473743
54.53742904

25.40209865
24.70301095
62.08010227
57.20885985
58.22055461
16.83335546
21.16302993

39.97159112
21.36665104
22.20945488
80.10407016
48.10469048
32.22882752
31.13019352
8.861879465

62.48209225
27.19961702
12.97990291

47.70725917
24.03845327
20.44175579
16.59448716
65.06681361

38.93680166
24.00602564
20.89662796
13.61058124

36.91462065
43.12512791
33.47762407
15.52494637
50.52027349
60.81473467
31.13748618
22.38035067

7.364232968
35.34937589
13.79659843
16.23360028
22.52070772
6.630296499
23.14490984
10.91246052
9.877821198
16.28176315

56.6690086
36.72842973
27.99814394
40.74363598
13.90028929
21.45388683
20.01416039
8.581136162
7.278938306];

% Create figure
figure('Position', [100, 100, 400, 500]);
hold on;

% Create positions for the two groups
positions = [1, 2];
data = {enh, MS2};

% Custom colors - red for crRNA, green for tracrRNA
boxColors = {[0.9, 0.3, 0.3], [0.3, 0.8, 0.4]};  % Red and green
pointColors = {[1, 0.3, 0.3], [0.3, 1, 0.3]};  % Lighter red and green for points

% Plot each boxplot with overlaid points
for i = 1:2
    % Get current data
    currData = data{i};
    pos = positions(i);
    
    % Create boxplot
    bp = boxplot(currData, 'Positions', pos, 'Width', 0.4, ...
                 'Colors', boxColors{i}, 'Symbol', '');
    
    % Fill the box with color
    h = findobj(gca, 'Tag', 'Box');
    patch(get(h(length(h)-i+1), 'XData'), get(h(length(h)-i+1), 'YData'), ...
          boxColors{i}, 'FaceAlpha', 0.6);
    
    % Make boxplot lines thicker
    set(bp, 'LineWidth', 1.5);
    
    % Make the median line bolder
    medianLines = findobj(gca, 'Tag', 'Median');
    set(medianLines(length(medianLines)-i+1), 'LineWidth', 3, 'Color', boxColors{i}*0.6);
    
    % Add individual points with jitter
    jitter = 0.08 * randn(length(currData), 1);
    scatter(pos + jitter, currData, 30, pointColors{i}, 'filled', ...
            'MarkerFaceAlpha', 0.6);
end

% Make whiskers (vertical lines) continuous, not dashed
whiskers = findobj(gca, 'Tag', 'Whisker');
set(whiskers, 'LineStyle', '-', 'LineWidth', 1.5);

% Add p-value annotation
% pval_text = ['{\itp} < 1.0 × 10^{-15}'];
y_max = max([enh; MS2]) + 5;
line([1, 2], [y_max, y_max], 'Color', 'k', 'LineWidth', 1.5);
line([1, 1], [y_max-2, y_max], 'Color', 'k', 'LineWidth', 1.5);
line([2, 2], [y_max-2, y_max], 'Color', 'k', 'LineWidth', 1.5);
% text(1.5, y_max + 15, pval_text, 'HorizontalAlignment', 'center', ...
%      'FontSize', 11);

% Formatting
set(gca, 'XTick', positions, 'XTickLabel', {'Klf4 enhancer', 'Klf4 MS2'}, ...
         'FontSize', 14, 'LineWidth', 1.2);
ylabel('SNR', 'FontSize', 20);
ylim([0, 90]);
xlim([0.5, 2.5]);

% Rotate x-axis labels
xtickangle(45);

% Remove top and right spines
box off;
set(gca, 'TickDir', 'out');

hold off;