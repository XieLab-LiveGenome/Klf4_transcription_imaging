clear
clc

file_path = '/Users/janaa/Desktop/MS2 transcription/4 color compiled/4 COLOR BRD4 MED14.xlsx';

MS2 = readmatrix(file_path,'Sheet','MS2 vs coloc', 'Range','b4:b252');
Prim = readmatrix(file_path,'Sheet','MS2 vs coloc', 'Range','d4:d252');
Sec = readmatrix(file_path,'Sheet','MS2 vs coloc', 'Range','e4:e252');
Tert = readmatrix(file_path,'Sheet','MS2 vs coloc', 'Range','f4:f252');

MS2 = MS2(~isnan(MS2));
Prim = Prim(~isnan(Prim));   % Manders coefficient for BRD4/MED14 colocalization of primary condensate 
Sec = Sec(~isnan(Sec));      % Manders coefficient for BRD4/MED14 colocalization of secondary condensate 
Tert = Tert(~isnan(Tert));   % Manders coefficient for BRD4/MED14 colocalization of tertiary condensate 

thr = 4000;  % MS2 intensity threhold for ON/OFF

Prim_ON = Prim(MS2>=thr);
Prim_OFF = Prim(MS2<thr);

Sec_ON = Sec(MS2>=thr);
Sec_OFF = Sec(MS2<thr);

Tert_ON = Tert(MS2>=thr);
Tert_OFF = Tert(MS2<thr);


%% ---- Plot: 1x3 Manders BRD4/MED14, ON vs OFF ----
figure('Position',[100 100 1200 400],'Color','w');
titles = {'Primary','Secondary','Tertiary'};
data_ON  = {Prim_ON,  Sec_ON,  Tert_ON};
data_OFF = {Prim_OFF, Sec_OFF, Tert_OFF};

% Define colors
cOFF = [0.6 0.6 0.6];       % Grey for OFF
cON  = [0.2 0.8 0.2];      % Original Red for ON scatter 
cON_med = [0.2 0.6 0.2];    % Green for ON median line
cOFF_med = [0.3 0.3 0.3];    % Green for ON median line


% Pre-allocate array to store Wilcoxon rank-sum p-values
p_values_wilcoxon = zeros(1, 3);

for k = 1:3
    ax = subplot(1,3,k); hold on;
    dOFF = data_OFF{k};
    dON  = data_ON{k};
    
    % open circle scatter
    jit = 0.15;
    scatter(ones(size(dOFF)) + jit*(rand(size(dOFF))-0.5), dOFF, 44, cOFF, 'LineWidth',1.2);
    scatter(2*ones(size(dON)) + jit*(rand(size(dON))-0.5), dON, 44, cON, 'LineWidth',1.2);
    
    % box overlay
    bp = boxplot([dOFF; dON], [ones(size(dOFF)); 2*ones(size(dON))], ...
        'Positions',[1 2],'Widths',0.35,'Symbol','','Colors','k');
    set(bp,{'linew'},{1.2});
    
% modify median lines (Row 6 of the boxplot handle matrix contains medians)
    set(bp(6,1), 'Color', cOFF_med, 'LineWidth', 5); % OFF median to grey and thicker
    set(bp(6,2), 'Color', cON_med, 'LineWidth', 5);  % ON median to green and thicker
    
    % stats: Wilcoxon rank-sum test
    p = ranksum(dOFF, dON);
    
    % Store the p-value
    p_values_wilcoxon(k) = p;
    
    ymax = max([dOFF; dON]);
    yt = ymax * 1.08;
    plot([1 2],[yt yt],'k-','LineWidth',1);
    
    if p < 0.001
        pstr = '***';
    elseif p < 0.01
        pstr = '**';
    elseif p < 0.05
        pstr = '*';
    else
        pstr = 'n.s.';
    end
    
    text(1.5, yt*1.03, pstr, 'HorizontalAlignment','center','FontSize',13,'FontWeight','bold');
    
    % n in x-tick labels
    lblOFF = sprintf('OFF');
    lblON  = sprintf('ON');
    set(ax,'XTick',[1 2],'XTickLabel',{lblOFF, lblON},...
        'FontSize',14,'FontWeight','bold','TickDir','out','Box','off','LineWidth',1.2);
    xlim([0.4 2.6]);
    ylim([0 yt*1.12]);
    ylabel('Manders coeff. (BRD4/MED14)','FontSize',16,'FontWeight','bold');
    title(titles{k},'FontSize',17,'FontWeight','bold');
end

% Optional: Display the stored p-values in the command window
disp('Wilcoxon rank-sum p-values (Primary, Secondary, Tertiary):');
disp(p_values_wilcoxon);

% sgtitle('BRD4/MED14 Manders: Gene ON vs OFF','FontSize',18,'FontWeight','bold');