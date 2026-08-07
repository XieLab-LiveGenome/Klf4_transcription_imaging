%% BRD4 hub composite figure  (enhancer -> BRD4 hub distance)
% Stacked ksdensity figure with a 25-75 IQR box + median marker above each
% distribution.
%
% Workbook layout expected:
%   - every sheet : PRIMARY distance in column B
%   - "control"   : SECONDARY in column I, TERTIARY in column J
% Requires Statistics & ML Toolbox (ksdensity, prctile), R2020a+ (exportgraphics).

clear; clc;

% ----------------------------- CONFIG -----------------------------
xlsxPath = '/Users/janaa/Desktop/MS2 transcription/primary_secondary_tertiary BRD4 hub_all.xlsx';   % set excel sheet path here

colPrimary   = 'B';    colSecondary = 'I';    colTertiary  = 'J';

sheetKeys   = {'control','g67','RAD21(-)','g67 + RAD21(-)','MED14(-)'};
panelLabels = {'Control','del OCT4/SOX2','RAD21(-)','RAD21(-) + del OCT4/SOX2','MED14(-)'};

xLabelText      = 'Enhancer-BRD4 hub distance (nm)';
yLabelText      = 'Probability density (nm^{-1})';
xLimits         = [0 3000];

% ---- Y-axis limits (one row per panel, top->bottom, matches sheetKeys) ----
%  Use NaN for autoscale on that side, e.g. [0 NaN] to keep zero baseline but
%  let MATLAB pick the top; use [NaN NaN] for full autoscale.
yLimits = [ 0  NaN;   % Control
            0  NaN;   % del OCT4/SOX2
            0  NaN;   % RAD21(-)
            0  NaN;   % RAD21(-) + del OCT4/SOX2
            0  NaN];  % MED14(-)
yTicksAuto = true;

showMedianValue = true;
medValFontSize  = 13;
curveLW         = 3;
outPng          = 'BRD4_hub_composite.png';

boxFrac = 0.24;
rectY   = 0.12;  rectH = 0.62;
medBar  = [0.05 0.80];
medTxtY = 0.88;

navy    = [ 31  45  61]/255;
boxEdge = [154 164 173]/255;
col.ctrlPrimary   = {[127 127 127]/255, [  0   0   0]/255};
col.ctrlSecondary = {[242 194  48]/255, [232 114  44]/255};
col.ctrlTertiary  = {[ 90 168  50]/255, [ 63 138  32]/255};
col.g67           = {[126  63 158]/255, [126  63 158]/255};
col.rad21         = {[232  56  44]/255, [232  56  44]/255};
col.g67rad21      = {[ 91 184 232]/255, [ 91 184 232]/255};
col.med14         = {[139  69  19]/255, [139  69  19]/255};   % brown (saddle brown)
% ------------------------------------------------------------------

realSheets = string(sheetnames(xlsxPath));
fprintf('Sheets found: %s\n', strjoin(realSheets, ', '));

panels = struct('label',{},'series',{});
for i = 1:numel(sheetKeys)
    sh  = resolveSheet(realSheets, sheetKeys{i});
    if sh == ""
        error('Sheet "%s" not found. Available: %s', ...
               sheetKeys{i}, strjoin(realSheets, ', '));
    end
    lbl = panelLabels{i};
    S   = {};
    if strcmp(lbl,'Control')
        S{end+1} = {readCol(xlsxPath,sh,colPrimary),   col.ctrlPrimary{:}};   
        S{end+1} = {readCol(xlsxPath,sh,colSecondary), col.ctrlSecondary{:}}; 
        S{end+1} = {readCol(xlsxPath,sh,colTertiary),  col.ctrlTertiary{:}};  
    else
        cc = pickColour(lbl, col);
        S{end+1} = {readCol(xlsxPath,sh,colPrimary), cc{:}};                   
    end
    panels(i).label  = lbl;
    panels(i).series = S;
end

xmin = xLimits(1);  xmax = xLimits(2);
xgrid = linspace(xmin, xmax, 512);

nCond = numel(panels);
fig = figure('Color','w','Units','pixels','Position',[100 100 1000 270*nCond]);

figLeft = 0.30; figRight = 0.97; figTop = 0.96; figBot = 0.075;
blockGap = 0.045; innerGap = 0.006;
W = figRight - figLeft;  Htot = figTop - figBot;
blockH = (Htot - (nCond-1)*blockGap)/nCond;
boxH   = blockH*boxFrac;
densH  = blockH - boxH - innerGap;

for i = 1:nCond
    blockTop   = figTop - (i-1)*(blockH + blockGap);
    boxBottom  = blockTop  - boxH;
    densTop    = boxBottom - innerGap;
    densBottom = densTop   - densH;

    boxAx  = axes('Parent',fig,'Position',[figLeft boxBottom  W boxH ]); hold(boxAx ,'on');
    densAx = axes('Parent',fig,'Position',[figLeft densBottom W densH]); hold(densAx,'on');

    for k = 1:numel(panels(i).series)
        vals = panels(i).series{k}{1};
        cc   = panels(i).series{k}{2};
        if numel(vals) < 2, continue; end
        f = ksdensity(vals, xgrid);
        plot(densAx, xgrid, f, 'Color',cc, 'LineWidth',curveLW);
    end

    xlim(densAx,[xmin xmax]);

    yl = ylim(densAx);
    yLo = yLimits(i,1);  if isnan(yLo), yLo = 0;      end
    yHi = yLimits(i,2);  if isnan(yHi), yHi = yl(2);  end
    ylim(densAx,[yLo yHi]);

    if yTicksAuto
        set(densAx,'YTickMode','auto');
    else
        set(densAx,'YTick',[]);
    end
    set(densAx,'XColor',navy,'YColor',navy,'LineWidth',1.4,'Layer','top');
    box(densAx,'on');
    if i == nCond
        xlabel(densAx, xLabelText, 'FontSize',12,'Color','k','Interpreter','none');
        set(densAx,'FontSize',10);
    else
        set(densAx,'XTick',[]);
    end

    if strcmp(panels(i).label,'Control')
        lblColor = [0 0 0];
    else
        lblColor = panels(i).series{1}{2};
    end
    text(densAx, -0.03, 0.5, panels(i).label, 'Units','normalized', ...
        'Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle', ...
        'FontSize',13,'Color',lblColor,'Interpreter','none');
    ylabel(densAx, yLabelText, 'FontSize',11,'Color','k');

    xlim(boxAx,[xmin xmax]);  ylim(boxAx,[0 1]);  axis(boxAx,'off');
    for k = 1:numel(panels(i).series)
        vals = panels(i).series{k}{1};
        boxC = panels(i).series{k}{3};
        if isempty(vals), continue; end
        q = prctile(vals, [25 50 75]);  q1 = q(1); med = q(2); q3 = q(3);
        rectangle(boxAx,'Position',[q1 rectY max(q3-q1,eps) rectH], ...
                  'EdgeColor',boxEdge,'LineWidth',1.2);
        line(boxAx,[med med],medBar,'Color',boxC,'LineWidth',6);
        if showMedianValue
            text(boxAx, med, medTxtY, sprintf('%.0f',med), ...
                'HorizontalAlignment','center','VerticalAlignment','bottom', ...
                'FontSize',medValFontSize,'Color',boxC,'Clipping','off');
        end
    end
end

% exportgraphics(fig, outPng, 'Resolution',200);
% fprintf('Saved %s\n', outPng);

%% ----------------------- local functions -----------------------
function v = readCol(xlsxPath, sheet, letter)
    v = readmatrix(xlsxPath, 'Sheet', sheet, 'Range', [letter ':' letter]);
    v = v(:);  v = v(~isnan(v));
end

function sh = resolveSheet(realSheets, wanted)
    key = lower(regexprep(string(wanted), '\s', ''));
    hit = find(lower(regexprep(realSheets, '\s', '')) == key, 1);
    if isempty(hit), sh = ""; else, sh = realSheets(hit); end
end

function cc = pickColour(lbl, col)
    switch lbl
        case 'del OCT4/SOX2',            cc = col.g67;
        case 'RAD21(-)',                 cc = col.rad21;
        case 'MED14(-)',                 cc = col.med14;
        case 'RAD21(-) + del OCT4/SOX2', cc = col.g67rad21;
        otherwise, error('No colour defined for "%s"', lbl);
    end
end