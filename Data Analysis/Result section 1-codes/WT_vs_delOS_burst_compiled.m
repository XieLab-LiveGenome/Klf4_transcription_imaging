%% ================= Bursting CDF panel figure del O/S =================
%%%cumulative-frequency plots (ON duration, OFF duration, burst amplitude)

clear
clc

%% ---------------- User settings ----------------
file_path = '/Volumes/xiel2lab/Xie_lab_Manuscripts/2026_MS2_Paper/Image analysis data/Results section 1/Excel data/ctrl vs g67 bursting % 7-14-2025.xlsx';

sheet1 = 'WT burst parameters pooled';   label1 = 'WT';             color1 = [0 0 0];
sheet2 = 'del OS burst parameters pooled';   label2 = 'del OCT4/SOX2';  color2 = [1 0 0];

lw     = 5;                % CDF curve width
fs     = 20;               % font size
axlw   = 2;                % axis box width
pcol   = [0.55 0.55 0.55]; % colour of the p-value text


legPos  = [0.32 0.68];     % top line of the two-line legend, panel 1 only
pPos    = [0.28 0.16];     % p-value text, every panel

%% ---------------- Load data ----------------
readcol = @(sheet,col) rmmissing(readmatrix(file_path, 'Sheet', sheet, ...
                                 'Range', [col '2:' col '200']));

A.ON  = readcol(sheet1,'D');   B.ON  = readcol(sheet2,'D');
A.OFF = readcol(sheet1,'E');   B.OFF = readcol(sheet2,'E');
A.AMP = readcol(sheet1,'H');   B.AMP = readcol(sheet2,'H');
A.SZ  = readcol(sheet1,'I');   B.SZ  = readcol(sheet2,'I');

%% ---------------- Panel definitions ----------------
% To add burst size back, just append a 4th entry (log-scaled example at the
% bottom of this block) — everything downstream is generic.
k = 0;
k = k+1; P(k) = mkpanel(A.ON,  B.ON,  'ON duration (min)',  [0 180], 0:60:180, 0.4, false);
k = k+1; P(k) = mkpanel(A.OFF, B.OFF, 'OFF duration (min)', [0 240], 0:60:240, 0.4, false);
k = k+1; P(k) = mkpanel(A.AMP, B.AMP, 'Burst amplitude',    [0 300], 0:100:300, 0.4, false);
% k = k+1; P(k) = mkpanel(A.SZ, B.SZ, 'Burst size', [20 20000], [10 100 1000 10000], 0.1, true);

nP = numel(P);

%% ---------------- Figure ----------------
figure('Color','w');
set(gcf,'Units','inches','Position',[1 1 3.6*nP 3.9]);

Stats_all = nan(1,nP);

for i = 1:nP
    x1 = P(i).d1;  x2 = P(i).d2;

    % --- KS test ---
    [~, pval] = kstest2(x1, x2, 'Alpha', 0.05);
    Stats_all(i) = pval;

    % --- smooth CDFs, extended across the full x range ---
    [f1, xi1] = ksdensity(x1, 'function','cdf', 'Bandwidth', P(i).bw);
    [f2, xi2] = ksdensity(x2, 'function','cdf', 'Bandwidth', P(i).bw);

    xq  = linspace(P(i).xlim(1), P(i).xlim(2), 2000);
    y1  = min(max(interp1(xi1, f1, xq, 'linear','extrap'), 0), 1);
    y2  = min(max(interp1(xi2, f2, xq, 'linear','extrap'), 0), 1);

    % --- plot ---
    ax = subplot(1, nP, i);
    plot(xq, y1, 'LineWidth', lw, 'Color', color1); hold on
    plot(xq, y2, 'LineWidth', lw, 'Color', color2);

    xlabel(P(i).xlabel, 'FontSize', fs, 'FontWeight','bold');
    if i == 1
        ylabel('Cumulative frequency', 'FontSize', fs, 'FontWeight','bold');
    else
        set(ax,'YTickLabel',[]);          % keep ticks, drop repeated numbers
    end

    xlim(P(i).xlim);  ylim([0 1]);
    xticks(P(i).xticks); yticks(0:0.2:1);
    if P(i).logx, set(ax,'XScale','log'); end

    set(ax, 'FontSize', fs, 'FontWeight','bold', 'LineWidth', axlw, ...
            'Box','on', 'TickDir','in', 'Layer','top');

    % --- condition labels (panel 1 only) ---
    if i == 1
        text(legPos(1), legPos(2),      label1, 'Units','normalized', ...
            'Color', color1, 'FontSize', fs, 'FontWeight','bold');
        text(legPos(1), legPos(2)-0.11, label2, 'Units','normalized', ...
            'Color', color2, 'FontSize', fs, 'FontWeight','bold');
    end

    % --- p-value ---
    text(pPos(1), pPos(2), pstring(pval), 'Units','normalized', ...
        'Color', pcol, 'FontSize', fs, 'FontWeight','bold', ...
        'Interpreter','tex');
end

% Tighten the row a little so panels sit close, as in the reference figure
set(findall(gcf,'Type','axes'), 'Units','normalized');
tightenRow(gcf, nP);

%% ---------------- Save (uncomment) ----------------
% exportgraphics(gcf, 'bursting_CDF_panel.pdf', 'ContentType','vector');
% exportgraphics(gcf, 'bursting_CDF_panel.tif', 'Resolution', 600);


%% ================= local functions =================
function s = mkpanel(d1, d2, lab, xl, xt, bw, logx)
    s.d1 = d1;  s.d2 = d2;  s.xlabel = lab;
    s.xlim = xl;  s.xticks = xt;  s.bw = bw;  s.logx = logx;
end

function str = pstring(p)
% Format a p-value as p=1.52*10^-10 with a TeX superscript.
    if p == 0
        str = 'p<10^{-300}';
        return
    end
    e = floor(log10(p));
    m = p / 10^e;
    if m >= 9.995      % rounding guard, e.g. 9.999e-5 -> 1.00e-4
        m = m/10;  e = e+1;
    end
    if e == 0
        str = sprintf('p=%.2f', m);
    else
        str = sprintf('p=%.2f*10^{%d}', m, e);
    end
end

function tightenRow(fh, nP)
% Evenly space the panels with a small gap, leaving room for labels.
    ax = flipud(findall(fh,'Type','axes'));
    left = 0.085; right = 0.015; gap = 0.035; bottom = 0.20; height = 0.72;
    w = (1 - left - right - gap*(nP-1)) / nP;
    for i = 1:numel(ax)
        set(ax(i), 'Position', [left + (i-1)*(w+gap), bottom, w, height]);
    end
end