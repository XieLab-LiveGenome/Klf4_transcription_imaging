%% ========================================================================
%  Spatial-scale panel: E-P separation / BRD4 hub size
%  LINEAR distance axis (0-1000 nm), one graded bar per species,
%  deepest colour at the mean.
%
%  Fill in MU and SD below (in NANOMETRES). Everything else is derived.
%  ========================================================================

clear; 
clc

endash = char(8211);          % en dash for 'E-P'
times  = char(215);           % multiplication sign
Delta  = char(916);           % capital delta, for the difference label

%% ---------------------- USER INPUT --------------------------------------
% mu and sd in NANOMETRES. lab is the text printed to the right of each bar.

S(1).name = ['E' endash 'P separation (control)'];
S(1).mu   = 283.3;                     % <-- mean distance (nm)
S(1).sd   = 142.3;                     % <-- SD (nm)
S(1).col  = [0.70 0.7 0.70];
S(1).lab  = '';                        % '' = auto-format from mu

S(2).name = ['E' endash 'P separation (cohesin(-))'];
S(2).mu   = 414.7;                     % <-- mean distance (nm)
S(2).sd   = 230.7;                     % <-- SD (nm)
S(2).col  = [0.80 0.208 0.80];
S(2).lab  = '';

S(3).name = 'BRD hub size';
S(3).mu   = 467.8;                     % <-- mean size (nm)
S(3).sd   = 173.2;                     % <-- SD (nm)
S(3).col  = [0.180 0.525 0.757];
S(3).lab  = '';

%% ---------------------- OPTIONS -----------------------------------------
XLIM      = [0 800];     % axis range (nm), linear
XTICKS    = 0:200:800;   % major ticks (nm)
NSD       = 1.5;            % bar extends mu +/- NSD*sd
NSEG      = 300;          % gradient segments per bar (higher = smoother)
EDGEW     = 0.05;         % colour weight at the very ends of the bar
                          %   0 = pure white tips, 1 = no fade at all
GAMMA     = 1.25;         % >1 pushes more of the bar toward the light end
                          %   (1 = plain linear ramp from mean to tip)
SHOWSEP   = true;         % draw the separation arrow between rows 1-2
SEPROWS   = [1 2];        % which two rows the arrow compares
SEPMODE   = 'fold';       % 'fold'  -> "~1.5x separation"
                          % 'delta' -> "D ~131 nm"
GUIDES    = [250 500];    % dashed vertical guides (nm), labelled automatically
FS        = 10;           % base font size
FIGSIZE   = [18 7.2];     % cm, width x height
OUTFILE   = 'distance_panel';    % written as .pdf and .png

n  = numel(S);
xr = diff(XLIM);          % axis span, used for all label offsets

%% ---------------------- FIGURE ------------------------------------------
fig = figure('Color','w','Units','centimeters', ...
             'Position',[2 2 FIGSIZE(1) FIGSIZE(2)]);
ax = axes('Parent',fig,'Units','normalized','Position',[0.20 0.13 0.63 0.62]);
hold(ax,'on'); box(ax,'off');

set(ax,'XScale','linear','XLim',XLIM, ...
       'YLim',[0.35 n+0.75],'YDir','reverse', ...
       'XAxisLocation','top','YTick',[], ...
       'XTick',XTICKS, ...
       'TickDir','out','TickLength',[0.012 0.012], ...
       'XMinorTick','on','LineWidth',1.4, ...
       'FontName','Arial','FontSize',FS,'Color','none');
ax.YAxis.Visible = 'off';
ax.Clipping = 'off';

% tick labels with thousands separators
lbl = cell(1,numel(XTICKS));
for k = 1:numel(XTICKS)
    if XTICKS(k) >= 1000
        lbl{k} = regexprep(sprintf('%d',XTICKS(k)),'(\d)(?=(\d{3})+$)','$1,');
    else
        lbl{k} = sprintf('%d',XTICKS(k));
    end
end
set(ax,'XTickLabel',lbl);

xlabel(ax,'Distance (nm)','FontName','Arial', ...
       'FontSize',FS+1.5,'FontWeight','bold');

%% ---------------------- GUIDES ------------------------------------------
for g = GUIDES
    if g > XLIM(1) && g < XLIM(2)
        plot(ax,[g g],[0.55 n+0.45],'--','Color',[0.79 0.83 0.86], ...
             'LineWidth',0.9,'HandleVisibility','off');
        text(ax,g,n+0.62,sprintf('%g nm',g),'HorizontalAlignment','center', ...
             'VerticalAlignment','top','FontName','Arial', ...
             'FontSize',FS-1.5,'Color',[0.54 0.58 0.62]);
    end
end

%% ---------------------- GRADED BARS -------------------------------------
halfH = 0.22;             % bar half-height, data units

for i = 1:n
    mu = S(i).mu;  sd = S(i).sd;  col = S(i).col;

    tlo = mu - NSD*sd;
    thi = mu + NSD*sd;
    tlo = max(tlo,XLIM(1));  thi = min(thi,XLIM(2));   % clamp to axis

    % sample evenly across the bar
    edges = linspace(tlo,thi,NSEG+1);
    ctr   = 0.5*(edges(1:end-1) + edges(2:end));

    % Gaussian profile, then rescaled so the tips fade out no matter what
    % NSD is (a raw Gaussian still sits at 61% colour at +/-1 sd).
    wRaw  = exp(-0.5*((ctr-mu)./sd).^2);
    wEdge = exp(-0.5*NSD^2);             % raw weight at the nominal bar end
    w     = (wRaw - wEdge)./(1 - wEdge); % 0 at the tips, 1 at the mean
    w     = max(w,0).^GAMMA;             % shape the falloff
    w     = EDGEW + (1-EDGEW)*w;         % lift the tips off pure white

    for k = 1:NSEG
        c = 1 - w(k)*(1-col);            % blend white -> full colour
        patch(ax,'XData',[edges(k) edges(k+1) edges(k+1) edges(k)], ...
                  'YData',[i-halfH i-halfH i+halfH i+halfH], ...
                  'FaceColor',c,'EdgeColor','none','FaceAlpha',1);
    end

    % mean tick
    plot(ax,[mu mu],[i-halfH-0.06 i+halfH+0.06],'k-','LineWidth',2.6);

    % row label, left of the axis
    text(ax,XLIM(1)-0.025*xr,i,S(i).name,'HorizontalAlignment','right', ...
         'VerticalAlignment','middle','FontName','Arial', ...
         'FontSize',FS+3,'FontWeight','bold','Color',col);

    % value label, right of the bar
    if isempty(S(i).lab), vl = fmtDist(mu); else, vl = S(i).lab; end
    text(ax,thi+0.018*xr,i,vl,'HorizontalAlignment','left', ...
         'VerticalAlignment','middle','FontName','Arial', ...
         'FontSize',FS,'FontWeight','bold','Color',[0.17 0.21 0.25]);
end

%% ---------------------- SEPARATION ARROW --------------------------------
if SHOWSEP && n >= max(SEPROWS)
    a = S(SEPROWS(1)).mu;  b = S(SEPROWS(2)).mu;
    yA = mean(SEPROWS) - 0.08;
    plot(ax,[a b],[yA yA],'-','Color',[0.37 0.42 0.47],'LineWidth',1.2);
    plot(ax,a,yA,'<','MarkerSize',5,'MarkerFaceColor',[0.37 0.42 0.47], ...
         'MarkerEdgeColor','none');
    plot(ax,b,yA,'>','MarkerSize',5,'MarkerFaceColor',[0.37 0.42 0.47], ...
         'MarkerEdgeColor','none');
    if strcmpi(SEPMODE,'delta')
        sepTxt = sprintf('%s ~%.0f nm',Delta,abs(b-a));
    else
        sepTxt = sprintf('~%.1f%s separation',b/a,times);
    end
    text(ax,0.5*(a+b),yA-0.10,sepTxt, ...
         'HorizontalAlignment','center','VerticalAlignment','bottom', ...
         'FontName','Arial','FontSize',FS,'FontWeight','bold', ...
         'Color',[0.17 0.21 0.25]);
end

uistack(ax,'top');

%% ---------------------- EXPORT ------------------------------------------
if exist('exportgraphics','file')
    exportgraphics(fig,[OUTFILE '.pdf'],'ContentType','vector');
    exportgraphics(fig,[OUTFILE '.png'],'Resolution',600);
else
    print(fig,[OUTFILE '.pdf'],'-dpdf','-painters');
    print(fig,[OUTFILE '.png'],'-dpng','-r600');
end
fprintf('Wrote %s.pdf and %s.png\n',OUTFILE,OUTFILE);

%% ---------------------- HELPER ------------------------------------------
function s = fmtDist(d)
% Auto-format a distance in nm for the label to the right of each bar.
if d < 1000
    s = sprintf('~%.0f nm',d);
else
    s = sprintf('~%.2f %sm',d/1000,char(956));   % micrometres
end
end