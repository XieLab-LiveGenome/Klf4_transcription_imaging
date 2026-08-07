%% ========================================================================
%  Timescale separation panel: E-P contact / BRD4 hub / burst lifetimes
%  Log time axis, one graded bar per species, deepest colour at the mean,
%  tips faded out toward white.
%  ========================================================================

clear;
clc

endash = char(8211);          % en dash for 'E-P'
times  = char(215);           % multiplication sign

%% ---------------------- USER INPUT --------------------------------------
% mu and sd in SECONDS. lab is the text printed to the right of each bar.

S(1).name = ['E' endash 'P contact'];
S(1).mu   = 7.1;                       % <-- mean lifetime (s)
S(1).sd   = 5;                         % <-- SD (s)
S(1).col  = [0.910 0.208 0.180];
S(1).lab  = '';                        % '' = auto-format from mu

S(2).name = 'BRD4 hub';
S(2).mu   = 790.2;                     % <-- mean lifetime (s)
S(2).sd   = 256.2;                     % <-- SD (s)
S(2).col  = [0.180 0.525 0.757];
S(2).lab  = '';

S(3).name = 'Burst';
S(3).mu   = 1632;                      % <-- mean lifetime (s)
S(3).sd   = 1009.8;                    % <-- SD (s)
S(3).col  = [0.180 0.620 0.357];
S(3).lab  = '';

%% ---------------------- OPTIONS -----------------------------------------
XLIM      = [1 1e4];      % axis range (s)
NSD       = 1.5;          % bar extends mu +/- NSD*sd
NSEG      = 300;          % gradient segments per bar (higher = smoother)
EDGEW     = 0.05;         % colour weight at the very ends of the bar
                          %   0 = pure white tips, 1 = no fade at all
GAMMA     = 1.25;         % >1 pushes more of the bar toward the light end
                          %   (1 = plain linear ramp from mean to tip)
DISTMODE  = 'linear';     % 'linear' : Gaussian in t      (mu, sd as entered)
                          % 'log'    : Gaussian in log10(t) (lognormal-like)
                          %   Use 'log' if sd >= mu, which is typical for
                          %   exponentially distributed lifetimes and makes
                          %   mu - NSD*sd go negative in linear mode.
SHOWSEP   = true;         % draw the fold-separation arrow between rows 1-2
SEPROWS   = [1 2];        % which two rows the arrow compares
MINUTEGUIDES = [60 600];  % dashed guides, labelled automatically
FS        = 10;           % base font size
FIGSIZE   = [18 7.2];     % cm, width x height
OUTFILE   = 'timescale_panel';   % written as .pdf and .png

n = numel(S);

%% ---------------------- FIGURE ------------------------------------------
fig = figure('Color','w','Units','centimeters', ...
             'Position',[2 2 FIGSIZE(1) FIGSIZE(2)]);
ax = axes('Parent',fig,'Units','normalized','Position',[0.20 0.13 0.63 0.62]);
hold(ax,'on'); box(ax,'off');

set(ax,'XScale','log','XLim',XLIM, ...
       'YLim',[0.35 n+0.75],'YDir','reverse', ...
       'XAxisLocation','top','YTick',[], ...
       'XTick',10.^(0:log10(XLIM(2))), ...
       'TickDir','out','TickLength',[0.012 0.012], ...
       'XMinorTick','on','LineWidth',1.4, ...
       'FontName','Arial','FontSize',FS,'Color','none');
ax.YAxis.Visible = 'off';
ax.Clipping = 'off';

% decade labels with thousands separators
dec = 10.^(0:log10(XLIM(2)));
lbl = cell(1,numel(dec));
for k = 1:numel(dec)
    if dec(k) >= 1000
        lbl{k} = regexprep(sprintf('%d',dec(k)),'(\d)(?=(\d{3})+$)','$1,');
    else
        lbl{k} = sprintf('%d',dec(k));
    end
end
set(ax,'XTickLabel',lbl);

xlabel(ax,'Mean lifetime (s)','FontName','Arial', ...
       'FontSize',FS+1.5,'FontWeight','bold');

%% ---------------------- MINUTE GUIDES -----------------------------------
for g = MINUTEGUIDES
    if g > XLIM(1) && g < XLIM(2)
        plot(ax,[g g],[0.55 n+0.45],'--','Color',[0.79 0.83 0.86], ...
             'LineWidth',0.9,'HandleVisibility','off');
        if g >= 60
            gl = sprintf('%g min',g/60);
        else
            gl = sprintf('%g s',g);
        end
        text(ax,g,n+0.62,gl,'HorizontalAlignment','center', ...
             'VerticalAlignment','top','FontName','Arial', ...
             'FontSize',FS-1.5,'Color',[0.54 0.58 0.62]);
    end
end

%% ---------------------- GRADED BARS -------------------------------------
halfH = 0.22;             % bar half-height, data units

for i = 1:n
    mu = S(i).mu;  sd = S(i).sd;  col = S(i).col;

    if strcmpi(DISTMODE,'log')
        % Gaussian in log10 space: sd converted to a multiplicative factor
        sLog = log10(1 + sd/mu);
        tlo  = 10^(log10(mu) - NSD*sLog);
        thi  = 10^(log10(mu) + NSD*sLog);
    else
        tlo = mu - NSD*sd;
        thi = mu + NSD*sd;
        if tlo <= 0
            tlo = mu/10;      % clamp: linear mode cannot reach <= 0 on a log axis
            warning(['Row %d: mu - %g*sd is <= 0, clamped to mu/10. ' ...
                     'Consider DISTMODE = ''log''.'],i,NSD);
        end
    end
    tlo = max(tlo,XLIM(1));  thi = min(thi,XLIM(2));

    % sample evenly on screen (i.e. evenly in log t)
    edges = logspace(log10(tlo),log10(thi),NSEG+1);
    ctr   = sqrt(edges(1:end-1).*edges(2:end));

    if strcmpi(DISTMODE,'log')
        wRaw = exp(-0.5*((log10(ctr)-log10(mu))./sLog).^2);
    else
        wRaw = exp(-0.5*((ctr-mu)./sd).^2);
    end

    % Rescale so the tips fade out regardless of NSD. A raw Gaussian is
    % still at 32% colour at +/-1.5 sd, which is why the old bars looked
    % blunt at the ends.
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
    text(ax,XLIM(1)*0.82,i,S(i).name,'HorizontalAlignment','right', ...
         'VerticalAlignment','middle','FontName','Arial', ...
         'FontSize',FS+3,'FontWeight','bold','Color',col);

    % value label, right of the bar
    if isempty(S(i).lab), vl = fmtTime(mu); else, vl = S(i).lab; end
    text(ax,thi*1.12,i,vl,'HorizontalAlignment','left', ...
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
    text(ax,sqrt(a*b),yA-0.10,sprintf('~%.0f%s separation',b/a,times), ...
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
function s = fmtTime(t)
% Auto-format a duration in seconds for the label to the right of each bar.
if t < 90
    s = sprintf('~%.0f s',t);
elseif t < 5400
    s = sprintf('~%.1f min',t/60);
else
    s = sprintf('~%.1f h',t/3600);
end
end