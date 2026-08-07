function [p_fit, binCenters, counts, hLine, stats] = fitGaussian_histogram(data, xmin, xmax, binWidth, color, linestyle, linewidth, varargin)
%FITGAUSSIAN_HISTOGRAM Fit a Gaussian to a histogram of DATA and overlay it.
%
%   [p_fit, binCenters, counts] = fitGaussian_histogram(data, xmin, xmax, binWidth, color)
%   [...] = fitGaussian_histogram(data, xmin, xmax, binWidth, color, linestyle, linewidth)
%
% Required inputs:
%   data      - vector of data (NaN/Inf are removed)
%   xmin      - lower edge of histogram
%   xmax      - upper edge of histogram
%   binWidth  - bin width
%   color     - line color ('r' or [1 0 0])
%   linestyle - line style ('-','--',':','-.')  [optional, default '-']
%   linewidth - curve line width                [optional, default 2]
%
% Name/value options:
%   'Normalization'  'pdf' (default) | 'probability' | 'count'
%                    'pdf' makes the fitted amplitude independent of binWidth.
%   'Method'         'lsq' (default) least-squares fit to bin heights
%                    'mle'           mean/std of the in-range data (unbiased,
%                                    amplitude set analytically to match the
%                                    chosen normalization)
%   'LineWidth'      scalar (default 2)
%   'Axes'           axes handle (default gca)
%   'PlotHistogram'  true|false (default false) - draw the bars too
%
% Outputs:
%   p_fit       [amplitude, mu, sigma]  (sigma always positive)
%   binCenters  bin centers
%   counts      bin heights under the chosen normalization
%   hLine       handle to the plotted curve
%   stats       struct with R2, RMSE, N, fracExcluded, area
%
% Notes:
%   - Least-squares fitting to bin heights is sensitive to binning and is
%     slightly biased. Use 'Method','mle' when you intend to *report* mu/sigma;
%     'lsq' is fine when the curve is only there to guide the eye.
%   - Requires Optimization Toolbox for 'lsq'; falls back to fminsearch if
%     lsqcurvefit is unavailable.

% ---------- arguments ----------
if nargin < 6 || isempty(linestyle), linestyle = '-'; end

% LINEWIDTH may be passed positionally (7th arg) or as a 'LineWidth' name/value
% pair. If the 7th arg is text, it's really the start of the name/value list.
if nargin < 7 || isempty(linewidth)
    linewidth = 2;
elseif ischar(linewidth) || isstring(linewidth)
    varargin  = [{linewidth}, varargin];
    linewidth = 2;
elseif ~(isnumeric(linewidth) && isscalar(linewidth) && linewidth > 0)
    error('fitGaussian_histogram:BadLineWidth', ...
          'linewidth must be a positive scalar.');
end

ip = inputParser;
ip.addParameter('Normalization', 'pdf', @(s) ischar(s) || isstring(s));
ip.addParameter('Method', 'lsq', @(s) ischar(s) || isstring(s));
ip.addParameter('LineWidth', linewidth, @(x) isnumeric(x) && isscalar(x) && x > 0);
ip.addParameter('Axes', [], @(h) isempty(h) || isgraphics(h, 'axes'));
ip.addParameter('PlotHistogram', false, @(x) islogical(x) && isscalar(x));
ip.parse(varargin{:});
opt = ip.Results;

normMode = lower(char(opt.Normalization));
method   = lower(char(opt.Method));
ax = opt.Axes;
if isempty(ax), ax = gca; end

% ---------- clean data ----------
data = data(:);
data = data(isfinite(data));
if numel(data) < 3
    error('fitGaussian_histogram:TooFewPoints', ...
          'Need at least 3 finite data points (got %d).', numel(data));
end

% ---------- edges ----------
edges = xmin:binWidth:xmax;
if edges(end) < xmax - 1e-12          % make sure xmax is actually covered
    edges(end+1) = edges(end) + binWidth;
end

inRange = data >= edges(1) & data <= edges(end);
fracOut = 1 - mean(inRange);
if fracOut > 0.01
    warning('fitGaussian_histogram:DataOutsideRange', ...
        '%.1f%% of the data lies outside [%g %g] and is excluded from the fit.', ...
        100*fracOut, edges(1), edges(end));
end
d = data(inRange);
if numel(d) < 3
    error('fitGaussian_histogram:EmptyRange', 'No data inside [%g %g].', edges(1), edges(end));
end

% ---------- histogram ----------
counts     = histcounts(d, edges, 'Normalization', normMode);
binCenters = edges(1:end-1) + diff(edges)/2;

gaussEq = @(p, x) p(1) .* exp(-((x - p(2)).^2) ./ (2 * p(3).^2));

% amplitude that makes a Gaussian of width sd integrate correctly for this
% normalization
ampFor = @(sd) ampForNorm(sd, normMode, binWidth, numel(d));

% ---------- fit ----------
switch method
    case 'mle'
        mu = mean(d);
        sd = std(d);
        p_fit = [ampFor(sd), mu, sd];

    case 'lsq'
        % moment-based initial guess taken from the BINNED data, so points
        % outside the range cannot pull it around
        w   = counts(:).' / sum(counts);
        mu0 = sum(w .* binCenters);
        sd0 = sqrt(sum(w .* (binCenters - mu0).^2));
        if ~isfinite(sd0) || sd0 <= 0, sd0 = binWidth; end

        init = [max(counts), mu0, sd0];
        lb   = [0,    edges(1),   eps];
        ub   = [Inf,  edges(end), Inf];

        if exist('lsqcurvefit', 'file') == 2
            o = optimset('Display', 'off');
            p_fit = lsqcurvefit(gaussEq, init, binCenters, counts, lb, ub, o);
        else
            o   = optimset('Display','off','TolX',1e-10,'TolFun',1e-10,'MaxFunEvals',1e4);
            obj = @(p) sum((gaussEq([abs(p(1)) p(2) abs(p(3))], binCenters) - counts).^2);
            p_fit    = fminsearch(obj, init, o);
            p_fit(1) = abs(p_fit(1));
        end

    otherwise
        error('fitGaussian_histogram:BadMethod', 'Method must be ''lsq'' or ''mle''.');
end
p_fit(3) = abs(p_fit(3));   % sigma enters only as sigma^2; force a positive report

% ---------- goodness of fit ----------
resid = counts - gaussEq(p_fit, binCenters);
ssRes = sum(resid.^2);
ssTot = sum((counts - mean(counts)).^2);
stats = struct( ...
    'R2',           1 - ssRes/ssTot, ...
    'RMSE',         sqrt(mean(resid.^2)), ...
    'N',            numel(d), ...
    'fracExcluded', fracOut, ...
    'area',         p_fit(1) * p_fit(3) * sqrt(2*pi));   % should be ~1 for 'pdf'

% ---------- plot ----------
wasHeld = ishold(ax);
hold(ax, 'on');

if opt.PlotHistogram
    bar(ax, binCenters, counts, 1, 'FaceColor', 'none', ...
        'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 1);
end

x_fit = linspace(edges(1), edges(end), 400);
hLine = plot(ax, x_fit, gaussEq(p_fit, x_fit), ...
             'Color', color, 'LineStyle', linestyle, 'LineWidth', opt.LineWidth);

if ~wasHeld, hold(ax, 'off'); end

end

% =====================================================================
function A = ampForNorm(sd, normMode, binWidth, n)
% Peak height of a Gaussian with std SD under the given histcounts normalization
switch normMode
    case 'pdf'
        A = 1 / (sd * sqrt(2*pi));
    case 'probability'
        A = binWidth / (sd * sqrt(2*pi));
    case 'count'
        A = n * binWidth / (sd * sqrt(2*pi));
    otherwise
        error('fitGaussian_histogram:BadNormalization', ...
              'Normalization must be ''pdf'', ''probability'', or ''count''.');
end
end