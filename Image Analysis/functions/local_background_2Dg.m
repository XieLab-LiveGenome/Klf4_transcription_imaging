function [background, sigma] = local_background_2Dg(img, center, max_radius, min_radius)
% Estimate local background by fitting a Gaussian + offset to the radial profile
%
% Model: I(r) = A * exp(-r^2 / (2*sigma^2)) + B
%
% INPUTS:
%   img        : 2D image (double)
%   center     : [y, x] center of the spot (row, col)
%   max_radius : outer radius for profile (pixels)
%   min_radius : inner radius to start fitting from (default = 2)
%
% OUTPUTS:
%   background : estimated background intensity (asymptotic offset B)
%   sigma      : fitted Gaussian sigma (pixels), useful for cross-checking spot size

    if nargin < 4
        min_radius = 2;
    end

    bin_size = 1;
    [r, I_r] = radialIntensityProfile(img, center, bin_size, max_radius);

    % Skip noisy center bins
    valid = r >= min_radius;
    r_fit = r(valid);
    I_fit = I_r(valid);

    if numel(r_fit) < 5
        warning('Too few bins for fitting. Falling back to outer median.');
        background = median(I_r(end-2:end));
        sigma = NaN;
        return;
    end

    % Model: I(r) = A * exp(-r^2 / (2*sigma^2)) + B
    model = @(p, r) p(1) * exp(-r.^2 / (2 * p(2)^2)) + p(3);

    % Initial guesses
    A0     = max(I_fit) - min(I_fit);
    sigma0 = max_radius / 4;
    B0     = min(I_fit);

    p0 = [A0, sigma0, B0];
    lb = [0,   0.5,   0];
    ub = [Inf, max_radius, max(I_fit)];

    opts = optimoptions('lsqcurvefit', 'Display', 'off');

    try
        p = lsqcurvefit(model, p0, r_fit, I_fit, lb, ub, opts);
        background = p(3);
        sigma = p(2);

        % Sanity check
        outer_med = median(I_r(end-2:end));
        if background > outer_med * 1.5 || background < 0
            background = outer_med;
        end
    catch
        warning('Gaussian fit failed. Falling back to outer median.');
        background = median(I_r(end-2:end));
        sigma = NaN;
    end
end