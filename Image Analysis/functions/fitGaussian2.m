function [B_n, fit_info] = fitGaussian2(Zin, b_max, TSS_x, TSS_y, spot_s, mask_r)
% Fits a 2D Gaussian with offset to a spot image using only unmasked pixels.
%
% Inputs:
%   Zin    - input image (ROI, background-subtracted)
%   b_max  - initial guess for peak amplitude
%   TSS_x  - initial guess for x center (local coords)
%   TSS_y  - initial guess for y center (local coords)
%   spot_s - initial guess for sigma (pixels)
%   mask_r - radius of circular mask around (TSS_x, TSS_y) used for fitting
%
% Outputs:
%   B_n      = [amplitude, x_center, y_center, sigma_x, sigma_y, offset]
%   fit_info = struct with fields:
%       .exitflag  - lsqcurvefit exit flag (>0 = converged)
%       .resnorm   - sum of squared residuals
%       .rmse      - root-mean-square residual (per-pixel)
%       .n_pixels  - number of pixels used in fit
%       .at_bound  - logical vector, true where parameter hit a bound

    Zin = double(Zin);
    [m, n] = size(Zin);
    [X, Y] = meshgrid(1:n, 1:m);

    % Circular mask around the initial guess center
    M = (X - TSS_x).^2 + (Y - TSS_y).^2 <= mask_r^2;

    % Vectorize: fit only the pixels inside the mask
    xdata = [X(M), Y(M)];        % [Npix x 2]
    zdata = Zin(M);              % [Npix x 1]

    % 2D Gaussian with constant offset
    %   p = [A, cx, cy, sx, sy, bg]
    f = @(p, xy) p(1) * exp( -( (xy(:,1) - p(2)).^2 / (2*p(4)^2) + ...
                                (xy(:,2) - p(3)).^2 / (2*p(5)^2) ) ) + p(6);

    % Initial guesses
    bg0 = median(zdata);                        % robust offset estimate
    A0  = max(double(b_max) - bg0, eps);        % amplitude above background
    B0  = double([A0, TSS_x, TSS_y, spot_s, spot_s, bg0]);

    % Bounds
    %   - center constrained to the mask region
    %   - sigma bounded below by ~PSF lower limit (tune sigma_min for your system)
    %   - sigma bounded above by mask radius (larger is unphysical for a spot)
    sigma_min = 0.7;   % pixels; calibrate against bead data for your setup
    lb = [0,   TSS_x - mask_r, TSS_y - mask_r, sigma_min, sigma_min, -Inf];
    ub = [Inf, TSS_x + mask_r, TSS_y + mask_r, mask_r,    mask_r,    Inf];

    % Clip B0 into bounds (in case b_max is weird or spot_s > mask_r)
    B0 = min(max(B0, lb), ub);

    opts = optimoptions('lsqcurvefit', ...
        'Display', 'off', ...
        'FunctionTolerance', 1e-8, ...
        'StepTolerance', 1e-8);

    [B_n, resnorm, residual, exitflag] = ...
        lsqcurvefit(f, B0, xdata, zdata, lb, ub, opts);

    if nargout > 1
        npix = numel(zdata);
        tol  = 1e-6;
        fit_info = struct( ...
            'exitflag', exitflag, ...
            'resnorm',  resnorm, ...
            'rmse',     sqrt(resnorm / max(npix, 1)), ...
            'n_pixels', npix, ...
            'at_bound', (B_n - lb) < tol | (ub - B_n) < tol);
    end
end