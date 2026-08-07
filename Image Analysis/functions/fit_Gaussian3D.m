function [x_fit, y_fit, z_fit, fit_params, R2] = fit_Gaussian3D(stack, guess, radius, spacing)
%REFINEGAUSSIAN3D Refines XYZ coordinates in a z-stack using 3D Gaussian fitting
%
% Inputs:
%   stack    - 3D array (Y x X x Z) from a CZI image stack
%   guess    - [x0, y0, z0] approximate coordinates
%   radius   - ROI radius around the guess (default = 5)
%   spacing  - [sx, sy, sz] voxel size in physical units (optional)
% Outputs:
%   x_fit, y_fit, z_fit - refined coordinates (NaN if fit fails)
%   fit_params - [A, x0, y0, z0, sigma_x, sigma_y, sigma_z, offset] (NaN if fit fails)
%   R2 - goodness of fit (NaN if fit fails)

% Initialize outputs as NaN (in case of failure)
x_fit = NaN;
y_fit = NaN;
z_fit = NaN;
fit_params = NaN(1, 8);
R2 = NaN;

if nargin < 3 || isempty(radius)
    radius = 5;
end
if nargin < 4 || isempty(spacing)
    spacing = [1, 1, 1];  % Assume isotropic voxels
end

% Validate inputs
if isempty(stack) || any(isnan(guess)) || any(isinf(guess))
    warning('fit_Gaussian3D: Invalid input data, skipping frame.');
    return;
end

try
    x0 = guess(1); 
    y0 = guess(2); 
    z0 = guess(3);
    [Ymax, Xmax, Zmax] = size(stack);
    
    % Validate guess is within bounds
    if x0 < 1 || x0 > Xmax || y0 < 1 || y0 > Ymax 
        warning('fit_Gaussian3D: Guess coordinates out of bounds, skipping frame.');
        return;
    end
    
    % Define subvolume bounds
    xRange = max(1, round(x0-radius)) : min(Xmax, round(x0+radius));
    yRange = max(1, round(y0-radius)) : min(Ymax, round(y0+radius));
    zRange = 1 : Zmax;
    
    % Check if ranges are valid
    if isempty(xRange) || isempty(yRange) || isempty(zRange)
        warning('fit_Gaussian3D: Invalid ROI range, skipping frame.');
        return;
    end
    
    % Extract subvolume
    subvol = double(stack(yRange, xRange, zRange));
    
    % Check for valid data
    if all(isnan(subvol(:))) || all(subvol(:) == 0)
        warning('fit_Gaussian3D: Subvolume contains no valid data, skipping frame.');
        return;
    end
    
    % Coordinate grid with scaling
    [Yg, Xg, Zg] = ndgrid(yRange, xRange, zRange);
    coords_scaled = [Xg(:)*spacing(1), Yg(:)*spacing(2), Zg(:)*spacing(3)];
    intensities = subvol(:);
    
    % Remove NaN and Inf values
    valid_idx = ~isnan(intensities) & ~isinf(intensities);
    if sum(valid_idx) < 10  % Need at least 10 points for meaningful fit
        warning('fit_Gaussian3D: Insufficient valid data points, skipping frame.');
        return;
    end
    intensities = intensities(valid_idx);
    coords_scaled = coords_scaled(valid_idx, :);
    
    % Initial guess: [A, x0, y0, z0, sx, sy, sz, B]
    A0 = max(intensities);
    B0 = min(intensities);
    
    % Validate A0 and B0
    if A0 <= B0 || isnan(A0) || isnan(B0)
        warning('fit_Gaussian3D: Invalid intensity range, skipping frame.');
        return;
    end
    
    p0 = [A0, x0*spacing(1), y0*spacing(2), z0*spacing(3), ...
          2*spacing(1), 2*spacing(2), 1*spacing(3), B0];
    
    % Gaussian model
    gaussian3D = @(p, coords) ...
        p(1) * exp( -((coords(:,1)-p(2)).^2/(2*p(5)^2) + ...
                      (coords(:,2)-p(3)).^2/(2*p(6)^2) + ...
                      (coords(:,3)-p(4)).^2/(2*p(7)^2)) ) + p(8);
    
    % Define model as a function of p (CORRECTED - now uses coords parameter)
    modelFun = @(p, coords) gaussian3D(p, coords);
    
    f = 1.0;
    % The allowed range for x_fit (in physical units)
    x_lower_bound_phys = (x0 - f*radius) * spacing(1);
    x_upper_bound_phys = (x0 + f*radius) * spacing(1);
    % The allowed range for y_fit (in physical units)
    y_lower_bound_phys = (y0 - f*radius) * spacing(2);
    y_upper_bound_phys = (y0 + f*radius) * spacing(2);
    
    lb = [0, x_lower_bound_phys, y_lower_bound_phys, 0, 0.1, 0.1, 0.1, 0];
    ub = [Inf, x_upper_bound_phys, y_upper_bound_phys, Zmax*spacing(3), Inf, Inf, Inf, max(intensities)];
    
    % Run lsqcurvefit with error handling
    opts = optimoptions('lsqcurvefit', 'Display', 'off', 'MaxIterations', 400);
    
    % Test initial function evaluation
    try
        test_val = modelFun(p0, coords_scaled);
        if any(isnan(test_val)) || any(isinf(test_val))
            warning('fit_Gaussian3D: Initial function evaluation produces invalid values, skipping frame.');
            return;
        end
    catch ME
        warning('fit_Gaussian3D: Initial function evaluation failed: %s');
        return;
    end
    
    % Perform the fit
    fit_params = lsqcurvefit(modelFun, p0, coords_scaled, intensities, lb, ub, opts);
    
    % Convert refined coords back to voxel indices
    x_fit = fit_params(2) / spacing(1);
    y_fit = fit_params(3) / spacing(2);
    z_fit = fit_params(4) / spacing(3);
    
    % Validate fitted parameters
    if any(isnan(fit_params)) || any(isinf(fit_params))
        warning('fit_Gaussian3D: Fit produced invalid parameters, skipping frame.');
        x_fit = NaN; y_fit = NaN; z_fit = NaN;
        fit_params = NaN(1, 8);
        R2 = NaN;
        return;
    end
    
    % --- Goodness of fit (R²) ---
    fit_vals = modelFun(fit_params, coords_scaled);
    residuals = intensities - fit_vals;
    SS_res = sum(residuals.^2);
    SS_tot = sum((intensities - mean(intensities)).^2);
    
    % Standard R2
    if SS_tot > 0
        R2 = 1 - SS_res/SS_tot;
    else
        R2 = NaN;
    end
    
catch ME
    % Catch any other errors and return NaN values
    warning('fit_Gaussian3D: Fitting failed with error: %s. Skipping frame.');
    x_fit = NaN;
    y_fit = NaN;
    z_fit = NaN;
    fit_params = NaN(1, 8);
    R2 = NaN;
end

end