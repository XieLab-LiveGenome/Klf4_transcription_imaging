function out = radial_spot_size2D(img_sub, center, max_radius, pixel_s)
% FWHM from radial intensity profile with interpolation
%
% INPUTS:
%   img_sub    : 2D background-subtracted image
%   center     : [y, x] center (row, col)
%   max_radius : max radius for profile (pixels)
%   pixel_s    : pixel size (um)
%
% OUTPUT:
%   out        : FWHM diameter in nm

    bin_size = 1;
    [r_s, I_r_s] = radialIntensityProfile(img_sub, center, bin_size, max_radius);

    % Guard against empty or flat profiles
    if isempty(I_r_s) || all(I_r_s == 0)
        out = NaN;
        return;
    end

    I0 = max(I_r_s);
    half_max = I0 / 2;

    % Find first bin that drops below half-max
    idx = find(I_r_s <= half_max, 1);

    % Handle edge cases
    if isempty(idx)
        % Profile never drops below half-max within search radius
        out = NaN;
        return;
    end

    if idx == 1
        % Peak is not at center — use the peak bin as reference instead
        [~, peak_idx] = max(I_r_s);
        idx = find(I_r_s(peak_idx:end) <= half_max, 1);
        if isempty(idx)
            out = NaN;
            return;
        end
        idx = idx + peak_idx - 1;
        if idx < 2
            out = NaN;
            return;
        end
    end

    % Linear interpolation between bins flanking half-max
    r1 = r_s(idx - 1);
    r2 = r_s(idx);
    I1 = I_r_s(idx - 1);
    I2 = I_r_s(idx);

    if I1 == I2
        r_half = r1;
    else
        r_half = r1 + (half_max - I1) * (r2 - r1) / (I2 - I1);
    end

    fwhm_est = 2 * r_half;
    out = fwhm_est * pixel_s * 1000;  % nm
end