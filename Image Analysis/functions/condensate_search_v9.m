function out = condensate_search_v9(stack_C, BW3, enh_pos, spacing, t_cond, mult, R_search_c, sz_cond, R2_thresh)
%%% Detects the closest condensate to an enhancer position.
%%% Strategy:
%%%   1) Try 2D Gaussian fitting + fit_Gaussian3D for 3D refinement
%%%   2) Compute R^2 of the 2D Gaussian fit
%%%   3) If R^2 < R2_thresh -> fall back to cntrd centroid + z intensity centroid
%%%   4) Return the closest condensate to the enhancer

% --- Default parameters ---
if nargin < 8 || isempty(sz_cond)
    sz_cond = 11;   % window size for cntrd (must be odd)
end
if mod(sz_cond, 2) == 0
    sz_cond = sz_cond + 1;
    warning('sz_cond must be odd. Adjusted to %d.', sz_cond);
end
if nargin < 9 || isempty(R2_thresh)
    R2_thresh = 0.2;
end

out    = cell(1, 1);
xpixel = spacing(1);
ypixel = spacing(2);
zpixel = spacing(3);
enh_x  = enh_pos(1);
enh_y  = enh_pos(2);
enh_z  = enh_pos(3);

C_cent_e_x = enh_x / xpixel;
C_cent_e_y = enh_y / ypixel;

max_iterations = 100;
iteration      = 1;
s_c            = 0;

% -----------------------------------------------------------------------
% 1. Iterative radius expansion until >= 1 condensate found
% -----------------------------------------------------------------------
while s_c < 1 && iteration <= max_iterations
    current_R_search = iteration * R_search_c;
    masked_img_cond  = imagemask(BW3, C_cent_e_x, C_cent_e_y, current_R_search);
    pk_cond          = pkfnd(masked_img_cond, mult * t_cond, 15);
    s_c              = size(pk_cond, 1);
    if s_c < 1
        iteration = iteration + 1;
    end
end

if s_c < 1
    out{1} = NaN;
    fprintf('No condensates found after %d iterations\n', max_iterations);
    return;
end

% -----------------------------------------------------------------------
% 2. Pre-compute cntrd for all spots (used as fallback)
% -----------------------------------------------------------------------
cntrd_out = cntrd(masked_img_cond, pk_cond, sz_cond);

num_spots_c = size(pk_cond, 1);
spotC_xyz   = zeros(num_spots_c, 3);
IC_fit      = zeros(num_spots_c, 1);
spot_size   = zeros(num_spots_c, 1);  % sigma (Gaussian) or sqrt(rg^2) (centroid)
spot_R2     = zeros(num_spots_c, 1);
method_used = cell(num_spots_c, 1);

bright_spots_x_c = pk_cond(:, 1);
bright_spots_y_c = pk_cond(:, 2);
R_c_3D_fit       = 8;

% -----------------------------------------------------------------------
% 3. Per-spot: Gaussian fit -> R^2 check -> route to Gaussian or centroid
% -----------------------------------------------------------------------
for i = 1:num_spots_c

    index1    = bright_spots_x_c(i);
    index2    = bright_spots_y_c(i);
    I_peak    = BW3(index2, index1);

    % --- Gaussian fit ---
    B_cond    = fitGaussian2(masked_img_cond, I_peak, index1, index2, 3, 6);
    % B_cond assumed: [amplitude, x0, y0, sigma_x, sigma_y, (offset)]

    % --- Compute R^2 of the 2D Gaussian fit ---
    R2 = compute_gaussian_R2(masked_img_cond, B_cond, index1, index2, 3);
    spot_R2(i) = R2;

    if R2 >= R2_thresh
        % ===============================================================
        % GOOD FIT: use Gaussian coordinates
        % ===============================================================
        method_used{i} = 'Gaussian';

        sigma_x       = B_cond(4);
        sigma_y       = B_cond(5);
        spot_size(i)  = sqrt(sigma_x * sigma_y);   % geometric mean of sigmas
        IC_fit(i)     = B_cond(1);

        guess_3D      = [B_cond(2), B_cond(3), enh_z / zpixel + 1];
        [x_ref_c, y_ref_c, z_ref_c, ~] = fit_Gaussian3D(stack_C, guess_3D, R_c_3D_fit, spacing);

        spotC_xyz(i, 1) = x_ref_c * xpixel;
        spotC_xyz(i, 2) = y_ref_c * ypixel;
        spotC_xyz(i, 3) = (z_ref_c - 1) * zpixel;

    else
        % ===============================================================
        % POOR FIT: fall back to intensity-weighted centroid (cntrd)
        % ===============================================================
        method_used{i} = 'Centroid';

        % Match this pkfnd peak to the cntrd output (nearest in xy)
        if ~isempty(cntrd_out)
            dists_to_pk   = sqrt((cntrd_out(:,1) - index1).^2 + ...
                                 (cntrd_out(:,2) - index2).^2);
            [~, ci]       = min(dists_to_pk);
            x_sub         = cntrd_out(ci, 1);
            y_sub         = cntrd_out(ci, 2);
            IC_fit(i)     = cntrd_out(ci, 3);
            spot_size(i)  = sqrt(cntrd_out(ci, 4));   % sqrt(rg^2)
        else
            % cntrd gave nothing; use raw pkfnd pixel position
            x_sub         = index1;
            y_sub         = index2;
            IC_fit(i)     = I_peak;
            spot_size(i)  = NaN;
        end

        % Intensity-weighted z centroid
        x_px      = max(1, min(size(stack_C, 2), round(x_sub)));
        y_px      = max(1, min(size(stack_C, 1), round(y_sub)));
        z_profile = double(squeeze(stack_C(y_px, x_px, :)));

        if sum(z_profile) > 0
            z_idx  = (1:numel(z_profile))';
            z_sub  = sum(z_idx .* z_profile) / sum(z_profile);
        else
            z_sub  = enh_z / zpixel + 1;
            warning('Flat z-profile at spot %d, falling back to enhancer z.', i);
        end

        spotC_xyz(i, 1) = x_sub  * xpixel;
        spotC_xyz(i, 2) = y_sub  * ypixel;
        spotC_xyz(i, 3) = (z_sub - 1) * zpixel;
    end
end

% -----------------------------------------------------------------------
% 4. Find the closest condensate to the enhancer
% -----------------------------------------------------------------------
EC_dist = sqrt((enh_x - spotC_xyz(:,1)).^2 + ...
               (enh_y - spotC_xyz(:,2)).^2 + ...
               (enh_z - spotC_xyz(:,3)).^2) * 1000;  % nm

[min_distance, closest_idx] = min(EC_dist);

fprintf('Closest condensate: %.1f nm away | R^2=%.3f | Method: %s\n', ...
        min_distance, spot_R2(closest_idx), method_used{closest_idx});

% Output: [x, y, z, spot_size, brightness, R2, distance_nm]
out{1} = [spotC_xyz(closest_idx, 1), ...
          spotC_xyz(closest_idx, 2), ...
          spotC_xyz(closest_idx, 3), ...
          spot_size(closest_idx),    ...
          IC_fit(closest_idx),       ...
          spot_R2(closest_idx),      ...
          min_distance];

end


% =======================================================================
% LOCAL HELPER: compute R^2 of a 2D Gaussian fit against the image patch
% =======================================================================
function R2 = compute_gaussian_R2(img, B_cond, x0_px, y0_px, half_win)
%%% Extracts the image patch used during fitting, reconstructs the Gaussian
%%% surface from B_cond, and computes R^2 = 1 - SS_res/SS_tot

try
    [nr, nc] = size(img);
    x_lo = max(1,  x0_px - half_win);
    x_hi = min(nc, x0_px + half_win);
    y_lo = max(1,  y0_px - half_win);
    y_hi = min(nr, y0_px + half_win);

    patch   = double(img(y_lo:y_hi, x_lo:x_hi));
    [yg, xg] = ndgrid(y_lo:y_hi, x_lo:x_hi);

    % Reconstruct Gaussian surface from fit parameters
    % B_cond: [amplitude, x0, y0, sigma_x, sigma_y, (optional offset)]
    amp     = B_cond(1);
    xc      = B_cond(2);
    yc      = B_cond(3);
    sx      = B_cond(4);
    sy      = B_cond(5);
    offset  = 0;
    if numel(B_cond) >= 6
        offset = B_cond(6);
    end

    fitted  = amp .* exp(-0.5 .* (((xg - xc)./sx).^2 + ((yg - yc)./sy).^2)) + offset;

    % R^2
    SS_res  = sum((patch(:) - fitted(:)).^2);
    SS_tot  = sum((patch(:) - mean(patch(:))).^2);

    if SS_tot == 0
        R2 = 0;   % flat patch — fitting is meaningless
    else
        R2 = 1 - SS_res / SS_tot;
    end

catch
    R2 = 0;   % if anything goes wrong, treat as poor fit
end

end
