function out = condensate_search1(stack_C, BW3, enh_pos, spacing, t_cond, mult, R_search_c)

out = cell(1, 1); % Change to hold 3 closest condensates

xpixel = spacing(1);
ypixel = spacing(2);
zpixel = spacing(3);

enh_x = enh_pos(1);
enh_y = enh_pos(2);
enh_z = enh_pos(3);

C_cent_e_x = enh_x/xpixel;
C_cent_e_y = enh_y/ypixel;

masked_img_cond = imagemask(BW3, C_cent_e_x, C_cent_e_y, R_search_c); % BW3 is the cleared 2d projected image
pk_cond = pkfnd(masked_img_cond, mult * t_cond, 5); % -------------set this-------------------

s_c = size(pk_cond, 1); % find available bright pixels in that region

if s_c >= 1
    bright_spots_y_c = pk_cond(:, 2);
    bright_spots_x_c = pk_cond(:, 1);
    num_spots_c = size(pk_cond, 1); % number of bright pixels

    I_spots_c = zeros(num_spots_c, 1); % intensity at those bright spots
    fit_spot_c = zeros(num_spots_c, 1);
    fit_spot_ar = zeros(num_spots_c, 1);
    IC_fit = zeros(num_spots_c, 1);

    spotC_xyz = zeros(num_spots_c, 3); % Now will store 3D coordinates

    R_c_3D_fit = 6; % Radius for 3D Gaussian fitting

    for i = 1:num_spots_c
        index1 = bright_spots_x_c(i, 1);
        index2 = bright_spots_y_c(i, 1);
        I_spots_c(i, 1) = BW3(index2, index1);

        % Fit 2D Gaussian first for initial x, y, and intensity
        B_cond = fitGaussian2(masked_img_cond, I_spots_c(i, 1), index1, index2, 3); % Assuming spot_cond=3

        % Store 2D fit parameters
        fit_spot_c(i, 1) = (B_cond(4) * B_cond(5))^0.5;      % approximate spot size
        fit_spot_ar(i, 1) = max(B_cond(4), B_cond(5)) / min(B_cond(4), B_cond(5));    % aspect ratio calculations
        IC_fit(i, 1) = B_cond(1); % Fitted intensity

        % Now, perform 3D Gaussian fit to get accurate x, y, z
        % Initial guess for 3D fit using 2D results and enhancer Z
        % The enhancer Z (enh_z) is a good starting point for the Z-coordinate guess
        guess_3D = [B_cond(2), B_cond(3), enh_z/zpixel + 1]; % Convert enh_z to pixel units and add 1 for 1-based indexing

        [x_ref_c, y_ref_c, z_ref_c, ~] = fit_Gaussian3D(stack_C, guess_3D, R_c_3D_fit, spacing);

        spotC_xyz(i, 1) = x_ref_c * xpixel;
        spotC_xyz(i, 2) = y_ref_c * ypixel;
        spotC_xyz(i, 3) = (z_ref_c - 1) * zpixel; % Convert back to original z-units
    end

    EC_dist = zeros(num_spots_c, 1);
    for i = 1:num_spots_c
        % Calculate 3D Euclidean distance
        EC_dist(i, 1) = (((enh_x - spotC_xyz(i, 1))^2 + (enh_y - spotC_xyz(i, 2))^2 + (enh_z - spotC_xyz(i, 3))^2)^0.5) * 1000;
    end

    % Sort distances and get indices of the closest 3
    [~, sorted_indices] = sort(EC_dist, 'ascend');

        current_index = sorted_indices(1);

        Cond_x = spotC_xyz(current_index, 1);
        Cond_y = spotC_xyz(current_index, 2);
        Cond_z = spotC_xyz(current_index, 3);
        Cond_s = fit_spot_c(current_index, 1);
        Cond_ar = fit_spot_ar(current_index, 1);
        Cond_int = IC_fit(current_index, 1);

        % Store parameters for the current closest condensate
        out{1} = [Cond_x Cond_y Cond_z Cond_s Cond_int Cond_ar];
    
else % s_c == 0 (no condensates found)
    out{1} = NaN;
   
end
end