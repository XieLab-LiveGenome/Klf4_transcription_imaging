function out = track_spot2D(inputImage_MS2, C_cent_e1, C_cent_e2, R_search, m_e)
% Tracks the brightest spot, performs local cropping, fits Gaussian, and returns global coordinates.

PEAK_SEP = 3;       % Minimum separation for pkfnd (in pixels, ~FWHM)
FIT_R = 8;         % Radius for cropping and fitting (e.g., 10 pixels -> 21x21 ROI)

% --- Step 1: Initial Search in the full image ---
% Note: Using the whole image for the initial threshold check is often fine.
t_MS2 = m_e * max(inputImage_MS2(:)); 

% Using the original imagemask for the initial *search* region is okay,
% but for simplicity, we'll run pkfnd on the whole image (or a larger ROI)
% and then filter the results based on the initial guess (C_cent_e1/e2, R_search).

% pkfnd returns approximate global pixel coordinates (usually integer indices)
pk_all = pkfnd(inputImage_MS2, t_MS2, PEAK_SEP); 

% Filter peaks to those within the R_search region defined by the guess (C_cent_e1, C_cent_e2)
pk_MS2_idx = hypot(pk_all(:,1) - C_cent_e1, pk_all(:,2) - C_cent_e2) <= R_search;
pk_MS2 = pk_all(pk_MS2_idx, :);

s = size(pk_MS2, 1); % Number of candidate spots found

if s >= 1
    % --- Step 2: Select the brightest spot among candidates ---
    I_spots_n = zeros(s, 1);
    for i = 1:s
        % Sample the intensity at the peak pixel location
        I_spots_n(i) = inputImage_MS2(pk_MS2(i, 2), pk_MS2(i, 1));
    end
    
    [peak_intensity, tss_n] = max(I_spots_n);
    
    % Global integer coordinates of the best peak
    pk_x_best_global = pk_MS2(tss_n, 1);
    pk_y_best_global = pk_MS2(tss_n, 2);
    
    % --- Step 3: CRITICAL: Crop the local region of interest (ROI) ---
    half_size = FIT_R;

    % Calculate ROI boundaries (using round() ensures integer indices)
    r_min = max(1, round(pk_y_best_global) - half_size);
    r_max = min(size(inputImage_MS2, 1), round(pk_y_best_global) + half_size);
    c_min = max(1, round(pk_x_best_global) - half_size);
    c_max = min(size(inputImage_MS2, 2), round(pk_x_best_global) + half_size);

    % Extract the small cropped image for fitting
    Z_local = double(inputImage_MS2(r_min:r_max, c_min:c_max));
    
    % --- Step 4: Define Local Center for Fitting ---
    % The peak is now at a new center relative to the crop origin (r_min, c_min)
    Local_Center_X = round(pk_x_best_global) - c_min + 1;
    Local_Center_Y = round(pk_y_best_global) - r_min + 1;
    
    % Use a small sigma guess (e.g., 1.5 - 2.5 pixels)
    spot_s = 3; 

    % --- Step 5: Fit a 6-parameter 2D Gaussian to the local image ---
    % NOTE: Assumes fitGaussian2 is the FIXED 6-parameter version
    B_n = fit2DGaussian(Z_local, peak_intensity, Local_Center_X, Local_Center_Y, spot_s);

    % --- Step 6: Convert Local Fit back to Global Coordinates ---
    % Global Refined X = Global X of crop start + (Local Fitted X - 1)
    global_x_refined = c_min + (B_n(2) - 1); 
    global_y_refined = r_min + (B_n(3) - 1);

    out(1) = global_x_refined; % Precise Global X-coordinate (sub-pixel)
    out(2) = global_y_refined; % Precise Global Y-coordinate (sub-pixel)
    out(3) = s;
else 
    % No spot found
    out(1) = C_cent_e1;
    out(2) = C_cent_e2;
    out(3) = s;
end
end