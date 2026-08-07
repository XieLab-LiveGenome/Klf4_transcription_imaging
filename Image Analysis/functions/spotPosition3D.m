function [out , mip_c] = spotPosition3D(MIP, stack_3D, C_cent_e1, C_cent_e2, slices, R_fit, m_enh, spacing)
    % Initialize output with NaNs to prevent crashes if the spot is lost
    out = [NaN, NaN, NaN];
    mip_c = [NaN, NaN];

    xpixel = spacing(1);
    ypixel = spacing(2);
    zpixel = spacing(3);

    MIP2 = SNR_inc2(MIP,m_enh);    % boost SNR of spots

    half_size = R_fit;

    % Calculate ROI boundaries (using round() ensures integer indices)
    r_min = max(1, round(C_cent_e2) - half_size);
    r_max = min(size(MIP, 1), round(C_cent_e2) + half_size);
    c_min = max(1, round(C_cent_e1) - half_size);
    c_max = min(size(MIP, 2), round(C_cent_e1) + half_size);

    % Extract the small cropped image for fitting
    Z_local = double(MIP2(r_min:r_max, c_min:c_max));

    peak_intensity = max(Z_local(:));
    
    % --- Step 4: Define Local Center for Fitting ---
    % The peak is now at a new center relative to the crop origin (r_min, c_min)
    Local_Center_X = round(C_cent_e1) - c_min + 1;
    Local_Center_Y = round(C_cent_e2) - r_min + 1;
    
    % Use a small sigma guess (e.g., 1.5 - 2.5 pixels)
    spot_s = 3; 

    % --- Step 5: Fit a 6-parameter 2D Gaussian to the local image ---
    % NOTE: Assumes fitGaussian2 is the FIXED 6-parameter version
    B_n = fit2DGaussian(Z_local, peak_intensity, Local_Center_X, Local_Center_Y, spot_s);

    % --- Step 6: Convert Local Fit back to Global Coordinates ---
    % Global Refined X = Global X of crop start + (Local Fitted X - 1)
    global_x_refined = c_min + (B_n(2) - 1); 
    global_y_refined = r_min + (B_n(3) - 1);

    C_cent_e1_g = round(global_x_refined);
    C_cent_e2_g = round(global_y_refined);
    % C3 = round(slices/2);        % assume middle of stack for Z guess
    
    I_enh_z = zeros(slices,1);

    for p=1:slices
        Im_enh=stack_3D(:,:,p);
        r= C_cent_e2_g;
        c= C_cent_e1_g;
        I_enh_z(p,1) = maskavg (Im_enh,c,r,3);    % this is correct order of c and r verified
    end

    z_corr = z_loc(I_enh_z(:,1),slices);
    C3=round(z_corr);

        
        % Guess latest 2D coordinates and 3D coordinate from last frame
        guess_e = [C_cent_e1_g, C_cent_e2_g, C3];     
        
        [x_ref_e, y_ref_e, z_ref_e, params_e, R2_e] = fit_Gaussian3D(stack_3D, guess_e, R_fit, spacing);
        
        % Check if 3D Gaussian fit was successful or if it crashed/returned NaN
        if ~isnan(x_ref_e) 
            out(1) = x_ref_e * xpixel;
            out(2) = y_ref_e * ypixel;
            out(3) = (z_ref_e - 1) * zpixel; 
        else
            % Fallback to 2D coordinates and middle slice Z guess
            % Optional: Log or display a warning when the spot isn't found
            warning('3D Gaussian fit failed. Falling back to 2D + midplane Z.');
            out(1) = C_cent_e1_g * xpixel;
            out(2) = C_cent_e2_g * ypixel;
            out(3) = (C3 - 1) * zpixel; % Applied -1 for consistency with successful fit
        end
mip_c(1) = global_x_refined;
mip_c(2) = global_y_refined;

        
end