function B_n = fit2DGaussian(Z_local, b_max, TSS_x_local, TSS_y_local, spot_s)
% Z_local: Small (e.g., 21x21) image matrix, ONLY containing the spot.
% TSS_x_local/TSS_y_local: Local center guess (e.g., 11, 11).

m = size(Z_local, 1);
n = size(Z_local, 2);

% --- Define the 6-Parameter 2D Gaussian Model ---
% B_n(1): Amplitude
% B_n(2): X_center
% B_n(3): Y_center
% B_n(4): Sigma_x
% B_n(5): Sigma_y
% B_n(6): Background (CRITICAL FOR ACCURACY)
f_n = @(B_n, XY_n) (B_n(1) * exp( -((XY_n(:,:,1) - B_n(2)).^2 / (2 * B_n(4)^2) + ...
                                    (XY_n(:,:,2) - B_n(3)).^2 / (2 * B_n(5)^2) ) ) + B_n(6)); % <--- Added B_n(6)

% Create the coordinate grids for the small local image
[X_n, Y_n] = meshgrid(1:n, 1:m);
XY_n(:,:,1) = X_n;
XY_n(:,:,2) = Y_n;

% --- Set Initial Guess B0 (6 Parameters) ---
B0_n = zeros(1, 6);
B0_n(1) = b_max;             % Amplitude guess (Max pixel value)
B0_n(2) = TSS_x_local;       % X center guess (Local center, e.g., 11)
B0_n(3) = TSS_y_local;       % Y center guess (Local center, e.g., 11)
B0_n(4) = spot_s;            % Sigma X guess (e.g., 1-2 pixels)
B0_n(5) = spot_s;            % Sigma Y guess
B0_n(6) = min(Z_local(:));   % Background guess (Min pixel value in the ROI)

% Perform the non-linear least squares fit
B_n = lsqcurvefit(f_n, B0_n, XY_n, Z_local);
end