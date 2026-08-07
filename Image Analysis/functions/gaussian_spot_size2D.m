function out = gaussian_spot_size2D(img, x, y, pixel_size, R_search_e)
% Crops an ROI around (x,y), then fits a 2D Gaussian.
% Inputs:
%   img         - 2D image (background-subtracted)
%   x, y        - rough centroid (x = column, y = row)
%   pixel_size  - physical pixel size (microns)
%   R_search_e  - half-width of ROI in pixels
% Outputs:
%   out(1) - spot size (nm), geometric mean of sigmas * 4
%   out(2) - aspect ratio
%   out(3) - fitted x center (global coords)
%   out(4) - fitted y center (global coords)

    m = size(img, 1);  % rows
    n = size(img, 2);  % cols

    % ---- Crop ROI around rough centroid ----
    x_min = max(1, round(x) - R_search_e);
    x_max = min(n, round(x) + R_search_e);
    y_min = max(1, round(y) - R_search_e);
    y_max = min(m, round(y) + R_search_e);

    ROI = double(img(y_min:y_max, x_min:x_max));

    % ---- Local coordinates of the rough centroid inside the ROI ----
    x_local = x - x_min + 1;
    y_local = y - y_min + 1;

    % ---- Apply circular mask within ROI ----
    [xx, yy] = meshgrid(1:size(ROI,2), 1:size(ROI,1));
    mask = hypot(xx - x_local, yy - y_local) <= R_search_e;
    ROI = ROI .* mask;

    % ---- 2D Gaussian fit on the ROI ----
    spot_s = 3;  % initial sigma guess
    b_max = max(ROI(:));
    B_e = fitGaussian2(ROI, b_max, x_local, y_local, spot_s,R_search_e);

    % ---- Convert fitted center back to global image coordinates ----
    x_global = B_e(2) + x_min - 1;
    y_global = B_e(3) + y_min - 1;

    % ---- Output ----
    out(1) = 2.355 * sqrt(B_e(4) * B_e(5)) * pixel_size * 1000;  % spot size in nm
    out(2) = max(B_e(4), B_e(5)) / min(B_e(4), B_e(5));      % aspect ratio
    out(3) = x_global;
    out(4) = y_global;
end