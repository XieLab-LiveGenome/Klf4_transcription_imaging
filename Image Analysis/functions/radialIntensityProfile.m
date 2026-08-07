function [r_bins, radial_mean] = radialIntensityProfile(img, center, bin_size, max_radius)
% Computes the radial intensity profile around a given center
%
% INPUTS:
%   img        : 2D image (double or uint)
%   center     : [y0, x0] coordinates (row, column)
%   bin_size   : bin width in pixels (e.g., 1)
%   max_radius : maximum radius in pixels (e.g., 10)
%
% OUTPUTS:
%   r_bins     : radius bin centers
%   radial_mean: mean intensity at each radius

    % Image size
    [rows, cols] = size(img);

    % Create coordinate grids
    [xGrid, yGrid] = meshgrid(1:cols, 1:rows);

    % Compute distance of each pixel to the center
    r = sqrt((xGrid - center(2)).^2 + (yGrid - center(1)).^2);

    % Flatten arrays
    r_flat = r(:);
    I_flat = double(img(:));

    % Bin edges
    edges = 0:bin_size:max_radius;
    r_bins = edges(1:end-1) + bin_size/2;
    radial_mean = zeros(size(r_bins));

    % Bin and average
    for i = 1:length(r_bins)
        in_bin = r_flat >= edges(i) & r_flat < edges(i+1);
        radial_mean(i) = mean(I_flat(in_bin));
    end
end
