function mask = CircularMask3D(m, n, z, center, radius, z_range)
    % m, n, z - dimensions of the full stack
    % center - [x0, y0, z0]
    % radius - circle radius in xy-plane
    % z_range - number of slices above and below z0 (can be scalar or [low, high])
    
    % Default to symmetric range
    if isscalar(z_range)
        z_range = [-z_range, z_range];
    end
    
    x0 = center(1);
    y0 = center(2);
    z0 = center(3);

    % Create output mask
    mask = false(m, n, z);

    % Compute valid Z slices
    z_start = max(1, round(z0 + z_range(1)));
    z_end   = min(z, round(z0 + z_range(2)));

    % Prepare XY grid
    [X, Y] = meshgrid(1:n, 1:m);
    dist2 = (X - x0).^2 + (Y - y0).^2;
    circle_mask = dist2 <= radius^2;

    % Insert circle into each Z-slice in the z_range
    for k = z_start:z_end
        mask(:,:,k) = circle_mask;
    end
end
