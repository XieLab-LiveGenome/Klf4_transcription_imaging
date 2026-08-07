function [maxFeret, minFeret] = feretDiameters(bw)
    B = bwboundaries(bw);
    if isempty(B)
        maxFeret = NaN;
        minFeret = NaN;
        return;
    end
    boundary = B{1};
    y = boundary(:,1);
    x = boundary(:,2);

    % Max Feret: largest pairwise boundary distance
    D = pdist2([x y], [x y]);
    maxFeret = max(D(:));

    % Min Feret: smallest caliper width across all orientations
    angles = linspace(0, pi, 180);
    minFeret = Inf;
    for theta = angles
        proj = x * cos(theta) + y * sin(theta);
        minFeret = min(minFeret, max(proj) - min(proj));
    end
end