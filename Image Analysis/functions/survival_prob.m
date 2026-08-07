function out = survival_prob(col)

% Sort the data

sortedCol = sort(col);
n = length(sortedCol);

% Survival probability: S(x) = P(X > x)
survivalProb = (n:-1:1)' / n;

% --- SMOOTHING SECTION ---
    % Generate a smooth evaluation grid
    gridPoints = linspace(min(col), max(col), 500);
    
    % Compute the smooth CDF using Kernel Density Estimation
    % 'Function','cdf' gives the smooth P(X <= x)
    [f, xi] = ksdensity(col, gridPoints, 'Function', 'cdf');
    
    % Convert CDF to Survival Probability: S(x) = 1 - CDF
    survivalProbSmooth = 1 - f;
    
    % Return both for comparison
    out{1} = survivalProb;       % Original y
    out{2} = sortedCol;          % Original x
    out{3} = survivalProbSmooth; % Smoothed y
    out{4} = xi;                 % Smoothed x

end
