function out = hist_gaussian(data, x_min, x_max, colr)

h = histogram(data, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'LineWidth', 1.0, 'BinWidth', 100, 'EdgeColor', 'none');


xlim([x_min x_max]); 

hold on

% Fit a Gaussian (normal distribution)
mu = mean(data);
sigma = std(data);

% Generate smooth x values across the range
x_fit = linspace(x_min, x_max, 200);
y_fit = normpdf(x_fit, mu, sigma);  % PDF of normal distribution
y_fit = y_fit * (h.BinWidth);       % scale to match histogram 'probability' mode

% Plot Gaussian fit
out = plot(x_fit, y_fit, '-','Color',colr,'LineWidth', 2);

end