function img_out = adjust_multichannel_image(img, thresholds, gammas, saturations)
% Adjust threshold, gamma, and saturation per channel.
% Works with uint8, uint16, or double images.

    % Determine type and max value
    if isa(img, 'uint8'),       max_val = 255;
    elseif isa(img, 'uint16'),  max_val = 65535;
    else,                       max_val = 1;
    end

    img_d = double(img);
    [~, ~, C] = size(img_d);

    % Expand scalars to match channel count
    thresholds  = repmat(thresholds(:),  C / numel(thresholds),  1);
    gammas      = repmat(gammas(:),      C / numel(gammas),      1);
    saturations = repmat(saturations(:), C / numel(saturations), 1);

    % Process each channel
    for c = 1:C
        ch = (img_d(:,:,c) - thresholds(c)) / (saturations(c) - thresholds(c));
        img_d(:,:,c) = max(0, min(1, ch)) .^ gammas(c) * max_val;
    end

    % Cast back to original type
    img_out = cast(img_d, class(img));
end