function out = SNR_inc2(I0, mul)
    % mul: threshold fraction of max intensity (0 < mul < 1)
    % Pixels above mul*max are brightened; below are dimmed

    if mul <= 0 || mul >= 1
        warning('SNR_inc2: mul should be between 0 and 1 (e.g. 0.5). Got %.2f', mul);
    end

    I2     = double(I0);
    thr    = mul * max(I2(:));  % threshold as fraction of max intensity

    bright = I2 > thr;
    I2( bright) = I2( bright) * 2;   % boost bright pixels
    I2(~bright) = I2(~bright) / 2;   % suppress dim pixels

    % Median filter applied after contrast enhancement to smooth dim background
    I2 = medfilt2(I2);
    I2 = medfilt2(I2);

    I2 = max(0,     I2);        % clip negative values
    I2 = min(65535, I2);        % clip overflow

    out = I2;
end