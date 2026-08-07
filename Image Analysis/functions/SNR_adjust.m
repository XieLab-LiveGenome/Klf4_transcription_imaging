function out = SNR_adjust(I0, mul, mul2)
    % SNR_adjust: Adjusts image contrast based on fractional intensity thresholds.
    % mul:  Lower threshold fraction of max intensity (e.g., 0.2)
    % mul2: Middle threshold fraction of max intensity (e.g., 0.5)
    
    % Validate inputs
    if mul <= 0 || mul >= 0.5
        warning('SNR_adjust: mul should ideally be between 0 and 0.5 to keep up_thr > thr. Got %.2f', mul);
    end
    if mul2 <= mul || mul2 >= (1-mul)
        warning('SNR_adjust: mul2 should be between mul and (1-mul) for sequential binning.');
    end

    % Convert to double for math operations
    I2 = double(I0);
    max_val = max(I2(:));
    
    % Define thresholds
    thr     = mul * max_val;         % Lower bound
    med_thr = mul2 * max_val;        % Median boundary
    up_thr  = (1 - mul) * max_val;   % Upper bound
    
    % 1. Pre-calculate all masks based on the original un-modified image
    % Using >= to ensure no pixels fall through the cracks
    BG          =  I2<thr;
    DIM         = (I2 >= thr)     & (I2 < med_thr);
    MED         = (I2 >= med_thr) & (I2 < up_thr);
    VERY_BRIGHT = (I2 >= up_thr);


    
    % 2. Apply transformations safely
    I2(BG) = I2(BG)/2; 

    %====== DIM CONDENSATES (increase intensity) =============
    I2(DIM) = I2(DIM) * 2;       
    
    %===== MED CONDENSATES (no change in intensity) ===========
    I2(MED) = I2(MED);      
    
    %===== BRIGHT CONDENSATES (decrease intensity) ===========
    I2(VERY_BRIGHT) = I2(VERY_BRIGHT) / 2;   
    
    % 3. Apply median filter to smooth background noise
    I2 = medfilt2(I2);
    
    % 4. Clip limits and restore original data type
    I2 = max(0, I2);           % clip negative values
    I2 = min(65535, I2);       % clip overflow for 16-bit
    
    % Cast back to original image type (e.g., uint16) if applicable
    if isa(I0, 'uint16')
        out = uint16(I2);
    elseif isa(I0, 'uint8')
        % Just in case a 8-bit image is passed
        out = uint8(min(255, I2));
    else
        out = I2;
    end
end