img = inputImage_cond;
% I2 = adjust_multichannel_image(I,4000,2.9,6000);


% adjusted = adjust_multichannel_image(img, 500, 1.8, 55000);
% threshold=500 (black point), gamma=1.8 (brightened), saturation=55000 (white point raised)

% --- Multichannel (e.g. 3 channels) ---
thresholds  = [3000];  % black points per channel
gammas      = [0.3];  % gamma > 1 brightens midtones
saturations = [8000]; % white points (higher = more dynamic range used)

adjusted = adjust_multichannel_image(img, thresholds, gammas, saturations);

imshow(adjusted, []);  % [] auto-scales display for 16-bit

% imshowpair(I, I2, 'montage')
