function imgOut = subtractBackground(img, sigma)
% subtractBackground  Gaussian illumination-based background subtraction.
%
%   imgOut = subtractBackground(img)
%   imgOut = subtractBackground(img, sigma)
%
%   Estimates the background as a heavily blurred copy of the image and
%   subtracts it.  Features much smaller than sigma (e.g. puncta) are
%   preserved; slow-varying gradients from uneven illumination, vignetting,
%   or autofluorescence are removed.
%
% INPUTS:
%   img   - 2D image (any numeric type; returned as double)
%   sigma - Gaussian blur sigma in pixels (default: 50)
%           Should be >> the diameter of the largest feature of interest.
%
% OUTPUT:
%   imgOut - Background-subtracted image (double, clipped to >= 0)

    if nargin < 2 || isempty(sigma)
        sigma = 50;
    end

    img    = double(img);
    bg     = imgaussfilt(img, sigma);
    imgOut = img - bg;
    imgOut(imgOut < 0) = 0;
end