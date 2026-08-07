function out = maskavg(img,TSS_x,TSS_y,R_c)

m=size(img,1);
n=size(img,2);

% R_c=5; % constant radius 
% circular mask with spot centroid and radius for average spot intensity
[xx_n,yy_n] = meshgrid(1:n,1:m);
mask2_n = false(m,n);
for ii = 1:numel(R_c)
	mask2_n = mask2_n | hypot(xx_n - TSS_x, yy_n - TSS_y) <= R_c;
end

masked_img2_n= img().*mask2_n();
I1 = mean(nonzeros(masked_img2_n)); % quantification of average spot intensity

if isnan(I1)
    out=0;
else
    out=I1;

end