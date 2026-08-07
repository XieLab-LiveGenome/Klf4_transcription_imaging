function out = imagemask(img,TSS_x,TSS_y,R_c)

m=size(img,1);
n=size(img,2);

% circular mask with spot centroid and radius
[xx_n,yy_n] = meshgrid(1:n,1:m);
mask1_n = false(m,n);
for ii = 1:numel(R_c)
	mask1_n = mask1_n | hypot(xx_n - TSS_x, yy_n - TSS_y) <= R_c;
end

I=double(img);

out=I().*mask1_n();