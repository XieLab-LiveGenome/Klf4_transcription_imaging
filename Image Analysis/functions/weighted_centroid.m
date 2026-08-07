% Create synthetic 3D image
% V = rand(50, 50, 50); % intensity volume

V =  extractZstack(Input_zstack, channel_cond, 54);

% Threshold and label
% bw = V > 0.95;
% L = bwlabeln(bw);  % label connected components
center = [398,137,9];
L = CircularMask3D(785,711,18,center,100,18);

% Use intensity as weight for centroid
stats = regionprops3(L, V, 'WeightedCentroid');

% Display weighted centroids
centroids = stats.WeightedCentroid;
disp(centroids);

com_centroids = spot_COM_centroid(V,398,137,100);
disp(com_centroids);