function out = point_location(BW, x, y)

CC = bwconncomp(BW);

index = NaN;
for i = 1:CC.NumObjects
    if ismember(sub2ind(size(BW), y, x), CC.PixelIdxList{i})
        index = i;
        break;
    end
end

out = index;

end