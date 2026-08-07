function out = point_location2D(label_matrix, x, y)
    % Ensure coordinates are integers (in case of sub-pixel coordinates)
    x = round(x);
    y = round(y);
    
    [rows, cols] = size(label_matrix);
    
    % 1. Boundary Check: Ensure the point is actually inside the image
    if x < 1 || x > cols || y < 1 || y > rows
        out = NaN;
        return;
    end
    
    % 2. Direct Lookup: MATLAB matrices are indexed as (row, column) -> (y, x)
    % This instantly grabs the ID without any loops.
    pixel_value = double(label_matrix(y, x));
    
    % 3. Determine the output
    if pixel_value > 0
        out = pixel_value; % This is the unique cell ID
    else
        out = NaN;         % Point landed on the background
    end
end