function stack = single_zstack(r, scene, time, channel)
    % r: 6D array from readimage6D
    % scene, time, channel: indices to extract
    % output: stack of shape (Y, X, Z)
    
    % Extract slice
    stack6D = squeeze(r(scene, time, :, channel, :, :));
    
    % Permute to (Y, X, Z)
    stack = permute(stack6D, [2 3 1]);
end
