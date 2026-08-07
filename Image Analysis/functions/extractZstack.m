function zStack = extractZstack(filename, channelIndex, timeIndex)
%EXTRACTZSTACKCZI Extract a single-channel z-stack from a 16-bit CZI file
%
% zStack = extractZStackCZI(filename, channelIndex, timeIndex)
%
% Inputs:
%   filename      - Path to the .czi file
%   channelIndex  - Channel index (1-based, e.g. 1 = first channel)
%   timeIndex     - Timepoint index (1-based, e.g. 1 = first timepoint)
%
% Output:
%   zStack        - 3D uint16 matrix (Y x X x Z) z-stack

    % Open file using Bio-Formats
    reader = bfGetReader(filename);

    % Get dimensions
    sizeX = reader.getSizeX();
    sizeY = reader.getSizeY();
    sizeZ = reader.getSizeZ();
    sizeC = reader.getSizeC();
    sizeT = reader.getSizeT();

    % Validate indices
    if channelIndex > sizeC || channelIndex < 1
        error('Invalid channel index. File has %d channels.', sizeC);
    end
    if timeIndex > sizeT || timeIndex < 1
        error('Invalid time index. File has %d timepoints.', sizeT);
    end

    % Allocate output (assuming 16-bit data)
    zStack = zeros(sizeY, sizeX, sizeZ, 'uint16');

    % Read z-stack
    for z = 0:sizeZ-1
        index = reader.getIndex(z, channelIndex-1, timeIndex-1);  % 0-based
        zStack(:,:,z+1) = bfGetPlane(reader, index + 1);          % 1-based
    end

    reader.close();
end
