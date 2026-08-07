
function out = ReadImage6D2(filename, useSeriesID, seriesID)
    switch nargin
        case 1
            useSeriesID = true;
            seriesID = 1;
        case 2
            seriesID = 1;
    end


    % --- TEMPORARILY SUPPRESS WARNINGS ---
    % Save the current warning state
    oldState = warning('off', 'all');
    % 
    % Get OME Meta-Information
    MetaData = GetOMEData(filename);

    % Initialize BioFormats Reader
    reader = bfGetReader(filename);

    % Determine which series (scene) to load
    if useSeriesID
        scenesToLoad = seriesID;
    else
        scenesToLoad = 1:MetaData.SeriesCount;
    end

    % Compute total number of frames for progress bar
    totalframes = numel(scenesToLoad) * MetaData.SizeC * MetaData.SizeZ * MetaData.SizeT;

    % Pre-allocate array:
    % Shape: (NumSeries, T, Z, C, Y, X)
    image6d = zeros(length(scenesToLoad), MetaData.SizeT, MetaData.SizeZ, ...
                    MetaData.SizeC, MetaData.SizeY, MetaData.SizeX);

    % Progress bar
    h = waitbar(0, 'Processing Data ...');
    framecounter = 0;

    % Loop over requested series only
    for i = 1:length(scenesToLoad)
        s = scenesToLoad(i);
        reader.setSeries(s - 1);  % 0-based for BioFormats

        for t = 1:MetaData.SizeT
            for z = 1:MetaData.SizeZ
                for c = 1:MetaData.SizeC
                    framecounter = framecounter + 1;

                    % Update waitbar
                    wstr = sprintf('Reading Images: %d of %d Frames', framecounter, totalframes);
                    waitbar(framecounter / totalframes, h, wstr);

                    % Compute plane index (0-based)
                    iplane = loci.formats.FormatTools.getIndex(reader, z - 1, c - 1, t - 1) + 1;

                    % Read and store the plane
                    image6d(i, t, z, c, :, :) = bfGetPlane(reader, iplane);
                end
            end
        end
    end

    % Clean up
    close(h);
    reader.close();

    % Store outputs
    out = cell(1, 2);
    out{1} = image6d;
    out{2} = MetaData;
end
