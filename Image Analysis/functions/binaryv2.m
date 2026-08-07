function [starts_1s, ends_1s, lengths_1s, starts_0s, ends_0s, lengths_0s] = binaryv2(v)
    % Pads input vector to find transitions
    v_padded = [NaN, v, NaN];
    edges = find(diff(v_padded) ~= 0);

    % Compute all run lengths and indices
    run_lengths = diff(edges);
    run_values = v(edges(1:end-1));
    starts_all = edges(1:end-1);
    ends_all = edges(2:end) - 1;

    % Logical indices for 1s and 0s
    idx_1s = find(run_values == 1);
    idx_0s = find(run_values == 0);

    % Extract runs of 1s
    starts_1s = starts_all(idx_1s);
    ends_1s = ends_all(idx_1s);
    lengths_1s = run_lengths(idx_1s);

    % Extract runs of 0s
    starts_0s = starts_all(idx_0s);
    ends_0s = ends_all(idx_0s);
    lengths_0s = run_lengths(idx_0s);
end
