%%%% ===== for PRIMARY/ SECONDARY/ TERTIARY condensates

function out = condensate_search_v6(stack_C, BW3, enh_pos, spacing, t_cond, mult, R_search_c)

out = cell(3, 1);

xpixel = spacing(1);
ypixel = spacing(2);
zpixel = spacing(3);

enh_x = enh_pos(1);
enh_y = enh_pos(2);
enh_z = enh_pos(3);

C_cent_e_x = enh_x/xpixel;
C_cent_e_y = enh_y/ypixel;

max_iterations = 200;
iteration = 1;
s_c = 0;

% Minimum distance between condensates (in nm) 
min_separation_nm = 500; % minimum separation between condensates

% Store all found condensates here
all_condensates = [];

% Iteratively increase search radius until at least 1 condensate is found
while s_c < 1 && iteration <= max_iterations
    current_R_search = iteration * R_search_c;
    
    masked_img_cond = imagemask(BW3, C_cent_e_x, C_cent_e_y, current_R_search);
    pk_cond = pkfnd(masked_img_cond, mult * t_cond, 15);
    s_c = size(pk_cond, 1);
    
    if s_c < 1
        fprintf('No condensates found with radius %.1f (iteration %d). Expanding search...\n', current_R_search, iteration);
        iteration = iteration + 1;
    end
end

if s_c >= 1
    fprintf('Found %d condensate(s) with search radius %.1f (%.1fx original)\n', s_c, current_R_search, current_R_search/R_search_c);
    
    bright_spots_y_c = pk_cond(:, 2);
    bright_spots_x_c = pk_cond(:, 1);
    num_spots_c = size(pk_cond, 1);
    I_spots_c = zeros(num_spots_c, 1);
    fit_spot_c = zeros(num_spots_c, 1);
    fit_spot_ar = zeros(num_spots_c, 1);
    IC_fit = zeros(num_spots_c, 1);
    spotC_xyz = zeros(num_spots_c, 3);
    R_c_3D_fit = 6;
    
    for i = 1:num_spots_c
        index1 = bright_spots_x_c(i, 1);
        index2 = bright_spots_y_c(i, 1);
        I_spots_c(i, 1) = BW3(index2, index1);
        
        B_cond = fitGaussian2(masked_img_cond, I_spots_c(i, 1), index1, index2, 1.5, 3);
        fit_spot_c(i, 1) = (B_cond(4) * B_cond(5))^0.5;
        fit_spot_ar(i, 1) = max(B_cond(4), B_cond(5)) / min(B_cond(4), B_cond(5));
        IC_fit(i, 1) = B_cond(1);
        
        guess_3D = [B_cond(2), B_cond(3), enh_z/zpixel + 1];
        [x_ref_c, y_ref_c, z_ref_c, ~] = fit_Gaussian3D(stack_C, guess_3D, R_c_3D_fit, spacing);
        spotC_xyz(i, 1) = x_ref_c * xpixel;
        spotC_xyz(i, 2) = y_ref_c * ypixel;
        spotC_xyz(i, 3) = (z_ref_c - 1) * zpixel;
    end
    
    EC_dist = zeros(num_spots_c, 1);
    for i = 1:num_spots_c
        EC_dist(i, 1) = (((enh_x - spotC_xyz(i, 1))^2 + (enh_y - spotC_xyz(i, 2))^2 + (enh_z - spotC_xyz(i, 3))^2)^0.5) * 1000;
    end
    
    [sorted_distances, sorted_indices] = sort(EC_dist, 'ascend');
    
    % NEW: Add condensates one by one, checking for minimum separation
    for i = 1:num_spots_c
        current_index = sorted_indices(i);
        current_xyz = spotC_xyz(current_index, :);
        current_distance = EC_dist(current_index); % Get the actual distance for THIS condensate
        
        % Check if this condensate is too close to any already-added condensate
        is_too_close = false;
        for j = 1:size(all_condensates, 1)
            existing_xyz = all_condensates(j, 1:3);
            distance_between = sqrt(sum((current_xyz - existing_xyz).^2)) * 1000; % in nm
            if distance_between < min_separation_nm
                is_too_close = true;
                fprintf('Rejecting condensate at distance %.1f nm (too close to existing condensate: %.1f nm apart)\n', ...
                        current_distance, distance_between);
                break;
            end
        end
        
        % Only add if sufficiently separated from all existing condensates
        if ~is_too_close
            new_condensate = [spotC_xyz(current_index, 1), spotC_xyz(current_index, 2), spotC_xyz(current_index, 3), ...
                              fit_spot_c(current_index, 1), IC_fit(current_index, 1), fit_spot_ar(current_index, 1), current_distance];
            all_condensates = [all_condensates; new_condensate];
            fprintf('Added condensate #%d at distance %.1f nm\n', size(all_condensates, 1), current_distance);
            
            if size(all_condensates, 1) >= 3
                break; % Found enough condensates
            end
        end
    end
    
    num_found = size(all_condensates, 1);
    
    % If we found fewer than 3, search in annular rings
    if num_found < 3
        fprintf('Only found %d condensate(s) initially. Searching annular rings...\n', num_found);
        
        inner_radius = current_R_search;
        annular_iteration = 1;
        max_annular_iterations = 20;
        
        % Search annular rings for remaining condensates
        while num_found < 3 && annular_iteration <= max_annular_iterations
            outer_radius = inner_radius + R_search_c;
            
            fprintf('Searching annular ring: inner=%.1f, outer=%.1f\n', inner_radius, outer_radius);
            
            % Create annular mask
            masked_img_outer = imagemask(BW3, C_cent_e_x, C_cent_e_y, outer_radius);
            masked_img_inner = imagemask(BW3, C_cent_e_x, C_cent_e_y, inner_radius);
            masked_img_annular = masked_img_outer - masked_img_inner;
            
            pk_cond_annular = pkfnd(masked_img_annular, mult * t_cond, 5);
            s_c_annular = size(pk_cond_annular, 1);
            
            if s_c_annular >= 1
                % Process new condensates found in annular ring
                bright_spots_y_ann = pk_cond_annular(:, 2);
                bright_spots_x_ann = pk_cond_annular(:, 1);
                num_spots_ann = size(pk_cond_annular, 1);
                
                spotC_xyz_ann = zeros(num_spots_ann, 3);
                fit_spot_c_ann = zeros(num_spots_ann, 1);
                fit_spot_ar_ann = zeros(num_spots_ann, 1);
                IC_fit_ann = zeros(num_spots_ann, 1);
                
                for i = 1:num_spots_ann
                    index1 = bright_spots_x_ann(i, 1);
                    index2 = bright_spots_y_ann(i, 1);
                    I_spot = masked_img_annular(index2, index1);
                    
                    B_cond = fitGaussian2(masked_img_annular, I_spot, index1, index2, 1.5, 3);
                    fit_spot_c_ann(i, 1) = (B_cond(4) * B_cond(5))^0.5;
                    fit_spot_ar_ann(i, 1) = max(B_cond(4), B_cond(5)) / min(B_cond(4), B_cond(5));
                    IC_fit_ann(i, 1) = B_cond(1);
                    
                    guess_3D = [B_cond(2), B_cond(3), enh_z/zpixel + 1];
                    [x_ref_c, y_ref_c, z_ref_c, ~] = fit_Gaussian3D(stack_C, guess_3D, R_c_3D_fit, spacing);
                    spotC_xyz_ann(i, 1) = x_ref_c * xpixel;
                    spotC_xyz_ann(i, 2) = y_ref_c * ypixel;
                    spotC_xyz_ann(i, 3) = (z_ref_c - 1) * zpixel;
                end
                
                EC_dist_ann = zeros(num_spots_ann, 1);
                for i = 1:num_spots_ann
                    EC_dist_ann(i, 1) = (((enh_x - spotC_xyz_ann(i, 1))^2 + (enh_y - spotC_xyz_ann(i, 2))^2 + (enh_z - spotC_xyz_ann(i, 3))^2)^0.5) * 1000;
                end
                
                [sorted_distances_ann, sorted_indices_ann] = sort(EC_dist_ann, 'ascend');
                
                % NEW: Add condensates from annular ring with separation check
                for i = 1:num_spots_ann
                    current_index = sorted_indices_ann(i);
                    current_xyz = spotC_xyz_ann(current_index, :);
                    
                    % Check if this condensate is too close to any already-added condensate
                    is_too_close = false;
                    for j = 1:size(all_condensates, 1)
                        existing_xyz = all_condensates(j, 1:3);
                        distance_between = sqrt(sum((current_xyz - existing_xyz).^2)) * 1000; % in nm
                        if distance_between < min_separation_nm
                            is_too_close = true;
                            fprintf('Rejecting annular condensate (too close: %.1f nm apart)\n', distance_between);
                            break;
                        end
                    end
                    
                    % Only add if sufficiently separated
                    if ~is_too_close
                        new_condensate = [spotC_xyz_ann(current_index, 1), spotC_xyz_ann(current_index, 2), spotC_xyz_ann(current_index, 3), ...
                                          fit_spot_c_ann(current_index, 1), IC_fit_ann(current_index, 1), fit_spot_ar_ann(current_index, 1), EC_dist_ann(current_index)];
                        all_condensates = [all_condensates; new_condensate];
                        num_found = num_found + 1;
                        fprintf('Found condensate #%d at distance %.1f nm\n', num_found, EC_dist_ann(current_index));
                        
                        if num_found >= 3
                            break; % Found enough
                        end
                    end
                end
            end
            
            if num_found >= 3
                break;
            end
            
            inner_radius = outer_radius;
            annular_iteration = annular_iteration + 1;
        end
    end
    
    % Sort all condensates by distance from enhancer (column 7)
    if ~isempty(all_condensates)
        [~, sort_idx] = sort(all_condensates(:, 7), 'ascend');
        all_condensates = all_condensates(sort_idx, :);
        fprintf('\nFinal sorted condensates by distance:\n');
        for j = 1:size(all_condensates, 1)
            fprintf('  Condensate #%d: %.1f nm from enhancer\n', j, all_condensates(j, 7));
        end
    end
    
    % Assign the found condensates to output
    for j = 1:min(size(all_condensates, 1), 3)
        out{j} = all_condensates(j, :);
    end
    
    % Fill remaining cells with NaN if fewer than 3 total found
    for j = (size(all_condensates, 1) + 1):3
        out{j} = NaN;
    end
    
else % s_c == 0 (no condensates found even after expanding search)
    out{1} = NaN;
    out{2} = NaN;
    out{3} = NaN;
end

end