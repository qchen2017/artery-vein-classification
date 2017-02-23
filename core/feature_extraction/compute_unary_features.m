

function [Gout] = compute_unary_features(Gout, RGB, segm)

    %[Gout] = compute_unary_features_intensities(Gout, RGB);
    %[Gout] = compute_unary_features_profiles(Gout, RGB);
    %[Gout] = compute_unary_features_contrast(Gout, RGB, segm);
    [Gout] = compute_unary_features_profiles_statistics(Gout, RGB, segm);

end






function [Gout] = compute_unary_features_profiles_statistics(Gout, RGB, segm)

    % Prepare auxiliar color planes
    color_planes = zeros(size(RGB,1), size(RGB,2), 5);
    
    % RGB color image
    color_planes(:,:,1) = normalize_intensities(RGB(:,:,1)); % red
    color_planes(:,:,2) = normalize_intensities(RGB(:,:,2)); % green
    color_planes(:,:,3) = normalize_intensities(RGB(:,:,3)); % blue
    
    % HSV color space
    HSV = rgb2hsv(RGB);
    color_planes(:,:,4) = normalize_intensities(HSV(:,:,1)); % hue
    color_planes(:,:,5) = normalize_intensities(HSV(:,:,2)); % saturation
    color_planes(:,:,6) = normalize_intensities(HSV(:,:,3)); % value
    
    % Equalize intensities
    color_planes(:,:,7:9) = normalize_intensities(contrastEqualization(RGB));
    
    num_profiles = 10; % to sample some normals
    profile_size = 5;
    
%     figure, imshow(RGB);
%     hold on
    
    
    % Compute intensity profile
    for i = 1 : length(Gout.node)

        sample_factor = floor(length(Gout.node(i).idx) / num_profiles);
        
        % identify points in the segment
        pts = zeros(length(Gout.node(i).idx),2);
        [pts(:,1), pts(:,2)] = ind2sub([Gout.w Gout.l], Gout.node(i).idx);
        
        if (isempty(pts(3:end-3,:)))
            displacement = 1;
        else
            displacement = 3;
        end
        
        % en next_pts tengo las mismas coordenadas, pero desplazadas de
        % acuerdo al valor de displacement
        next_pts = circshift(pts',[0 -displacement])'; 
        % by default, we capture the profile size
        vector_norm = profile_size;

        % we remove the points in each corner to avoid wrong normals
        pts = pts(displacement:end-displacement,:);
        next_pts = next_pts(displacement:end-displacement,:);

        % in case we want less than n profiles we can use sample_factor
        % to recover n/sample_factor elements
        %v_sampled = next_pts(mod(1:length(pts),sample_factor)==1,:) - pts(mod(1:length(pts),sample_factor)==1,:);
        if (sample_factor <= 1)
            v_sampled = next_pts - pts;
        else
            pts = pts(1:sample_factor:size(pts,1),:);
            next_pts = next_pts(1:sample_factor:size(next_pts,1),:);
            v_sampled = next_pts - pts;
        end

        % compute vectors orthogonal to the skeleton
        norm_ = ones(size(v_sampled,1),1) * vector_norm;
        u_sampled = zeros(size(v_sampled));
        u_sampled(:,1) = norm_ ./ sqrt(1 + (v_sampled(:,1) ./ v_sampled(:,2)).^2);
        u_sampled(:,2) = -1 * (v_sampled(:,1) ./ v_sampled(:,2)) .* u_sampled(:,1);

        % en caso de que haya alg?n delta_y=0, se corrige:
        u_sampled(v_sampled(:,1)==0, 2) = 0;
        u_sampled(v_sampled(:,2)==0, 2) = vector_norm;

        % calculo los puntos ortogonales al actual
        orthogonal_a = round(pts + u_sampled);
        orthogonal_b = round(pts - u_sampled);

%         quiver(pts(:,2), pts(:,1), u_sampled(:,2), u_sampled(:,1), 'color', [1 1 1]);
%         hold on
%         quiver(pts(:,2), pts(:,1), -u_sampled(:,2), -u_sampled(:,1), 'color', [1 1 1]);
%         hold on

        Gout.node(i).features = [];
        
        % initialize an array with all the border trails
        profiles = cell(size(orthogonal_a, 1), 1);

        % for each profile
        for j = 1 : size(orthogonal_a,1)
            
            % recover each profile and estimate which points are inside the
            % vessel
            [idx_x, idx_y, vessel_profile] = improfile(segm, [orthogonal_a(j,2) orthogonal_b(j,2)], [orthogonal_a(j,1) orthogonal_b(j,1)], double(vector_norm*2+1));
            idx_x = round(idx_x(logical(vessel_profile)));
            idx_y = round(idx_y(logical(vessel_profile)));
            
            % assign profile
            profiles{j} = cat(2, idx_y, idx_x);

        end
        
        % initialize feature arrays
        max_intensities = zeros(size(orthogonal_a, 1), size(color_planes, 3));
        min_intensities = zeros(size(max_intensities));
        avg_intensities = zeros(size(max_intensities));
        std_intensities = zeros(size(max_intensities));
        width_intensities = zeros(size(orthogonal_a, 1), 1);
        
        
        % for each profile
        for j = 1 : size(orthogonal_a, 1)
            
            coordinates = profiles{j};
            
            % for each color band
            for color_idx = 1 : size(color_planes, 3)
                
                color_plane = color_planes(:,:,color_idx);
                
                max_intensities(j, color_idx) = max(color_plane(sub2ind(size(color_plane), coordinates(:,1), coordinates(:,2))));
                min_intensities(j, color_idx) = min(color_plane(sub2ind(size(color_plane), coordinates(:,1), coordinates(:,2))));
                avg_intensities(j, color_idx) = mean(color_plane(sub2ind(size(color_plane), coordinates(:,1), coordinates(:,2))));
                std_intensities(j, color_idx) = std(color_plane(sub2ind(size(color_plane), coordinates(:,1), coordinates(:,2))));
            end
            % get vessel width
            width_intensities(j) = size(profiles{j}, 1);
        
        end
        
        if size(orthogonal_a, 1)==1
            Gout.node(i).features = cat(2, max_intensities, min_intensities, avg_intensities, std_intensities, width_intensities)';
        else
            Gout.node(i).features = cat(2, mean(max_intensities), mean(min_intensities), mean(avg_intensities), mean(std_intensities), mean(width_intensities))';
        end
        if any(isnan(Gout.node(i).features))
            disp('A');
        end
        
    end
    
    % Also include the dimensionality of the feature vector
    Gout.properties.unary_dim = length(Gout.node(1).features);
    
end



function [Gout] = compute_unary_features_contrast(Gout, RGB, segm)

    % Prepare auxiliar color planes
    color_planes = zeros(size(RGB,1), size(RGB,2), 9);
    
    % RGB color image
    color_planes(:,:,1) = normalize_intensities(RGB(:,:,1)); % red
    color_planes(:,:,2) = normalize_intensities(RGB(:,:,2)); % green
    color_planes(:,:,3) = normalize_intensities(RGB(:,:,3)); % blue
    
    % HSV color space
    HSV = rgb2hsv(RGB);
    color_planes(:,:,4) = normalize_intensities(HSV(:,:,1)); % hue
    color_planes(:,:,5) = normalize_intensities(HSV(:,:,2)); % saturation
    color_planes(:,:,6) = normalize_intensities(HSV(:,:,3)); % value
    
    % Equalize intensities
    color_planes(:,:,7:9) = normalize_intensities(contrastEqualization(RGB));
    
    num_profiles = 10; % to sample some normals
    profile_size = 5;
    
%     figure, imshow(color_planes(:,:,7:9));
%     hold on
    
    % Compute intensity profile
    for i = 1 : length(Gout.node)

        sample_factor = floor(length(Gout.node(i).idx) / num_profiles);
        
        % identify points in the segment
        pts = zeros(length(Gout.node(i).idx),2);
        [pts(:,1), pts(:,2)] = ind2sub([Gout.w Gout.l], Gout.node(i).idx);
        
        if (isempty(pts(3:end-3,:)))
            displacement = 1;
        else
            displacement = 3;
        end
        
        % en next_pts tengo las mismas coordenadas, pero desplazadas de
        % acuerdo al valor de displacement
        next_pts = circshift(pts',[0 -displacement])'; 
        % by default, we capture the profile size
        vector_norm = profile_size;

        % we remove the points in each corner to avoid wrong normals
        pts = pts(displacement:end-displacement,:);
        next_pts = next_pts(displacement:end-displacement,:);

        % in case we want less than n profiles we can use sample_factor
        % to recover n/sample_factor elements
        %v_sampled = next_pts(mod(1:length(pts),sample_factor)==1,:) - pts(mod(1:length(pts),sample_factor)==1,:);
        if (sample_factor <= 1)
            v_sampled = next_pts - pts;
        else
            pts = pts(1:sample_factor:size(pts,1),:);
            next_pts = next_pts(1:sample_factor:size(next_pts,1),:);
            v_sampled = next_pts - pts;
        end

        % compute vectors orthogonal to the skeleton
        norm_ = ones(size(v_sampled,1),1) * vector_norm;
        u_sampled = zeros(size(v_sampled));
        u_sampled(:,1) = norm_ ./ sqrt(1 + (v_sampled(:,1) ./ v_sampled(:,2)).^2);
        u_sampled(:,2) = -1 * (v_sampled(:,1) ./ v_sampled(:,2)) .* u_sampled(:,1);

        % en caso de que haya alg?n delta_y=0, se corrige:
        u_sampled(v_sampled(:,1)==0, 2) = 0;
        u_sampled(v_sampled(:,2)==0, 2) = vector_norm;

        % calculo los puntos ortogonales al actual
        orthogonal_a = round(pts + u_sampled);
        orthogonal_b = round(pts - u_sampled);

%         quiver(pts(:,2), pts(:,1), u_sampled(:,2), u_sampled(:,1), 'color', [1 1 1]);
%         hold on
%         quiver(pts(:,2), pts(:,1), -u_sampled(:,2), -u_sampled(:,1), 'color', [1 1 1]);
%         hold on

        Gout.node(i).features = [];
        
        % initialize an array with all the border trails
        border_trails = uint32(zeros(2, 2, size(orthogonal_a,1)));
        type_of_vessel = zeros(1, size(orthogonal_a, 1));

        % for each profile
        for j = 1 : size(orthogonal_a,1)
            
            % recover each profile and estimate which points are inside the
            % vessel
            [idx_x, idx_y, vessel_profile] = improfile(segm, [orthogonal_a(j,2) orthogonal_b(j,2)], [orthogonal_a(j,1) orthogonal_b(j,1)], double(vector_norm*2+1));
            vessel_profile_idx = find(vessel_profile);
            
            % independently of the vessel width, we will always have one
            % point
            border_trails(1,:,j) = [idx_y(vessel_profile_idx(1)); idx_x(vessel_profile_idx(1))];
            
            % vessel type 1: just 1 pixel width
            if length(vessel_profile_idx) == 1
                type_of_vessel(j) = 1;
                
            % vessel type 2: 2 pixel width
            elseif length(vessel_profile_idx) == 2
                type_of_vessel(j) = 2;
                % if is the same coordinate than that the one in the
                % skeleton, replace it
                if pts(j,:) == [idx_y(vessel_profile_idx(1)); idx_x(vessel_profile_idx(1))]
                    border_trails(1,:,j) = [idx_y(vessel_profile_idx(end)); idx_x(vessel_profile_idx(end))];
                end
            % vessel type 3: vessel width > 2
            else
                type_of_vessel(j) = 3;
                % in this case, we have also another point
                border_trails(2,:,j) = [idx_y(vessel_profile_idx(end)); idx_x(vessel_profile_idx(end))];
            end

        end
        
        % Initialize the feature vector
        Gout.node(i).features = zeros(size(color_planes, 3), 1);
        
        % For each color plane
        for color_idx = 1 : size(color_planes, 3)
        
            % initialize a matrix of difference in contrast
            differences_in_contrast = zeros(size(border_trails, 3), 1);
            
            % for each border trail coordinate
            for border_trail_idx = 1 : size(border_trails, 3)
                
                switch type_of_vessel(border_trail_idx)
                    
                    case 1
                        differences_in_contrast(border_trail_idx) = 0;
                    case 2
                        differences_in_contrast(border_trail_idx) = ...
                            color_planes(pts(border_trail_idx, 1), pts(border_trail_idx, 2), color_idx) - ...
                                color_planes(border_trails(1,1,border_trail_idx), border_trails(1,2,border_trail_idx), color_idx);
                    case 3
                        my_differences = zeros(2,1);
                        my_differences(1) = color_planes(pts(border_trail_idx, 1), pts(border_trail_idx, 2), color_idx) - ...
                                color_planes(border_trails(1,1,border_trail_idx), border_trails(1,2,border_trail_idx), color_idx);
                        my_differences(2) = color_planes(pts(border_trail_idx, 1), pts(border_trail_idx, 2), color_idx) - ...
                                color_planes(border_trails(2,1,border_trail_idx), border_trails(2,2,border_trail_idx), color_idx);
                        differences_in_contrast(border_trail_idx) = mean(my_differences);
                end
                
            end
            
            Gout.node(i).features(color_idx) = mean(differences_in_contrast);
            
            if any(isnan(Gout.node(i).features))
                disp('A');
            end
%             figure, plot(current_profiles', '--or');
%             ylim([-1 1]);
%             hold on
%             plot(mean(current_profiles)', 'LineWidth', 3);
%             grid on
            
        end
        
        % Add a bias term
        %Gout.node(i).features = cat(1, Gout.node(i).features, 1);
        
    end
    
    % Also include the dimensionality of the feature vector
    Gout.properties.unary_dim = length(Gout.node(1).features);
    
end


function [Gout] = compute_unary_features_intensities(Gout, RGB)%, segm)

    % UNARY FEATURES
    
    % Prepare auxiliar color planes
    
    % RGB color image
    red = normalize_intensities(RGB(:,:,1));
    green = normalize_intensities(RGB(:,:,2));
    blue = normalize_intensities(RGB(:,:,3));
    
    % HSV color space
    HSV = rgb2hsv(RGB);
    H = normalize_intensities(HSV(:,:,1));
    S = normalize_intensities(HSV(:,:,2));
    V = normalize_intensities(HSV(:,:,3));
    
    
    % --------------------------------------
    
    
    % --------------------------------------
    
    % for each node in the graph
    for i = 1 : length(Gout.node)
        
        % I have 12 unary features so far
        current_node_unary_features = zeros(12, 1);
        
        % Retrieve node pixels
        pxs = Gout.node(i).idx;
        
        % get mean intensities on each color plane
        current_node_unary_features(1) = mean(red(pxs));
        current_node_unary_features(2) = mean(green(pxs));
        current_node_unary_features(3) = mean(blue(pxs));
        current_node_unary_features(4) = mean(H(pxs));
        current_node_unary_features(5) = mean(S(pxs));
        current_node_unary_features(6) = mean(V(pxs));
        
        % get std intensities on each color plane
        current_node_unary_features(7)  = std(red(pxs));
        current_node_unary_features(8)  = std(green(pxs));
        current_node_unary_features(9)  = std(blue(pxs));
        current_node_unary_features(10) = std(H(pxs));
        current_node_unary_features(11) = std(S(pxs));
        current_node_unary_features(12) = std(V(pxs));
        
        % assign current features
        Gout.node(i).features = current_node_unary_features;
        
    end
    
    % Also include the dimensionality of the feature vector
    Gout.properties.unary_dim = length(current_node_unary_features);
    
end


function [Gout] = compute_unary_features_profiles(Gout, RGB)

    % Prepare auxiliar color planes
    color_planes = zeros(size(RGB,1), size(RGB,2), 9);
    
    % RGB color image
    color_planes(:,:,1) = normalize_intensities(RGB(:,:,1)); % red
    color_planes(:,:,2) = normalize_intensities(RGB(:,:,2)); % green
    color_planes(:,:,3) = normalize_intensities(RGB(:,:,3)); % blue
    
    % HSV color space
    HSV = rgb2hsv(RGB);
    color_planes(:,:,4) = normalize_intensities(HSV(:,:,1)); % hue
    color_planes(:,:,5) = normalize_intensities(HSV(:,:,2)); % saturation
    color_planes(:,:,6) = normalize_intensities(HSV(:,:,3)); % value
    
    % Equalize intensities
    color_planes(:,:,7:9) = normalize_intensities(contrastEqualization(RGB));
    
    num_profiles = 10; % to sample some normals
    profile_size = 5;
    
%     figure, imshow(RGB);
%     hold on
    
    
    % Compute intensity profile
    for i = 1 : length(Gout.node)

        sample_factor = floor(length(Gout.node(i).idx) / num_profiles);
        
        % identify points in the segment
        pts = zeros(length(Gout.node(i).idx),2);
        [pts(:,1), pts(:,2)] = ind2sub([Gout.w Gout.l], Gout.node(i).idx);
        
        if (isempty(pts(3:end-3,:)))
            displacement = 1;
        else
            displacement = 3;
        end
        
        % en next_pts tengo las mismas coordenadas, pero desplazadas de
        % acuerdo al valor de displacement
        next_pts = circshift(pts',[0 -displacement])'; 
        % by default, we capture the profile size
        vector_norm = profile_size;

        % we remove the points in each corner to avoid wrong normals
        pts = pts(displacement:end-displacement,:);
        next_pts = next_pts(displacement:end-displacement,:);

        % in case we want less than n profiles we can use sample_factor
        % to recover n/sample_factor elements
        %v_sampled = next_pts(mod(1:length(pts),sample_factor)==1,:) - pts(mod(1:length(pts),sample_factor)==1,:);
        if (sample_factor <= 1)
            v_sampled = next_pts - pts;
        else
            pts = pts(1:sample_factor:size(pts,1),:);
            next_pts = next_pts(1:sample_factor:size(next_pts,1),:);
            v_sampled = next_pts - pts;
        end

        % compute vectors orthogonal to the skeleton
        norm_ = ones(size(v_sampled,1),1) * vector_norm;
        u_sampled = zeros(size(v_sampled));
        u_sampled(:,1) = norm_ ./ sqrt(1 + (v_sampled(:,1) ./ v_sampled(:,2)).^2);
        u_sampled(:,2) = -1 * (v_sampled(:,1) ./ v_sampled(:,2)) .* u_sampled(:,1);

        % en caso de que haya alg?n delta_y=0, se corrige:
        u_sampled(v_sampled(:,1)==0, 2) = 0;
        u_sampled(v_sampled(:,2)==0, 2) = vector_norm;

        % calculo los puntos ortogonales al actual
        orthogonal_a = round(pts + u_sampled);
        orthogonal_b = round(pts - u_sampled);

%         quiver(pts(:,2), pts(:,1), u_sampled(:,2), u_sampled(:,1), 'color', [1 1 1]);
%         hold on
%         quiver(pts(:,2), pts(:,1), -u_sampled(:,2), -u_sampled(:,1), 'color', [1 1 1]);
%         hold on

        Gout.node(i).features = [];
        
        % For each color plane
        for color_idx = 1 : size(color_planes, 3)
        
            % initialize an array for each intensity profile
            current_profiles = zeros(size(u_sampled, 1), profile_size * 2 + 1);
            
            % for each profile
            for j = 1 : size(orthogonal_a,1)
                % recover each profile
                current_profiles(j,:) = improfile(color_planes(:,:,color_idx), [orthogonal_a(j,2) orthogonal_b(j,2)], [orthogonal_a(j,1) orthogonal_b(j,1)], double(vector_norm*2+1));
            end
            % get mean profile and assign to features
            if size(current_profiles, 1)==1
                Gout.node(i).features = cat(1, Gout.node(i).features, current_profiles');
            else
                Gout.node(i).features = cat(1, Gout.node(i).features, mean(current_profiles)');
            end
            if any(isnan(Gout.node(i).features))
                disp('A');
            end
%             figure, plot(current_profiles', '--or');
%             ylim([-1 1]);
%             hold on
%             plot(mean(current_profiles)', 'LineWidth', 3);
%             grid on
            
        end
        
        % Add a bias term
        %Gout.node(i).features = cat(1, Gout.node(i).features, 1);
        
    end
    
    % Also include the dimensionality of the feature vector
    Gout.properties.unary_dim = length(Gout.node(1).features);
    
end




function output_I = normalize_intensities(I)
    % initialize the output parameter
    output_I = zeros(size(I));
    for i = 1 : size(I, 3)
        % turn the image into doubles
        output_I(:,:,i) = im2double(I(:,:,i));
        % and now standardize
        %output_I(:,:,i) = (current_band - mean(current_band(:))) / std(current_band(:));
    end
end



function [I_out] = contrastEqualization(I)

    % get mask
    mask = get_fov_mask(I, 0.5);
    
    % initialize the new image
    I_out = uint8(zeros(size(I)));
    w = floor(3*(size(I,1))/30);
    
    % for each color plane
    for i = 1 : size(I,3)
        
        % fakepad current color band
        [I_extended, mask_extended] = fakepad(I(:,:,i), mask, 5, w);
        % apply gaussian filter
        G = imfilter(I_extended, fspecial('gaussian', [w, w], (size(I,1))/30));
        % rebuild image
        I_extended = 4 * double(I_extended) - 4 * double(G) + 128;
        % rebuild current color band
        I_current = zeros(size(I(:,:,i)));
        I_current(mask) = I_extended(mask_extended>0);
        % assign current band to the new image
        I_out(:,:,i) = uint8(I_current);
        
    end
end
