function [ G ] = vesselsCalibre(G, RGB)

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
    
    displacement = 3; % to compute vessel normals
    num_profiles = 10; % to sample some normals
    profile_size = 5;
    
%     figure, imshow(RGB);
%     hold on
    
    
    % Compute intensity profile
    for i = 1 : length(G.node)

        sample_factor = floor(length(G.node(i).idx) / num_profiles);
        if (sample_factor < 1)
            sample_factor = 1;
        end
        
        % identify points in the segment
        pts = zeros(length(G.node(i).idx),2);
        [pts(:,1), pts(:,2)] = ind2sub([G.w G.l], G.node(i).idx);
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
        v_sampled = next_pts(mod(1:length(pts),sample_factor)==1,:) - pts(mod(1:length(pts),sample_factor)==1,:);
        if (size(v_sampled,1) < 3)
            v_sampled = next_pts - pts;
        else
            pts = pts(mod(1:length(pts),sample_factor)==1,:);
            next_pts = next_pts(mod(1:length(next_pts),sample_factor)==1,:);
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

        G.node(i).features = [];
        
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
            G.node(i).features = cat(1, G.node(i).features, mean(current_profiles));

%             figure, plot(current_profiles', '--or');
%             ylim([-1 1]);
%             hold on
%             plot(mean(current_profiles)', 'LineWidth', 3);
%             grid on
            
        end
    end

end    
    
    
function output_I = normalize_intensities(I)
    % initialize the output parameter
    output_I = zeros(size(I));
    for i = 1 : size(I, 3)
        % turn the image into doubles
        current_band = im2double(I(:,:,i));
        % and now standardize
        output_I(:,:,i) = (current_band - mean(current_band(:))) / std(current_band(:));
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

