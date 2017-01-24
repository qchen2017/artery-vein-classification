
function [Gout] = compute_unary_features(Gout, RGB)

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
    
end


function current_band = normalize_intensities(current_band)
    % turn the image into doubles
    current_band = im2double(current_band);
    % and now standardize
    current_band = (current_band - mean(current_band(:))) / std(current_band(:));
end