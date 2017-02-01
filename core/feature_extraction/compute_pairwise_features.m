
function [Gout] = compute_pairwise_features(Gout, segm)

    % PAIRWISE FEATURES
    
    % ---------------------------------------------------------------------
    % Compute differences in vessel calibers
    % ---------------------------------------------------------------------

    % we use bwdist on the complement of the vessel segmentation 
    % to aproximate the vessel diameter.
    vessel_calibres_map = bwdist(imcomplement(segm));

    % initialize an array of average profiles, with the same number of
    % nodes that we have
    vessel_calibres = zeros(length(Gout.node), 1);
    
    % save each vessel calibre
    for i = 1 : length(Gout.node)
        % retrieve the average vessel calibre
        vessel_calibres(i) = mean(vessel_calibres_map(Gout.node(i).idx));
    end
    
    % lets compute the differences in vessel calibres
    difference_in_vessel_calibres = zeros(length(Gout.link), 1);
    % estimate the difference in vessel calibres
    for i = 1 : length(Gout.link)
        difference_in_vessel_calibres(i) = ...
            abs(vessel_calibres(Gout.link(i).n1) - vessel_calibres(Gout.link(i).n2)) / ...
                max(vessel_calibres(Gout.link(i).n1), vessel_calibres(Gout.link(i).n2));
    end
    
    % ---------------------------------------------------------------------
    % Compute angles
    % ---------------------------------------------------------------------
    
    % initialize the array of differences in angles
    difference_in_angles = zeros(length(Gout.link), 1);
    
    % for each link in our new graph of segments
    for i = 1 : length(Gout.link)
        
        % Get points of each segment
        segment_1 = get_last_pixels(Gout, Gout.node(Gout.link(i).n1), Gout.link(i));
        segment_2 = get_last_pixels(Gout, Gout.node(Gout.link(i).n2), Gout.link(i));
        % Normalize them by their norm
        segment_1 = segment_1;% / norm(segment_1);
        segment_2 = segment_2;% / norm(segment_2);
        
        % Identify the central point of the branching point / vessel
        % crossing
        point_coordinate = [ Gout.link(i).comx Gout.link(i).comy ];
        % Normalize it by its norm
        point_coordinate = point_coordinate;% / norm(point_coordinate);
        
        % Move the segments to the origin
        segment_1 = segment_1 - point_coordinate;
        segment_2 = segment_2 - point_coordinate;
        
        % and now take the angle between them
        difference_in_angles(i) = abs(segment_1 * segment_2') / (norm(segment_1) * norm(segment_2));
        
    end
    
    
    % ---------------------------------------------------------------------
    % Assign vessel calibres and angles
    % ---------------------------------------------------------------------
    
    % For each link in the graph
    for i = 1 : length(Gout.link)
        
        % Initialize current feature vector
        current_features = zeros(3, 1);
        
        % Determine if it corresponds to a crossing or a branching point
        if Gout.link(i).is_branching_point
            
            % Assign a 1 in the first coordinate because is a branching point
            current_features(1) = 1;
            % and 0 if it is not a crossing
            current_features(2) = 0;
            current_features(3) = 0;
            
        else
            
            % Assign a 0 because is a crossing
            current_features(1) = 0;
            % Assign the difference in vessel calibres
            current_features(2) = difference_in_vessel_calibres(i);
            % Assign the difference in angles
            current_features(3) = difference_in_angles(i);
            
        end
        
        % Assign current features
        Gout.link(i).features = current_features;
        
    end
    
    % Also include the dimensionality of the feature vector
    Gout.properties.pairwise_dim = length(Gout.link(1).features);

end


function last_pxs = get_last_pixels(Gout, node, link)

    % Check if the length of the segment is smaller than 10
    if length(node.idx) > 10
        max_dist = 10;
    else
        max_dist = length(node.idx) - 1;
    end
    
    % Get (x,y) coordinates of the segment
    segment_pxs = zeros(length(node.idx), 2);
    [segment_pxs(:,1), segment_pxs(:,2)] = ind2sub([Gout.w Gout.l], node.idx);

    % Estimate the distance with respect to current node
    distances = bsxfun(@minus, segment_pxs, [link.comx, link.comy]);
    % Take the Euclidean distance
    distances = sqrt(distances(:,1).^2 + distances(:,2).^2);
    % Identify the smallest value
    [~, idx] = min(distances);

    % Get the mean coordinate
    if idx==1
        last_pxs = segment_pxs(idx+max_dist, :);
    else
        last_pxs = segment_pxs(end-max_dist, :);
    end

end