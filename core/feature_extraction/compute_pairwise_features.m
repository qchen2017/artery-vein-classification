
function [Gout] = compute_pairwise_features(Gout, segm)

    % PAIRWISE FEATURES
    
    % --------------------------------------
    
    % 1. Mean vessel calibre

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
    
    % lets use the calibres now
    for i = 1 : length(Gout.link)
        % we compute the squared difference on the calibres as a pairwise
        % feature
        Gout.link(i).features = (vessel_calibres(Gout.link(i).n1) - vessel_calibres(Gout.link(i).n2))^2;
    end
    
    % --------------------------------------
    
    % 2. Bifurcation angles
    
    % for each link in our new graph of segments
    for i = 1 : length(Gout.link)
        
        % Get points of each segment
        segment_1 = [ Gout.node(Gout.link(i).n1).comx, Gout.node(Gout.link(i).n1).comy ];
        segment_2 = [ Gout.node(Gout.link(i).n2).comx, Gout.node(Gout.link(i).n2).comy ];
        % Normalize them by their norm
        segment_1 = segment_1 / norm(segment_1);
        segment_2 = segment_2 / norm(segment_2);
        
        % Identify the central point of the branching point / vessel
        % crossing
        point_coordinate = [ Gout.link(i).comx Gout.link(i).comy ];
        % Normalize it by its norm
        point_coordinate = point_coordinate / norm(point_coordinate);
        
        % Move the segments to the origin
        segment_1 = segment_1 - point_coordinate;
        segment_2 = segment_2 - point_coordinate;
        
        % and now take the angle between them
        Gout.link(i).features = [ Gout.link(i).features, 180 - abs(acosd(segment_1 * segment_2')) ];
        
    end
     
    % Also include the dimensionality of the feature vector
    Gout.properties.pairwise_dim = length(Gout.link(1).features);

end