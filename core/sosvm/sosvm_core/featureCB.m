
function [phi] = featureCB(config, X, Y)
% featureCB  Compute the feature map. 
% [phi] = featureCB(config, x, y)
% OUTPUT: phi: feature map
% INPUT: config: configuration structure
%        x: cell-array with the training data
%        y: cell-array with a labeling.


    % Collect the unary features
    phi_u = get_features(X, Y);
    
    % Collect the pairwise features
    phi_p = get_pairwise_features_with_potts(X, Y);
    
    % Concatenate both the unary and the pairwise features and sum for all
    % the nodes
    phi = sparse(double(sum(cat(2, phi_u, phi_p))'));% / length(Y);   
    
end