
function [phi] = get_features(patterns, labels)
% GET_FEATURES Computes the Kroneker product between the unary features and
% the label vector

    % Get the feature vectors
    unary_features = patterns{1};
    
    % Compute the unary features
    phi_u = zeros(size(unary_features, 1), size(unary_features, 2) * 2);
    
    % Take the Kronecker product of the features with the corresponding
    % binary vector, according to the given labeling y
    
    % veins (arteries are 1, unknown are 0s)
    phi_u(labels == -1, :) = kron(unary_features(labels == -1, :), [1 0]);
    % arteries
    phi_u(labels == 1, :) = kron(unary_features(labels == 1, :), [0 1]);
    
    % Return the unary features
    phi = phi_u;
    
end
