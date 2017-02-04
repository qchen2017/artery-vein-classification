
function [yhat] = constraintCB(config, model, X, Y)
% constraintCB Compute the most violated constraint
% [yhat] = constraintCB(config, model, x, y)
% OUTPUT: yhat: estimated labelling
% INPUT: config: configuration structure
%        model: learned model
%        X: a cell array containing the unary and pairwise features, the
%        pairwise ids and the original graph
%        Y: ground truth labelling

    % Separate the weights between unaries and pairwises
    [w_unary, w_pairwise] = get_weights(model.w, config);
        
    % Find the most violated constraint
    % 1) Get the unary features

    % Encode the unary features
    unary_features = get_features(X, Y); % Take only the unary features and the bias, and ignore the pairwise
    unary_features = reshape(unary_features, size(unary_features, 1), 2, config.size_w_u / 2);
    unary_features = permute(unary_features, [3 1 2]);
    
    % 2) Get the pairwise features
    pairwise_features = X{2};
    
    % 3) Include the penalizations into the unary features
    
    % As we want to get the most violated constraint, we will penalize each time
    % we do it well, so:
    % If the class has to be -1 and we classify it as 1, then we penalize (for veins)
    penalization_for_veins = zeros(size(Y));
    penalization_for_veins(Y==-1) = 1;
    % If the class has to be 1 but we classify it differently, then we penalize (for arteries)
    penalization_for_arteries = zeros(size(Y));
    penalization_for_arteries(Y==1) = 1;
       
    % 4) Compute the scores and multiply them by -1 to get the unary energy
    unary_potentials = zeros(size(unary_features, 2), 2);
    
    % scores for veins
    unary_potentials(:,1) = -(w_unary(1, :) * unary_features(:,:,1))' - penalization_for_veins;%sum(penalization_for_veins(:)); % veins
    % scores for arteries
    unary_potentials(:,2) = -(w_unary(2, :) * unary_features(:,:,2))' - penalization_for_arteries;%sum(penalization_for_arteries(:)); % arteries    

    % retrieve the number of elements
    N = size(pairwise_features, 1);
    % reshape the pairwise features so we can compute the pairwise
    % potentials easily
    pairwise_features = reshape(pairwise_features, N * N, config.size_w_p);
    % take the product between the pairwise features and their weights
    pairwise_potentials = pairwise_features * w_pairwise;
    % reshape the pairwise potentials to their original shape
    pairwise_potentials = reshape(pairwise_potentials, N, N);
    
    % 5) Get the most violated prediction   
    yhat = optimize_with_maxflow(sparse(unary_potentials), sparse(pairwise_potentials));
    
end
