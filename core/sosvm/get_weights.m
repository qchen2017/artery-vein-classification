
function [w_unaries, w_pairwises] = get_weights(W, config)
% get_weights Separates the weights for unary, pairwise and bias
% [W_unaries, W_pairwises, bias] = getWeights(W, config)
% OUTPUT: W_unaries: weights for the unary potentials
%         W_pairwises: weights for the pairwise potentials
            
    % Separate the pairwise weights
    w_pairwises = W(end - config.size_w_p + 1 : end);
    w_no_pairwises = W(1 : end - config.size_w_p);

    % separate the weights of each class in such a way that i indicates the
    % number of class and j the number of feature
    w_unaries = reshape(w_no_pairwises, 2, config.size_w_u/2);
    
end