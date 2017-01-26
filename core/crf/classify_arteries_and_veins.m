
function [Gout, labelling, ground_truth] = classify_arteries_and_veins(Gout, model)

    % if the graph has labels, lets copy them to the ground truth variable
    if isfield(Gout, 'labels')
        ground_truth = Gout.labels;
    else
        ground_truth = [];
    end
    
    % separate the configuration file
    config = model.config;
    
    % initialize a new labelling
    labelling = zeros(Gout.w, Gout.l);
    
    % Encode the graph in a different way (patterns and labels 
    % have size = 1)
    [patterns, labels] = encodeTrainingData({Gout});
    patterns = patterns{1};
    labels = labels{1};
    
    % separate the weights between unaries and pairwises
    [w_unary, w_pairwise] = get_weights(model.w, config);
    

    % 1) Get the unary features

    % Encode the unary features
    unary_features = get_features(patterns, labels); % Take only the unary features and the bias, and ignore the pairwise
    unary_features = reshape(unary_features, size(unary_features, 1), 2, config.size_w_u / 2);
    unary_features = permute(unary_features, [3 1 2]);
    
    % 2) Get the pairwise features
    pairwise_features = patterns{2};
       
    % 4) Compute the scores and multiply them by -1 to get the unary energy
    unary_potentials = zeros(size(unary_features, 2), 2);
    
    % scores for veins
    unary_potentials(:,1) = (w_unary(1, :) * unary_features(:,:,1))'; % veins
    % scores for arteries
    unary_potentials(:,2) = (w_unary(2, :) * unary_features(:,:,2))'; % arteries    

    % retrieve the number of elements
    N = size(pairwise_features, 1);
    % reshape the pairwise features so we can compute the pairwise
    % potentials easily
    pairwise_features = reshape(pairwise_features, N * N, config.size_w_p);
    % take the product between the pairwise features and their weights
    pairwise_potentials = pairwise_features * w_pairwise;
    % reshape the pairwise potentials to their original shape
    pairwise_potentials = reshape(pairwise_potentials, N, N);
    
    % 5) Optimize the energy  
    yhat = optimize_with_maxflow(sparse(unary_potentials), sparse(pairwise_potentials));
    
    
    % assign labels to the graph and prepare the image
    for i = 1 : length(Gout.node)
        % assign to the node in the graph
        Gout.node(i).label = yhat(i);
        % and also to the labelling
        labelling(Gout.node(i).idx) = yhat(i);
    end

end