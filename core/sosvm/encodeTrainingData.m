
function [patterns, labels] = encodeTrainingData(training_data) 

    % Preallocate arrays for patterns and labels
    patterns = cell(size(training_data));
    labels = cell(size(training_data));
    
    % For each graph in training_data
    for i = 1 : length(training_data)
        
        % Encode unary features (including a bias term) and labels
        [unary_features, y] = encode_unary_features_and_labels(training_data{i});
        
        % Encode pairwise features
        [pairwise_features, pairwise_idx] = encode_pairwise_features(training_data{i});
        
        % Assign current pattern
        patterns{i} = {unary_features, pairwise_features, pairwise_idx, training_data{i}};
        % And current labels
        labels{i} = y;
        
    end

end




function [unary_features, labels] = encode_unary_features_and_labels(Gout)
    
    % Our matrix of unary features will have the size of the number of
    % features x number of segments (nodes)
    unary_features = ones(length(Gout.node), Gout.properties.unary_dim + 1) * 1;

    % Initialize the array of labels
    labels = zeros(length(Gout.node), 1);
    
    % For each node in the graph
    for i = 1 : length(Gout.node)
        
        % Retrieve the unary features and assign to the matrix
        unary_features(i, 1:end-1) = Gout.node(i).features;
        % Get current label
        labels(i) = Gout.node(i).label;
        if Gout.node(i).label == 0
            labels(i) = -1; % by default, we assume that segments labeled with 0 are veins
        end
        
    end
    
    % CAREFUL WITH THIS
    %unary_features = cat(2, labels, ones(size(labels)));
    
    % Standardize unary features
    unary_features(:,1:end-1) = standardizeCols(unary_features(:,1:end-1), mean(unary_features(:,1:end-1)), std(unary_features(:,1:end-1)));

end


function [pairwise_features, pairwise_idx] = encode_pairwise_features(Gout)

    % Our matrix of pairwise features will be like an adjacency matrix,
    % with |V| x |V| x |pairwise_feature| size
    pairwise_features = zeros(length(Gout.node), length(Gout.node), Gout.properties.pairwise_dim);
    % We will also have a similar matrix to preserve pairwise idx
    pairwise_idx = zeros(length(Gout.node), length(Gout.node));
    
    % For each link in the graph
    for i = 1 : length(Gout.link)
        
       % Assign features
       pairwise_features(Gout.link(i).n1, Gout.link(i).n2, :) = Gout.link(i).features;
       pairwise_features(Gout.link(i).n2, Gout.link(i).n1, :) = Gout.link(i).features;
        
       % Copy the idx of the link
       pairwise_idx(Gout.link(i).n1, Gout.link(i).n2) = i;
       pairwise_idx(Gout.link(i).n2, Gout.link(i).n1) = i;
        
    end

end

