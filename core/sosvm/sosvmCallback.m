
function [model, config, state] = sosvmCallback(config, training_data)
% sosvmCallback Configure the SOSVM and call it to learn the model
% [model, config, state] = sosvmCallback(config, trainingdata)
% OUTPUT: model: learned model
%         config: configuration structure
%         state: last state
% INPUT: config: configuration structure
%        trainingdata: training data (is a cell array of graphs)
    
    % Assign functions
    config.findMostViolatedMarginFn = @findMostViolatedConstraint;
    config.lossFn = @lossComputing;
    config.psiFn = @featureComputing;
    
    % % Normalize the value of C by the number of pixels
    %config.C = config.C / something; % We need to do this because of our SOSVM implementation
    
    % Encode the training data
    fprintf('Encoding training data\n');
    [patterns, labels] = encodeTrainingData(training_data);
    fprintf('Training data encoded\n');
    
    % Train the SOSVM
    [model, config, state] = sosvm2(config, patterns, labels);
    
end

% callback to the most violated constraint
function [yhat] = findMostViolatedConstraint(param, model, x, y)
    yhat = constraintCB(param, model, x, y);
end

% callback to the loss function
function [loss] = lossComputing(param, y, tildey)
    loss = lossCB(param, y, tildey);
end

% callback to the feature map function
function [psi] = featureComputing(sparm, x, y)
    [psi] = featureCB(sparm, x, y);
end