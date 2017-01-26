function [model, config, state] = sosvm2(config, patterns, labels, oldstate)
% sosvm Learn a model using a SOSVM
% [model, config, state] = sosvm(config, patterns, labels, oldstate)
% OUTPUT: model: learned model
%         config: configuration structure, updated with learning
%         information
%         state: last state
% INPUT: config: configuration structure
%        patterns: cell-array with the training data
%        labels: training labels
%        oldstate: previous state
    
    state = bundler(); % initialize state
    
    % lambda = 1 / C
    state.lambda = 1 ./ (config.C);
    
    % Include additional parameters
    config.convergenceThreshold = 0.001;
    config.formulationType = 'margin';
    config.maxIterations = 40;
    config.minIterations = 10;
    
    % Initially, w has 0s
    config.w = zeros(config.sizePsi, 1);
    state.w = config.w;
    
    % Set the positivity contraints
    for i = 1 : length(config.posindx);
        % add in the positivity (submodularity) constraints
        phi = zeros(size(config.w));
        phi(config.posindx(i)) = 1;
        % call bundler with hard constraint argument
        state = bundler(state, phi, 0, false);
    end
    % Copy the weights learned from bundler with hard positivity contraints
    % (if we don't have any, we will copy an array of 0s)
    config.w = state.w;
    model.w = state.w;

    % if there is an old state...
    if (exist('oldstate','var'))
        for i=1:length(oldstate.b)
            if(oldstate.softVariables(i))
                state = bundler(state,oldstate.a(:,i),oldstate.b(i));
            end
        end
    end
    
    numIterations = 0;
    bestPrimalObjective = Inf;
    
    % repeat until convergence
    while (((bestPrimalObjective - state.dualObjective) / state.dualObjective > config.convergenceThreshold ...
            || config.minIterations>0 ) && numIterations < config.maxIterations)
        
        numIterations = numIterations + 1;
        config.minIterations = config.minIterations - 1;

        % Compute using one slack formulation of SOSVM with margin
        % rescaling
        [phi, b] = computeOneslackMarginConstraint(config, model, patterns, labels);

        % Get the primal objective (we use it to check convergence) (as psi
        % is higher than b - dot(state.w, phi), then this is an upper bound
        primalobjective = (state.lambda / 2) * (state.w' * state.w) + b - dot(state.w, phi);
        if ((primalobjective < bestPrimalObjective) || true)
            bestPrimalObjective = primalobjective;
            bestState = state;
        end
        
        % Print values on creen
        fprintf([' %d primal objective: %f, best primal: %f, dual objective: %f, gap: %f\n'], ...
                   numIterations, primalobjective, bestPrimalObjective, state.dualObjective, ...
                   (bestPrimalObjective - state.dualObjective) / state.dualObjective);

        % Call bundler again to retrieve the state
        state = bundler(state, phi, b);
        config.w = state.w;
        model.w = state.w;
        
        % Assert positivity constraint
        assertPositivity(config, model);
    
    end
    
    % Return values
    config.w = bestState.w;
    model.w = bestState.w;

end


function assertPositivity(param, model)
% ASSERTPOSITIVITY This function checks if the positivity constraints are
% violated or not
    if ~isempty(param.posindx)
        assert(sum(model.w(param.posindx) >= ones(length(param.posindx), 1) * -1.0e-6) == length(param.posindx), 'Positivity contraint violated by ');
    end
end


function [phi, b] = computeOneslackMarginConstraint(config, model, X, Y)
% COMPUTEONESCLACKMARGINCONSTRAINT This function returns phi and b as
% obtained by cutting planes

    phi = 0;
    b = 0;
    
    % For each pattern
    for i = 1 : length(X);
        
        % Our separation oracle...
        [tildeY] = config.findMostViolatedMarginFn(config, model, X{i}, Y{i});
        % Our loss (Delta)
        delta = config.lossFn(config, Y{i}, tildeY);
        % Our loss as the differences...
        deltaPsi = config.psiFn(config, X{i}, Y{i}) - config.psiFn(config, X{i}, tildeY);
        % Accumulate b and phi
        if (dot(model.w,deltaPsi) < delta)
            b = b + delta;
            phi = phi + deltaPsi;
        end 

    end
    
end
