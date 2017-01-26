
function [classification, flow] = optimize_with_maxflow(unary_potentials, pairwise_potentials)

    % initialize matrix of unary and pairwise potentials
    % A (smoothness term): pairwise
    % T (data term): unary
    A = sparse(pairwise_potentials);
    T = sparse(unary_potentials);
    
    % Optimize with maxflow
    [flow, classification] = maxflow(A, T);
    classification = 2 * classification - 1;

end