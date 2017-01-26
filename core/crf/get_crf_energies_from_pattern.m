
function [A, T] = get_crf_energies_from_pattern(X, w_u, w_p)

    % initialize matrix of unary and pairwise potentials
    % A (smoothness term): pairwise
    % T (data term): unary
    
    % T will be the unary potentials
    T = sparse(X{1} * w_u);
    
    % A will be the pairwise potentials
    A = sparse(X{2} * w_p);

end