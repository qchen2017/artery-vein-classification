
function [A, T] = get_crf_energies_from_graph(Gout, w_u, w_p)

    % initialize matrix of unary and pairwise potentials
    % A (smoothness term): pairwise
    % T (data term): unary
    A = sparse(length(Gout.node), length(Gout.node));
    T = sparse(length(Gout.node), 2);
    
    % for each node in Gout
    for i = 1 : length(Gout.node)
        T(i, 1) = w_u(:,1)' * Gout.node(i).features;
        T(i, 2) = w_u(:,2)' * Gout.node(i).features;
    end
    
    % for each link
    for i = 1 : length(Gout.link)
        A(Gout.link(i).n1, Gout.link(i).n2) = w_p * Gout.link(i).features';
        A(Gout.link(i).n2, Gout.link(i).n1) = A(Gout.link(i).n1, Gout.link(i).n2);
    end

end