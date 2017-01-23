
function Gmatlab = generate_matlab_graph(Gin)

    adjacency_matrix = logical(length(Gin.node));
    
    for i = 1 : length(Gin.node)
        for j = 1 : length(Gin.node(i).conn)
            adjacency_matrix(i, Gin.node(i).conn(j)) = true;
            adjacency_matrix(Gin.node(i).conn(j), i) = true;
        end
    end
    
    Gmatlab = graph(adjacency_matrix);

end