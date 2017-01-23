
function [Gout] = solveSelfCycles(Gin)

    % Initialize the new graph
    Gout.node = [];
    Gout.link = [];
    Gout.w = Gin.w;
    Gout.l = Gin.l;
    Gout.onh_perimeter = Gin.onh_perimeter;    
    
    % Copy all the nodes
    for i = 1 : length(Gin.node)
        Gout.node(i).idx = Gin.node(i).idx;
        Gout.node(i).conn = [];
        Gout.node(i).links = [];
        Gout.node(i).comx = Gin.node(i).comx;
        Gout.node(i).comy = Gin.node(i).comy;
        Gout.node(i).numLinks = 0;
    end
    
    % Initialize the aliases
    last_alias = 0;
    % For each link
    for i = 1 : length(Gin.link)
        % If the link goes from one node to a different one
        if (Gin.link(i).n1 ~= Gin.link(i).n2)
            % Add the link to the graph
            new_alias = last_alias + 1;
            [Gout] = addNewLinkToGraph(Gout, new_alias, Gin.link(i));
            last_alias = new_alias;
        
        end
    end


end