function [Gout] = solveIddleLinksAndNodes(Gin)

    % Initialize the new graph
    Gout.node = [];
    Gout.link = [];
    Gout.w = Gin.w;
    Gout.l = Gin.l;
    Gout.onh_perimeter = Gin.onh_perimeter;    
    
    % IDENTIFY IDDLE NODES
    % nodesAliases is used to identify the equivalent index of each new
    % node
    % ----> nodes=-1 are fake nodes
    % ----> nodes=0 were not named yet
    nodesAliases = (extractfield(Gin.node, 'numLinks')==0);
    % Add nodes without conflicts
    [Gout, nodesAliases] = addNodesWithoutConflicts(Gin, Gout, nodesAliases);
    
    % IDENTIFY IDDLE LINKS
    % Initialize the aliases
    last_link_alias = 0;
    last_node_alias = max(nodesAliases(:));
    % Find onh_perimeter foreground pixels
    onh_pixels = find(Gin.onh_perimeter);
    
    % For each link
    for i = 1 : length(Gin.link)
        % If the link goes from one node to a different one
        if (Gin.link(i).n1 == -1) && (Gin.link(i).n2 == -1)
            
            % Identify which points of the segment are connected to the ONH
            % perimeter (or at least are close)
            are_member = ismember(Gin.link(i).point, onh_pixels);
            
            % Let's create a new node
            last_node_alias 
            
            
            
        else
            % Add the link to the graph
            new_alias = last_link_alias + 1;
            % Create the new link
            if (Gin.link(i).n1==-1)
                new_link.n1 = -1;
            else
                new_link.n1 = nodesAliases(Gin.link(i).n1);
            end
            if (Gin.link(i).n2==-1)
                new_link.n2 = -1;
            else
                new_link.n2 = nodesAliases(Gin.link(i).n2);
            end
            new_link.point = Gin.link(i).point;
            % Add the link and update the graph
            [Gout] = addNewLinkToGraph(Gout, new_alias, new_link);
            % Update the last link alias
            last_link_alias = new_alias;
        end
    end


end