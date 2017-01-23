function [Gout, nodesAliases, visited] = getLiveLinks(currentNode, visited, Gin, Gout, nodesAliases)

    % if the node was never visited...
    if (~visited(currentNode))
        
        % mark the node as visited
        visited(currentNode) = 1;
        
        % check if it belongs to the list of fake nodes
        if (nodesAliases(currentNode)~=-1)
            
            % it is not a fake node, check if it exists
            if (nodesAliases(currentNode)==0)
                % if it don't exist, add it to the graph
                % set the alias
                nodesAliases(currentNode) = length(Gout.node) + 1;
                nodeNewName = nodesAliases(currentNode);
                % add the node to the graph structure
                Gout.node(nodeNewName).idx = Gin.node(currentNode).idx;
                Gout.node(nodeNewName).links = [];
                Gout.node(nodeNewName).conn = [];
                Gout.node(nodeNewName).comx = Gin.node(currentNode).comx;
                Gout.node(nodeNewName).comy = Gin.node(currentNode).comy;
                Gout.node(nodeNewName).numLinks = 0;
            else
                nodeNewName = nodesAliases(currentNode);
            end
            
            % recover the links
            link_idxs = Gin.node(currentNode).links;
            links = Gin.link(link_idxs);
            % for each link
            for l = 1 : length(links)
                
                % if the link goes/comes to/from a node alive or a non active node...
                if ((links(l).n2==-1) || (links(l).n1==-1) || ( (nodesAliases(links(l).n1)~=-1) && (nodesAliases(links(l).n2)~=-1) ))
                    
                    % get the link index
                    linkNewName = length(Gout.link) + 1;
                    
                    % add the link
                    new_link.n1 = nodeNewName;
                    if (links(l).n2==-1)
                        endNode = -1; % set the end node as -1
                    else
                        if (links(l).n2==currentNode)
                            endNode = nodesAliases(links(l).n1); 
                            links(l).point = fliplr(links(l).point);
                        else
                            endNode = nodesAliases(links(l).n2); % set the end node to the (alias of the) node
                        end
                    end
                    new_link.n2 = endNode;
                    new_link.point = links(l).point;
                    [Gout] = addNewLinkToGraph(Gout, linkNewName, new_link);
                    
                end
                
            end
                        
        end
        
        % analyze the neighbors of the current node
        nbs = Gin.node(currentNode).conn(Gin.node(currentNode).conn~=-1);

        for i = 1 : length(nbs)
            if (~visited(nbs(i))) % if they were not visited
                [Gout, nodesAliases, visited] = getLiveLinks(nbs(i), visited, Gin, Gout, nodesAliases);
            end
        end
                    
    end
        
    
end