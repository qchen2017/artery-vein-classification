
function [Gout, nodesAliases] = addNodesWithoutConflicts(Gin, Gout, nodesAliases)

    nonDeadNodes = find(nodesAliases==0);
    firstAlias = 0;
    for i = 1 :length(nonDeadNodes)
        % generate the alias
        nodeNewName = firstAlias + i;
        % add the new node
        Gout.node(nodeNewName).idx = Gin.node(nonDeadNodes(i)).idx;
        Gout.node(nodeNewName).links = [];
        Gout.node(nodeNewName).conn = [];
        Gout.node(nodeNewName).comx = Gin.node(nonDeadNodes(i)).comx;
        Gout.node(nodeNewName).comy = Gin.node(nonDeadNodes(i)).comy;
        Gout.node(nodeNewName).numLinks = 0;
        % save the alias
        nodesAliases(nonDeadNodes(i)) = nodeNewName;
    end
    
    % Reconstruct the graph including only edges and vertices that are
    % clean in the original version
    visitedNonDeadNodes = [];
    for i = 1 : length(nonDeadNodes)       
        % mark the node as visited
        visitedNonDeadNodes = [visitedNonDeadNodes nonDeadNodes(i)];
        % get the old id of the current node
        currentNode = nonDeadNodes(i);
        % recover the links
        link_idxs = Gin.node(currentNode).links;
        links = Gin.link(link_idxs);
        % for each link
        for l = 1 : length(links)
            % identify the node different than the current one
            otherNode = links(l).n1;
            if (otherNode == currentNode)
                otherNode = links(l).n2;
            end
            % if the link goes/comes to/from a non active node or both extremes are active nodes or the active end node is not visited...
            if (links(l).n2==-1) || (links(l).n1==-1) || ( (nodesAliases(links(l).n1)~=-1) && (nodesAliases(links(l).n2)~=-1) && ~ismember(otherNode,visitedNonDeadNodes))
                
                % get the link index
                linkNewName = length(Gout.link) + 1;
                % prepare the new link
                new_link.n1 = nodesAliases(currentNode);% = nodesAliases(currentNode);;
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
                % add the new link
                [Gout] = addNewLinkToGraph(Gout, linkNewName, new_link);

            end
        end
    end
    
end