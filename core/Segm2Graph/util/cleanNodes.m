
function Gout = cleanNodes(Gin)

    % Initialize the new graph
    Gout.node = [];
    Gout.link = [];
    Gout.w = Gin.w;
    Gout.l = Gin.l;
    Gout.onh_perimeter = Gin.onh_perimeter;
    Gout.segm = Gin.segm;
    
    % nodesAliases is used to identify the equivalent index of each new
    % node
    % ----> nodes=-1 are fake nodes
    % ----> nodes=0 were not named yet
    nodesAliases = ((extractfield(Gin.node, 'numLinks')<=2) .* ~cell2mat(extractfield(Gin.node, 'isroot'))) * -1; 
    
    % Add nodes without conflicts
    [Gout, nodesAliases] = addNodesWithoutConflicts(Gin, Gout, nodesAliases);
    
    % Now we need to identify the indexes of the nodes to be removed from
    % the original graph: vessel crossings
    deadNodes = find(nodesAliases==-1);
    
    % We iterate for each node to be removed until the array of dead nodes
    % is empty
    iterator = 1;
    last_node_idx = sum(find(nodesAliases~=-1)>0);
    while ~(isempty(deadNodes))

        % get one dead node
        i = deadNodes(iterator);
        
        % recover the new link (link without the dead node i)
        [newLink, deadNodes] = getRefinedLinks(Gin, nodesAliases, i, deadNodes);
        
        if (isstruct(newLink))
        
            % generate the id of the new link
            newLink_idx = length(Gout.link) + 1;

            % swap nodes every time a node source is -1, so every node begins
            % with an active node
            if (newLink.n1 == -1)
                newLink.n1 = newLink.n2;
                newLink.n2 = -1;
                % flip the array of points so the points of the link begins in
                % an active node
                newLink.point = fliplr(newLink.point);
            end
            % update the id of the n1 node according to nodes aliases in the
            % new graph
            if (newLink.n1~=-1)

                % update the id to the new alias
                newLink.n1 = nodesAliases(newLink.n1);
                % update the id of the n2 node
                if (newLink.n2~=-1)
                    % update the id to the new alias
                    newLink.n2 = nodesAliases(newLink.n2);
                end
                % add the new link to the graph
                [Gout] = addNewLinkToGraph(Gout, newLink_idx, newLink);

            end
            
        else
            
            % add the dead node because actually it is not one of them
            last_node_idx = last_node_idx + 1;
            
            % add the new node
            Gout.node(last_node_idx).idx = Gin.node(i).idx;
            Gout.node(last_node_idx).links = [];
            Gout.node(last_node_idx).conn = [];
            Gout.node(last_node_idx).comx = Gin.node(i).comx;
            Gout.node(last_node_idx).comy = Gin.node(i).comy;
            Gout.node(last_node_idx).numLinks = 0;
            Gout.node(last_node_idx).isroot = Gin.node(i).isroot;
            % save the alias
            nodesAliases(i) = last_node_idx;
            
            for non_dead_link_idx = 1 : length(Gin.node(i).links)
                new_link = Gin.link(Gin.node(i).links(non_dead_link_idx));
                if (new_link.n1~=-1)
                    new_link.n1 = nodesAliases(new_link.n1);
                end
                if (new_link.n2~=-1)
                    new_link.n2 = nodesAliases(new_link.n2);
                end
                [Gout] = addNewLinkToGraph(Gout, length(Gout.link)+1, new_link);
            end
            
            
        end
        
        % update the iterator
        iterator = 1;
        
    end

end



% Obtain clean links (without nodes to be removed)
function [link_out, remainingDeadNodes] = getRefinedLinks(Gin, nodesAliases, currentNode, deadNodes)

    % get nodes connected to the current node
    nodesConnected = Gin.node(currentNode).conn;
    
    % if two nodes are connected to the node...
    if (length(nodesConnected)==2)
        
        % identify each node so it is possible to identify each case
        n1 = nodesConnected(1);
        n2 = nodesConnected(2);

        if ((n1==-1) && (n2==-1))
            
            new_link = [];
            
        
        % case A: both nodes are alive
        elseif (((n1==-1) || (nodesAliases(n1)~=-1)) && ((n2==-1) || (nodesAliases(n2)~=-1)))

            % get l1
            l1 = Gin.link(Gin.node(currentNode).links(nodesConnected==n1));
            if (length(l1)>1)
                l1 = l1(1);
            end
            if (l1.n2~=currentNode) 
                % flip the points so the dead node is at the end
                l1.point = fliplr(l1.point);
            end
            % get l2
            l2 = Gin.link(Gin.node(currentNode).links(nodesConnected==n2));
            if (length(l2)>1)
                l2 = l2(1);
            end
            if (l2.n1~=currentNode)
                % flip the points so the dead node is at the beginning
                l2.point = fliplr(l2.point);
            end
            % the link l is such that l=[l1,N,l2]
            new_link.n1 = n1;
            new_link.n2 = n2;
            new_link.point = [l1.point, l2.point];

        % case B: n1 alive, n2 dead
        elseif (((n1==-1) || (nodesAliases(n1)~=-1)) && ((n2~=-1) && (nodesAliases(n2)==-1)))

            % get l1
            l1 = Gin.link(Gin.node(currentNode).links(nodesConnected==n1));
            if (l1.n2~=currentNode) 
                % flip the points so the dead node is at the end
                l1.point = fliplr(l1.point);
            end
            % get l2=l(n2)
            [l2, deadNodes] = rebuildEdge(Gin, nodesAliases, currentNode, n2, deadNodes);
            % the link l is such that l=[l1,N,l(n2)]
            new_link.n1 = n1;
            if (l2.n1==currentNode)
                new_link.n2 = l2.n2;
            else
                new_link.n2 = l2.n1;
                % flip the points so the dead node is at the beginning
                l2.point = fliplr(l2.point);
            end
            new_link.point = [l1.point, l2.point];

        % case C: n1 dead, n2 alive
        elseif (((n1~=-1) && (nodesAliases(n1)==-1)) && ((n2==-1) || (nodesAliases(n2)~=-1)))

            % swap n1 y n2 so we can use exactly the same code than in case
            % B
            aux = n2;
            n2 = n1;
            n1 = aux;
            
            % get l1
            l1 = Gin.link(Gin.node(currentNode).links(nodesConnected==n1));
            if (l1.n2~=currentNode) 
                % flip the points so the dead node is at the end
                l1.point = fliplr(l1.point);
            end
            % get l2=l(n2)
            [l2, deadNodes] = rebuildEdge(Gin, nodesAliases, currentNode, n2, deadNodes);
            % the link l is such that l=[l1,N,l(n2)]
            new_link.n1 = n1;
            if (l2.n1==currentNode)
                new_link.n2 = l2.n2;
            else
                new_link.n2 = l2.n1;
                % flip the points so the dead node is at the beginning
                l2.point = fliplr(l2.point);
            end
            new_link.point = [l1.point, l2.point];
            
            new_link.point = fliplr(new_link.point);
            aux = new_link.n1;
            new_link.n1 = new_link.n2;
            new_link.n2 = aux;

        % case D: n1 y n2 dead
        elseif (((n1~=-1) && (nodesAliases(n1)==-1)) && ((n2~=-1) && (nodesAliases(n2)==-1)))

            % get l1=l(n1)
            [l1, deadNodes] = rebuildEdge(Gin, nodesAliases, currentNode, n1, deadNodes);
            if (l1.n1==currentNode)
                new_link.n1 = l1.n2;
                l1.point = fliplr(l1.point);
            else
                new_link.n1 = l1.n1;
            end
            % get l2 =l(n2)
            [l2, deadNodes] = rebuildEdge(Gin, nodesAliases, currentNode, n2, deadNodes);
            if (l2.n1==currentNode)
                new_link.n2 = l1.n2;
            else
                new_link.n2 = l1.n1;
                l2.point = fliplr(l2.point);
            end
            % the link l is such that l=[l(n1),N,l(n2)]
            new_link.point = [l1.point, l2.point];

        end
        
    % if only one node is connected
    else
        
        if isempty(nodesConnected)
            
            new_link = -1;
            
        else
            
            % identify the node
            nn = nodesConnected(1); 

            % case E: if the node connected to the current one must die
            if ((nn~=-1) && (nodesAliases(nn)==-1))

                % get l(nn)
                [l, deadNodes] = rebuildEdge(Gin, nodesAliases, currentNode, nn, deadNodes);

            % case F: if the node connected to the current one is alive
            else

                % get the link
                l = Gin.link(Gin.node(currentNode).links(nodesConnected==nn));

            end

            % include the link
            if (l.n2~=currentNode)
                new_link.n1 = l.n2;
                new_link.n2 = -1;
                l.point = fliplr(l.point);
            else
                new_link.n1 = l.n1;
                new_link.n2 = -1;
            end
            new_link.point = [l.point];
            
        end
        
        
    end
    
    % assign the new link
    link_out = new_link;
    
    % remove the current node from the list of dead nodes
    remainingDeadNodes = deadNodes(deadNodes~=currentNode);
    
    
end



function [link_out, remainingDeadNodes] = rebuildEdge(Gin, nodesAliases, previousNode, currentNode, deadNodes)

    % identifico el link de ingreso
    nodesConnected = Gin.node(currentNode).conn;
    %link_in_idx = find(nodesConnected==previousNode);
    link_in = Gin.link(Gin.node(currentNode).links(nodesConnected==previousNode));
    if (link_in.n2~=currentNode)
        link_in.point = fliplr(link_in.point);
        link_in.n1 = link_in.n2;
        link_in.n2 = currentNode;
    end
        
    % si el nodo es un nodo que deba morir
    if (nodesAliases(currentNode)==-1)
        
        % identifico el próximo nodo
        nextNode = nodesConnected(nodesConnected~=previousNode);
        
        % inicializo los nodos muertos restantes con la actual lista de
        % nodos muertos
        remainingDeadNodes = deadNodes;
        
        % obtengo el resto del link
        if (~isempty(nextNode) && (nextNode~=-1)) % si su próximo nodo existe o no es el vacio...
            
            % obtengo el resto del arco
            [link_out, remainingDeadNodes] = rebuildEdge(Gin, nodesAliases, currentNode, nextNode, deadNodes);

            % construyo el nuevo link
            new_link.n1 = (previousNode);            
            if (link_out.n1~=currentNode)
                new_link.n2 = (link_out.n1);
                new_link.point = fliplr(link_out.point);
            else
                new_link.n2 = (link_out.n2);
            end
            new_link.point = [link_in.point, link_out.point];
            
        elseif isempty(nextNode)  % si su próximo nodo no existe, es porque es un nodo de un solo vecino;
            
            % agrego el link
            new_link.n1 = (previousNode);
            new_link.n2 = -1;
            new_link.point = [link_in.point];
            
        else % si el próximo nodo es -1
            
            % obtengo el link hacia el final
            %link_out_idx = find(nodesConnected==nextNode);
            link_out = Gin.link(Gin.node(currentNode).links(nodesConnected==nextNode));
            if (link_out.n1~=currentNode)
                link_out.point = fliplr(link_out.point);
                link_out.n1 = currentNode;
                link_out.n2 = -1;
            end
            
            % agrego el link
            new_link.n1 = (previousNode);
            new_link.n2 = -1;
            new_link.point = [link_in.point, link_out.point];
            
        end
        
        % eliminar el nodo muerto de la lista de links muertos 
        remainingDeadNodes(remainingDeadNodes==currentNode) = [];
            
        % retorno el link
        link_out = new_link;
        
    else % si no debe morir
        
        % devuelvo el link anterior
        link_out = link_in;
        remainingDeadNodes = deadNodes;
        
    end
    
end