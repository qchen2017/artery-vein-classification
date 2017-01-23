
function [Gout] = solveProximityConflicts(Gin)

    % Initialize the new graph
    Gout.node = [];
    Gout.link = [];
    Gout.w = Gin.w;
    Gout.l = Gin.l;
    Gout.onh_perimeter = Gin.onh_perimeter;

    % nodesAliases is used to identify the equivalent index of each new
    % node
    % ----> nodes=-1 are nodes with proximity conflicts
    % ----> nodes=0 were not named yet
    node_idxs_four_neighbors = find(extractfield(Gin.node, 'numLinks')==4); 
    nodesAliases = zeros(size(Gin.node));
    for i = 1 : length(node_idxs_four_neighbors)
        % get the index of the potential node with proximity conflict
        candidate_idx = node_idxs_four_neighbors(i);
        % get nodes connected to the current candidate
        nodesConnected = Gin.node(candidate_idx).conn;
        % if there are two links to the same node, the node must be marked
        % to be removed
        nodesAliases(candidate_idx) = ...
            (length(unique(nodesConnected(nodesConnected~=-1))) ...
                < length(nodesConnected(nodesConnected~=-1))) * -1;
    end
    
    % identificar las raíces de los nodos con conflictos de proximidad
    node_idxs_three_neighbors = find(extractfield(Gin.node, 'numLinks')==3);
    rootNodes = [];
    count = 1;
    for i = 1 : length(node_idxs_three_neighbors)
        % get the index of the candidate
        candidate_idx = node_idxs_three_neighbors(i);
        % get connected nodes
        nodesConnected = Gin.node(candidate_idx).conn;
        % if there are two links to the same node, the node must be marked
        % to be removed
        if (length(unique(nodesConnected(nodesConnected~=-1))) < length(nodesConnected(nodesConnected~=-1)))
            
            n = mode(Gin.node(candidate_idx).conn);
            % si su vecino está en la lista de rootNodes, son nodos
            % próximos que hay que eliminar porque son de otro tipo
            if (~ismember(n,rootNodes))
                rootNodes(count) = candidate_idx;
                count = count + 1;
            else
                rootNodes(rootNodes==n)=[];
                count = count - 1;
            end
            
        end
    end
    
    % obtener el grafo considerando nodos y links que no tengan conflictos
    % de proximidad
    %visited = zeros(size(Gin.node));
    conflictiveNodes = nodesAliases;
    conflictiveNodes(rootNodes) = 0;
    [Gout, nodesAliases] = addNodesWithoutConflicts(Gin, Gout, conflictiveNodes);
    %[Gout, nodesAliases, ~] = getLiveLinks(1, visited, Gin, Gout, nodesAliases);
    
    % resolver para cada raíz sus conflictos de proximidad
    for i = 1 : length(rootNodes)
        
        % get the root node
        n0 = rootNodes(i);
        % identifico los nodos conectados repetidos
        nodesConnected = Gin.node(n0).conn;
        n1 = mode(nodesConnected);
        
        % identifico los links conectados a n1
        next_links = find(nodesConnected==n1);
        link_begin_1 = swapLink(Gin.link(Gin.node(n0).links(next_links(1))), n1);
        link_begin_2 = swapLink(Gin.link(Gin.node(n0).links(next_links(2))), n1);
        
        % obtengo el resto del link
        [link_end_1, link_end_2] = solveLinksWithProximityConflicts(Gin, n0, n1, nodesAliases);
        
        n0 = nodesAliases(n0);
        
        % armo el link nuevo 1
        newlink_1.n1 = n0;
        if (link_end_1.n2==-1)
            newlink_1.n2 = -1;
        else
            newlink_1.n2 = nodesAliases(link_end_1.n2);
        end
        newlink_1.point = [link_begin_1.point link_end_1.point];
        % armo el link nuevo 2
        newlink_2.n1 = n0;
        if (link_end_2.n2==-1)
            newlink_2.n2 = -1;
        else
            newlink_2.n2 = nodesAliases(link_end_2.n2);
        end
        newlink_2.point = [link_begin_2.point link_end_2.point];
        
        % get the new link index
        newlink_idx_1 = length(Gout.link) + 1;
        newlink_idx_2 = newlink_idx_1 + 1;
        
        % add links to the graph structure
        [Gout] = addNewLinkToGraph(Gout, newlink_idx_1, newlink_1);
        [Gout] = addNewLinkToGraph(Gout, newlink_idx_2, newlink_2);
        
    end



end




function [swaped_link] = swapLink(link, n2)
    if (link.n2~=n2)
        swaped_link.n1 = link.n2;
        swaped_link.n2 = n2;
        swaped_link.point = fliplr(link.point);
    else
        swaped_link = link;
    end   
end




function [link1, link2] = solveLinksWithProximityConflicts(Gin, previousNode, currentNode, nodesAliases)


    % si el nodo actual tiene proximity conflict
    if (nodesAliases(currentNode)==-1)

        % identifico los nodos conectados repetidos
        nodesConnected = unique(Gin.node(currentNode).conn);
        nextNode = nodesConnected(nodesConnected~=previousNode);
        
        % si los próximos nodos son 2, cortamos:
        % son dos links que se separan
        if (length(nextNode)==2)
            
            % identifico los 2 nodos destino
            next_node_1 = nextNode(1);
            next_node_2 = nextNode(2);
            
            link_begin_1 = swapLink(Gin.link(Gin.node(currentNode).links(find(Gin.node(currentNode).conn==next_node_1))), next_node_1);
            link_begin_2 = swapLink(Gin.link(Gin.node(currentNode).links(find(Gin.node(currentNode).conn==next_node_2))), next_node_2);
            
            % el principio es el fin en este caso
            link1 = link_begin_1;
            link2 = link_begin_2;
          
        % si hay uno solo, es porque hay que seguir
        else
            
            % identifico los links conectados a nextNode
            next_links = find(Gin.node(currentNode).conn==nextNode);
            link_begin_1 = swapLink(Gin.link(Gin.node(currentNode).links(next_links(1))), nextNode);
            link_begin_2 = swapLink(Gin.link(Gin.node(currentNode).links(next_links(2))), nextNode);
            
            % si el próximo nodo es -1, corto
            if (nextNode==-1)
                
                % el principio es el fin en este caso
                link1 = link_begin_1;
                %link1.point = [Gin.node(currentNode).idx' link1.point];
                link2 = link_begin_2;
                %link2.point = [Gin.node(currentNode).idx' link2.point];
                
            % si no, sigo    
            else
            
                % obtengo los links posteriores
                [link_end_1, link_end_2] = solveLinksWithProximityConflicts(Gin, currentNode, nextNode, nodesAliases);
                % armo los links nuevos
                link1.n2 = link_end_1.n2;
                %link1.point = [Gin.node(currentNode).idx' link_begin_1.point link_end_1.point];
                link1.point = [link_begin_1.point link_end_1.point];
                link2.n2 = link_end_2.n2;
                %link2.point = [Gin.node(currentNode).idx' link_begin_2.point link_end_2.point];
                link2.point = [link_begin_2.point link_end_2.point];
                
            end
        
        end

    % si no tiene proximity conflict, es porque es un root de cierre
    else
        
        % identifico el nodo conectado que no está repetido
        % (ya sé que el previo es un nodo de conflicto)
        nodesConnected = unique(Gin.node(currentNode).conn);
        nextNode = nodesConnected(nodesConnected~=previousNode);
        
        if length(nextNode)==1
        
            link1 = swapLink(Gin.link(Gin.node(currentNode).links(find(Gin.node(currentNode).conn==nextNode))), nextNode);
            link2.n2 = -1;
            link2.point = [];
            
        elseif isempty(nextNode)
            
            link1.n2 = -1;
            link1.point = [];
            link2.n2 = -1;
            link2.point = [];
            
        else
            
            if sum(Gin.node(currentNode).conn == -1)==1
                
                link1 = swapLink(Gin.link(Gin.node(currentNode).links(find(Gin.node(currentNode).conn==-1))), -1);
                link2.n2 = currentNode;
                link2.point = [];
                
            elseif sum(Gin.node(currentNode).conn == -1)==20
                
                links = Gin.link(Gin.node(currentNode).links(find(Gin.node(currentNode).conn==-1)));
                link1 = swapLink(links(1), -1);
                link2 = swapLink(links(2), -1);
                
            else
                
                nodesConnected = Gin.node(currentNode).conn;
                nodesConnected(nodesConnected == previousNode) = [];
                
                if length(unique(nodesConnected)) == length(nodesConnected)
                                    
                    nodesConnected = nextNode;

                    [~, links_id] = ismember(nodesConnected,Gin.node(currentNode).conn);
                    links_crossing = Gin.link(Gin.node(currentNode).links(links_id));

                    % for each link, compute the individual (xi,yi) points
                    vs = zeros(length(nodesConnected), 2);
                    for j = 1 : length(links_crossing)
                        % if the link connected to the current node has less than 5
                        % points
                        if length(links_crossing(j).point) <= 5
                            [vs(j,1), vs(j,2)]=ind2sub([Gin.w Gin.l],links_crossing(j).point(end));
                        else
                            [vs(j,1), vs(j,2)]=ind2sub([Gin.w Gin.l],links_crossing(j).point(end-5));
                        end
                    end

                    % move all vectors to (x0,y0) and normalize them by their
                    % Euclidean norm
                    v0 = [Gin.node(currentNode).comx Gin.node(currentNode).comy];
                    vs = normr(bsxfun(@minus, vs, v0));
                    % compute the inner product between the first vector and all the other
                    % ones, and evaluate their difference with respect to -1: the one with
                    % the minor difference corresponds to the connected link to v
                    v = repmat(vs(1,:),size(vs,1)-1,1);
                    us = vs(2:end, :);
                    differences = abs(-1 - dot(v,us,2));
                    [~, nodeConn] = min(differences);
                    nodeConn = nodeConn + 1;

                    % rebuilding link 1
                    link1 = swapLink(Gin.link(Gin.node(currentNode).links(find(Gin.node(currentNode).conn==nodesConnected(nodeConn)))), nodesConnected(nodeConn));
                    % rebuilding link 2
                    link2.n2 = currentNode;
                    link2.point = [];
                    
                else
                    
                    link1.n2 = nodesConnected(1);
                    link1.point = [];
                    link2.n2 = nodesConnected(2);
                    link2.point = [];
                    
                end
                
            end
            
        end


    end


end




