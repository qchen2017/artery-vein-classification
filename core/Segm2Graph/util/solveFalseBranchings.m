
function [Gout] = solveFalseBranchings(Gin, threshold)

    % Initialize the new graph
    Gout.node = [];
    Gout.link = [];
    Gout.w = Gin.w;
    Gout.l = Gin.l;
    
    % nodesAliases is used to identify the equivalent index of each new
    % node
    % ----> nodes=-1 are nodes with false branchings
    % ----> nodes=0 were not named yet
    node_idxs_three_neighbors = find((extractfield(Gin.node, 'numLinks')==3) .* not(cell2mat(extractfield(Gin.node, 'isroot'))));
    toRemove = [];
    for i = 1 : length(node_idxs_three_neighbors)
        % remove nodes with 2 branches to the same node (nodes with 2
        % terminal links)
        if length(unique(Gin.node(node_idxs_three_neighbors(i)).conn)) ~= length(Gin.node(node_idxs_three_neighbors(i)).conn)
            toRemove = [toRemove node_idxs_three_neighbors(i)];
        end
    end
    node_idxs_three_neighbors(ismember(node_idxs_three_neighbors, toRemove)) = [];
    nodesAliases = zeros(size(Gin.node));
    iteratorThreeNeighbors = 1;
    iteratorFalseBranchings = 1;
    falseBranchings ={};
    visited_nodes = false(length(node_idxs_three_neighbors), 1);
    
    while (iteratorThreeNeighbors < length(node_idxs_three_neighbors))
        
        currentNode = node_idxs_three_neighbors(iteratorThreeNeighbors);

        % si el nodo est? conectado a otro que est? en la lista de nodos de 3 arcos
        if (~visited_nodes(iteratorThreeNeighbors)) && (sum(ismember(Gin.node(currentNode).conn, node_idxs_three_neighbors))>=1)
            
            visited_nodes(iteratorThreeNeighbors) = true;
            % obtengo los links conectados a nodos de 3 arcos
            links_idxs = Gin.link(Gin.node(currentNode).links);
            links_idxs = links_idxs(ismember(Gin.node(currentNode).conn, node_idxs_three_neighbors));
            % me fijo si tiene alg?n link corto
            for j = 1 : length(links_idxs)
                % si la longitud del link es muy corta
                if (length(links_idxs(j).point) <= threshold)
                    % la lista de false branchings
                    if links_idxs(j).n1 ~= currentNode
                        otherNodeToRemove = links_idxs(j).n1;
                    else
                        otherNodeToRemove = links_idxs(j).n2;
                    end
                    
                    
                    % identifico n1 y n2
                    n1 = Gin.node(currentNode);
                    n2 = Gin.node(otherNodeToRemove);

                    % controlo si los links que en el futuro tendr?a que
                    % unir porque son crossings en realidad van ambos a -1

                    % identify links crossing
                    links_crossing = [n1.links(n1.conn~=links_idxs(j).n2), n2.links(n2.conn~=links_idxs(j).n1)];
                    nodesConnected = [Gin.node(currentNode).conn, Gin.node(otherNodeToRemove).conn];
                    nodesConnected(nodesConnected==currentNode) = [];
                    nodesConnected(nodesConnected==otherNodeToRemove) = [];                   
                    % for each link, compute the individual (xi,yi) points
                    vs = zeros(length(nodesConnected), 2);
                    for l_c = 1 : length(links_crossing)
                        % if the link connected to the current node has less than 5
                        % points
                        if length(Gin.link(links_crossing(l_c)).point) <= 5
                            [vs(l_c,1), vs(l_c,2)]=ind2sub([Gin.w Gin.l],Gin.link(links_crossing(l_c)).point(end));
                        else
                            [vs(l_c,1), vs(l_c,2)]=ind2sub([Gin.w Gin.l],Gin.link(links_crossing(l_c)).point(end-5));
                        end
                    end
                    % move all vectors to the point at which will be located the new node,
                    % and normalize them by their Euclidean norm
                    v0 = [mean([n1.comx n2.comx]) mean([n1.comy n2.comy])];
                    vs = normr(bsxfun(@minus, vs, v0));
                    % compute the inner product between the first vector and all the other
                    % ones, and evaluate their difference with respect to -1: the one with
                    % the minor difference corresponds to the connected link to v
                    v = repmat(vs(1,:),size(vs,1)-1,1);
                    us = vs(2:end, :);
                    differences = abs(-1 - dot(v,us,2));
                    [~, link_connected_to_vs_1] = min(differences);
                    link_connected_to_vs_1 = link_connected_to_vs_1 + 1;
                    % now, check if both vs_1 and vs_link_connected_to_vs_1
                    % are connected to -1. In that case, the nodes do not
                    % need to be transformed to a vessel crossing
                    the_other_joined_links_are = ones(4,1);
                    the_other_joined_links_are(1) = 0; the_other_joined_links_are(link_connected_to_vs_1) = 0;
                    the_other_joined_links_are = find(the_other_joined_links_are);
                    if ~(((Gin.link(links_crossing(1)).n2==-1) && (Gin.link(links_crossing(link_connected_to_vs_1)).n2==-1)) ...
                        || ((Gin.link(links_crossing(the_other_joined_links_are(1))).n2==-1) && (Gin.link(links_crossing(the_other_joined_links_are(2))).n2==-1)))
                        
                        % lo marco para eliminarlo
                        nodesAliases(currentNode) = -1;
                        % elimino el nodo conectado al actual por este 
                        % arco, porque tambi?n es un nodo con false
                        % branching
                        node_idxs_three_neighbors(node_idxs_three_neighbors==otherNodeToRemove) = [];
                        nodesAliases(otherNodeToRemove)=-1;
                        % agrego el arco a la lista de false branchings
                        falseBranchings{iteratorFalseBranchings} = links_idxs(j);
                        iteratorFalseBranchings = iteratorFalseBranchings + 1;
                        
                    else
                        
                        % lo marco para no eliminarlo
                        nodesAliases(currentNode) = 0;
                        % y marco tambi?n a su vecino para no borrarlo
                        nodesAliases(otherNodeToRemove)=0;
                        
                    end
                    
                    [~, visited_idx] = ismember(otherNodeToRemove,node_idxs_three_neighbors);
                    if (visited_idx~=0)
                        visited_nodes(visited_idx) = 1;
                    end
                    
                end
            end
            
        end       
        iteratorThreeNeighbors = iteratorThreeNeighbors + 1;
        
    end
    
    % obtener el grafo considerando nodos y links que no tengan false
    % branchings
    [Gout, nodesAliases] = addNodesWithoutConflicts(Gin, Gout, nodesAliases);
    
    % agrego los nodos en la ubicaci?n promedio, y conservo sus alias para
    % cuando reconstruya los arcos
    newNodesAliases = zeros(size(falseBranchings));
    falseBranchingNodes = [];
    for i = 1 : length(falseBranchings)
        
        % calculo el nuevo alias
        newNodesAliases(i) = length(Gout.node) + 1;

        % obtengo los dos nodos del false branching
        n1 = Gin.node(falseBranchings{i}.n1);
        n2 = Gin.node(falseBranchings{i}.n2);
        
        % agrego el nodo
        comx = mean([n1.comx n2.comx]);
        comy = mean([n1.comy n2.comy]);
        posx = round(comx);
        posy = round(comy);
        Gout.node(newNodesAliases(i)).idx = sub2ind([Gin.w Gin.l], posx, posy);
        Gout.node(newNodesAliases(i)).conn = [];
        Gout.node(newNodesAliases(i)).links = [];
        Gout.node(newNodesAliases(i)).comx = posx;
        Gout.node(newNodesAliases(i)).comy = posy;
        Gout.node(newNodesAliases(i)).numLinks = 0;
        Gout.node(newNodesAliases(i)).isroot = false;
        
    end
    
    % ahora resuelvo los arcos
    for i = 1 : length(falseBranchings)
              
        % identifico n1 y n2
        n1 = Gin.node(falseBranchings{i}.n1);
        n2 = Gin.node(falseBranchings{i}.n2);
 
        % obtengo los links de entrada a n1 y n2
        [l11, l12] = getInputLinks(n1, falseBranchings{i}.n2, Gin);
        [l21, l22] = getInputLinks(n2, falseBranchings{i}.n1, Gin);
         
        % genero los nuevos links
        [new_l11] = getNewLink(Gin, l11, falseBranchings{i}.n1, Gout.node(newNodesAliases(i)), newNodesAliases(i), nodesAliases, falseBranchings, newNodesAliases, Gout);
        node_connected_11 = [new_l11.n1 new_l11.n2];
        node_connected_11 = node_connected_11(node_connected_11~=newNodesAliases(i));
        
        [new_l12] = getNewLink(Gin, l12, falseBranchings{i}.n1, Gout.node(newNodesAliases(i)), newNodesAliases(i), nodesAliases, falseBranchings, newNodesAliases, Gout);
        node_connected_12 = [new_l12.n1 new_l12.n2];
        node_connected_12 = node_connected_12(node_connected_12~=newNodesAliases(i));
        [new_l21] = getNewLink(Gin, l21, falseBranchings{i}.n2, Gout.node(newNodesAliases(i)), newNodesAliases(i), nodesAliases, falseBranchings, newNodesAliases, Gout);
        node_connected_21 = [new_l21.n1 new_l21.n2];
        node_connected_21 = node_connected_21(node_connected_21~=newNodesAliases(i));
        [new_l22] = getNewLink(Gin, l22, falseBranchings{i}.n2, Gout.node(newNodesAliases(i)), newNodesAliases(i), nodesAliases, falseBranchings, newNodesAliases, Gout);
        node_connected_22 = [new_l22.n1 new_l22.n2];
        node_connected_22 = node_connected_22(node_connected_22~=newNodesAliases(i));

        nodes_already_connected = Gout.node(newNodesAliases(i)).conn;
        nodes_already_connected = nodes_already_connected(nodes_already_connected~=-1);
        
        % agrego cada link
        if ~ismember(node_connected_11, nodes_already_connected)
            [Gout] = addNewLinkToGraph(Gout, length(Gout.link)+1, new_l11);
        end
        if ~ismember(node_connected_12, nodes_already_connected)
            [Gout] = addNewLinkToGraph(Gout, length(Gout.link)+1, new_l12);
        end
        if ~ismember(node_connected_21, nodes_already_connected)
            [Gout] = addNewLinkToGraph(Gout, length(Gout.link)+1, new_l21);
        end
        if ~ismember(node_connected_22, nodes_already_connected)
            [Gout] = addNewLinkToGraph(Gout, length(Gout.link)+1, new_l22);
        end
        
    end

end



function [newlink] = getNewLink(Gin, link, ni_idx, nn, nn_idx, nodesAliases, falseBranchings, newNodesAliases, Gout)

    % obtengo el punto (x1,y1)
    if (link.n2 ~= ni_idx)
        link.point = fliplr(link.point);
        link.n1 = link.n2;
        link.n2 = ni_idx;
    end
    if (length(link.point) > 5)
        [x1, y1] = ind2sub([Gin.w Gin.l], link.point(end - 5));
        startingpoints = link.point(1 : end-5);
    else
        [x1, y1] = ind2sub([Gin.w Gin.l], link.point(1));
        startingpoints = link.point(1);
    end
    
    % obtengo el (x2,y2), que es el centro del nuevo nodo
    x2 = nn.comx;
    y2 = nn.comy;
    
    % actualizo los ?ltimos 5 puntos del link anterior
    [x, y] = bresenham(x1,y1,x2,y2);
    endpoints = sub2ind([Gin.w Gin.l], x, y);
    link.point = [startingpoints, endpoints'];
    
    % create the new link
    if (link.n1 ~= -1)
        newlink.n1 = nodesAliases(link.n1);
        newlink.n2 = nn_idx;
    else
        newlink.n1 = nn_idx;
        newlink.n2 = -1;
        link.point = fliplr(link.point);
    end
    newlink.point = link.point;
    
    % si el nodo est? conectado a un nodo que tambi?n tiene un false
    % branching, pongo como nodo n1 al alias correspondiente al nodo del
    % false branching
    if (newlink.n1==-1)
        
        % detecto cu?l es el false branching donde est? involucrado n1
        i = 1;
        while ((i<length(falseBranchings)) && ~((falseBranchings{i}.n1 == link.n1) || (falseBranchings{i}.n2 == link.n1)))
            i = i + 1;
        end
        
        % asigno el alias de acuerdo al alias nuevo
        newlink.n1 = newNodesAliases(i);
        
        % corrijo el otro extremo
        if (length(newlink.point) > 5)
            [x2, y2] = ind2sub([Gin.w Gin.l], newlink.point(5));
            endpoints = newlink.point(5 : end);
        else
            [x2, y2] = ind2sub([Gin.w Gin.l], newlink.point(end));
            endpoints = newlink.point(end);
        end
        
        x1 = Gout.node(newlink.n1).comx;
        y1 = Gout.node(newlink.n1).comy;
        
        % actualizo los ?ltimos 5 puntos del link anterior
        [x, y] = bresenham(x1,y1,x2,y2);
        startingpoints = sub2ind([Gin.w Gin.l], x, y);
        newlink.point = [startingpoints', endpoints];
        
    end
end



function [li1, li2] = getInputLinks(ni, nj, Gin)

    inlinks = ni.links(ni.conn~=nj);
    li1 = Gin.link(inlinks(1));
    li2 = Gin.link(inlinks(2));

end

