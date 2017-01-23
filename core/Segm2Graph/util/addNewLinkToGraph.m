
function [Gout] = addNewLinkToGraph(Gout, link_idx, new_link)

    % add node information
    Gout.node(new_link.n1).links = [Gout.node(new_link.n1).links, link_idx]; % update links that begin/end in the node
    Gout.node(new_link.n1).numLinks = length(Gout.node(new_link.n1).links); % update the number of links connected to the node
    Gout.node(new_link.n1).conn = [Gout.node(new_link.n1).conn, new_link.n2]; % update the connections
    if (new_link.n2~=-1) % if the sink node is an active node...
        Gout.node(new_link.n2).links = [Gout.node(new_link.n2).links, link_idx]; % update the links that begin/end in the node
        Gout.node(new_link.n2).numLinks = length(Gout.node(new_link.n2).links); % update the number of links connected to the node
        Gout.node(new_link.n2).conn = [Gout.node(new_link.n2).conn, new_link.n1]; % update the connections
    end

    % assign new link to the graph structure
    Gout.link(link_idx).n1 = new_link.n1;
    Gout.link(link_idx).n2 = new_link.n2;
    Gout.link(link_idx).point = new_link.point;
        
end