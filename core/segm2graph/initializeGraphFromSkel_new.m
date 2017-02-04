
function [Gout] = initializeGraphFromSkel_new(skel)

    % Compute the intersecting pts
    intersecting_pts = find_skel_intersection_2(skel);

    % Get indices
    %idx = sub2ind(size(skel), intersecting_pts(:,2), intersecting_pts(:,1));

    % Generate only branching points
    branching_points = false(size(skel));
    branching_points(intersecting_pts) = true;

    % Generate identified junctions
    identified_junctions = skel .* branching_points;

    % Remove junctions
    skel_without_junctions = skel .* imcomplement(branching_points);

    
    
    % get image dimensions
    w=size(skel,1);
    l=size(skel,2);

    % change image boundaries to zero (in the 3D version, the entire image is
    % set to 0 because the Z=1 and Z=end positions of the matrix are set to 0)
    skel(1,:)=0;
    skel(:,1)=0;
    skel(end,:)=0;
    skel(:,end)=0;

    
    
    % this variable will be used to label each node afterwards
    pixel_labels = double(skel);

    % get all foreground pixels
    foreground_pixels = find(skel);
    
    % get 8-neighborhoods of all canal pixels
    nh = logical(pk_get_nh(skel,foreground_pixels));
    % get 8-neighborhoods indices of all canal pixels
    nhi = pk_get_nh_idx(skel,foreground_pixels);

    % sum the number of foreground pixel neighbors in the 8-neighborhood
    sum_nh = sum(logical(nh),2);
    
    %intersecting_pts = find_skel_intersection(skel);
    %nodes = sub2ind(size(skel), intersecting_pts(:,2), intersecting_pts(:,1));
    nodes = find(identified_junctions);

    % all pixels with exactly 3 nb are part of what we call a canal
    % 0 0 0 0 0
    % - - X - -
    % 0 0 0 0 0
    canal_pixels = find(skel_without_junctions);
    
    % can_nh contains the neighbors of pixels in the canal
    % can_nh_idx contains indexes of neighbors of pixels in the canal
    can_nh_idx = pk_get_nh_idx(skel,canal_pixels);
    can_nh = pk_get_nh(skel,canal_pixels);
    % the central pixel of the 3x3 nbh is removed
    can_nh_idx(:,5)=[];
    can_nh(:,5)=[];
    % keep only the two existing foreground pixels
    can_nb = sort(logical(can_nh).*can_nh_idx,2);
    % remove zeros
    can_nb(:,1:end-2) = [];
    % add neighbours to canalicular pixel list (this might include nodes)
    canal_pixels = [canal_pixels can_nb];

    % group clusters of node pixels to nodes
    node=[];
    link=[];

    tmp=false(w,l);
    tmp(nodes)=1; % true node pixels 
    cc2=bwconncomp(tmp); % number of unique nodes

    % create node structure
    % for each node structure
    for i=1:cc2.NumObjects
        % get pixels involved
        node(i).idx = cc2.PixelIdxList{i};
        node(i).links = [];
        node(i).conn = [];
        node(i).numLinks = 0;
        [x y]=ind2sub([w l],node(i).idx);
        % the node is assumed to be located at the average point of all pixels
        % that are part of the node structure
        node(i).comx = mean(x);
        node(i).comy = mean(y);
        % assign index to node pixels
        pixel_labels(node(i).idx) = i+1;
    end;
    
    
    % link iterator
    last_link_idx = 1;
    last_node_idx = length(node);
    % for each node
    for i=1:length(node)

        % find all the neighbors indices of all the pixels in the node
        link_idx = find(ismember(nhi(:,5),node(i).idx));

        % for each canal pixel element on the nbh of the node pixels
        for j=1:length(link_idx)
            %j
            % visit all pixel of this node

            % all potential unvisited links emanating from this pixel:
            % get all neighbors of the pixel link_idx(j) such that they are part
            % of the skeleton
            link_cands = nhi(link_idx(j),nh(link_idx(j),:)==1); % get link candidates
            % but get only those that are not nodes (pixel_label has 1 in
            % all node pixels)
            link_cands = link_cands(pixel_labels(link_cands)==1); % remove pixels that are not part of the canal

            % for each pixel candidate
            for k=1:length(link_cands)
                % follow the link and obtain
                % - edge: list of pixels of the edge
                % - end_node_idx: id of the last node (-1 if it terminates without reaching a node)
                [edge, end_node_idx] = pk_follow_link(pixel_labels, node, i, j, link_cands(k), canal_pixels);
                   
                % remove from skel2 the entire edge (marking as visited)
                pixel_labels(edge(2:end-1))=0;
                if(end_node_idx<0) % for endpoints, also remove last pixel
                    pixel_labels(edge(end))=0;
                end; 
                % only large branches or non-loops
                if((end_node_idx<0 && length(edge)>3) || (i~=end_node_idx && end_node_idx>0))
                    
                    % encode link information
                    link(last_link_idx).n1 = i;
                    link(last_link_idx).n2 = end_node_idx; % node number
                    link(last_link_idx).point = edge;
                    node(i).links = [node(i).links, last_link_idx];
                    node(i).numLinks = length(node(i).links);
                    node(i).conn = [node(i).conn, end_node_idx];
                    % if it ends with a node, we encode also the other
                    % direction of the link
                    if(end_node_idx>0)
                        node(end_node_idx).links = [node(end_node_idx).links, last_link_idx];
                        node(end_node_idx).numLinks = length(node(end_node_idx).links);
                        node(end_node_idx).conn = [node(end_node_idx).conn, i];
                    end;
                    % update the node iterator
                    last_link_idx = last_link_idx + 1;
                end;

            end;
        end;

    end;

    % assign information to the graph structure
    Gout.node = node;
    Gout.link = link;
    Gout.w = w;
    Gout.l = l;

    
end




function nhoods = pk_get_nh(skelImage,foregroundPixels)
% pk_get_nh return the 8-pixels neighborhood for each foreground pixel in
% foregroundPixels according to the skelImage skeleton
%
% Input
% - skelImage = skeletonization of the given segmentation
% - foregroundPixels = set of foreground pixels in skelImage
%
% Output
% - nhoods = set of neighbors for each given pixel
%

    width = size(skelImage,1);
    height = size(skelImage,2);

    [x y]=ind2sub([width height],foregroundPixels);

    nhoods = false(length(foregroundPixels),8);

    for xx=1:3
        for yy=1:3
            w=sub2ind([3 3],xx,yy);
            idx = sub2ind([width height],x+xx-2,y+yy-2);
            nhoods(:,w)=skelImage(idx);
        end;
    end;
end


function nhoods = pk_get_nh_idx(skel,foregroundPixels)
% pk_get_nh_idx return the indexes of each pixel in the 
% 8-pixels neighborhood for each foreground pixel in
% foregroundPixels according to the skelImage skeleton
%
% Input
% - skelImage = skeletonization of the given segmentation
% - foregroundPixels = set of foreground pixels in skelImage
%
% Output
% - nhoods = set of neighbors for each given pixel

    width = size(skel,1);
    height = size(skel,2);

    [x y]=ind2sub([width height],foregroundPixels);

    nhoods = zeros(length(foregroundPixels),9);

    for xx=1:3
        for yy=1:3
            w=sub2ind([3 3],xx,yy);
            nhoods(:,w) = sub2ind([width height],x+xx-2,y+yy-2);
        end;
    end;
end


function [edge, end_node_idx] = pk_follow_link(pixel_labels, nodesStructure, sourceNode, firstPixel, canalCand, cans)
% pk_follow_link

    edge = [];
    end_node_idx = [];

    % assign start node to first pixel
    edge(1) = nodesStructure(sourceNode).idx(firstPixel);

    i=1;
    isdone = false;
    % while no node reached
    while(~isdone) 
        i=i+1; % next pixel
        current = find(cans(:,1)==canalCand,1);
        % if it is not the last pixel of the canal
        if(~isempty(current))

            nextCand = cans(current,2);
            %if it goes back or any of the next candidate neighbors already
            %is in the edge
            if nextCand==edge(i-1)      % <----------- EL PROBLEMA ANDA POR AC?
               % switch direction to avoid loops
               nextCand = cans(current,3);
            end;

            % if it is a node
            if (nextCand<1)
                end_node_idx = -1;
                isdone = 1;
            else
                if(pixel_labels(nextCand)>1)
                    % update the edge information and return
                    edge(i) = canalCand;
                    edge(i+1) = nextCand; % first node
                    end_node_idx = pixel_labels(nextCand)-1; % node #
                    isdone = 1;
                else
%                     if ~isempty(find(nextCand==edge)) % <----------- EL PROBLEMA ANDA POR AC?
%                         edge(i) = canalCand;
%                         end_node_idx = -1;
%                         isdone = 1;
%                     else
                        % then is not a node and is not in one of our cicles, then continue 
                        edge(i) = canalCand;
                        canalCand = nextCand;
%                     end
                end;
            end
        else
            edge(i) = canalCand;
            end_node_idx = -1;
            isdone = 1;
        end;
    end;
end