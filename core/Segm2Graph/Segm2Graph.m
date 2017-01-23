
function [G] = Segm2Graph(segm, thresholds, verbose)

    % default arguments
    if (nargin<2)
        thresholds = [10 10]; % threshold to remove false branches
        verbose = 0; % verbose option
    else
        if length(thresholds)<2
            thresholds = ones(length(thresholds),1) .* thresholds;
        end
    end
    % prevent 'RGB binary masks'
    if (size(segm,3)>1)
        segm = segm(:,:,1);
    end
    
    % remove isolated segments
    

    % skeletonize given segmentation
    if (verbose), disp('Start skeletonization process...'); end
    skel = bwmorph(segm, 'skel', Inf);
    if (verbose), disp('Skeletonization process finished.'); end
     
    % convert skeleton to graph structure
    if (verbose), disp('Start initialization of the graph from skeleton...'); end
    G = initializeGraphFromSkel(skel, segm, thresholds(1));
    if (verbose), disp('Initialization of the graph from skeleton finished.'); end
    
    % clean false positive nodes
    %if (verbose), disp('Start cleaning of false positive nodes...'); end
    %G = cleanNodes(G);
    %if (verbose), disp('Cleaning of false positive nodes finished.'); end

    % solve false branchings due to vessel crossings
    if (verbose), disp('Start solving false branchings...'); end
    [G] = solveFalseBranchings(G, thresholds(2));
    if (verbose), disp('Start solving false branchings...'); end
    
    % Identify vessel crossings 
    if (verbose), disp('Identifying vessel crossings...'); end
    vesselCrossingNodes = extractfield(G.node, 'numLinks');
    if (iscell(vesselCrossingNodes))
        vesselCrossingNodes = cell2mat(vesselCrossingNodes);
    end
    G.vesselCrossings = find(vesselCrossingNodes==4);
    if (verbose), disp('Vessel crossings identified.'); end
    
    % Identify root nodes
    if (verbose), disp('Identifying roots...'); end
    roots = extractfield(G.node, 'isroot');
    if (iscell(roots))
        roots = cell2mat(roots);
    end
    G.roots = find(roots==true);
    if (verbose), disp('Roots identified.'); end
    
    % Compute the chain codes
    if (verbose), disp('Assigning chain codes to all segment links...'); end
    for i = 1 : length(G.link)
        G.link(i).code = getChainCode(G, G.link(i).point);
    end
    if (verbose), disp('Chain codes to all segment links assigned.'); end
    
    % Assign the ONH mask
    G.onh_mask = onh_mask;
    
    % Generate adjacency matrix
    if (verbose), disp('Generating adjacency matrix...'); end
    [G.adjacencyMatrix] = generateAdjacencyMatrix(G);
    if (verbose), disp('Adjacency matrix generated.'); end
    
    
    
    % Identify root nodes
    
%     %solve false branchings due to self cycles
%     if (verbose), disp('Start solving self cycles...'); end
%     [G] = solveSelfCycles(G);
%     if (verbose), disp('Start solving self cycles...'); end
%    
%     %solve iddle segments and nodes
%     if (verbose), disp('Start solving iddle segments and nodes...'); end
%     [G] = solveIddleLinksAndNodes(G);
%     if (verbose), disp('Start solving iddle segments and nodes...'); end
    
end


function [A] = generateAdjacencyMatrix(G)

    % create adjacency matrix
    A = zeros(length(G.node));
    % for each node
    for i=1:length(G.node)
        % get nodes connected to the current node 
        idx1=find(G.node(i).conn>0);
        % get all links to all nodes
        idx2=find(G.node(i).links>0);
        % get the intersection
        idx=intersect(idx1,idx2);
        % for each edge, save the connection and the length of the link
        for j=1:length(idx)
            if(i==G.link(G.node(i).links(idx(j))).n1)
                A(i,G.link(G.node(i).links(idx(j))).n2)=length(G.link(G.node(i).links(idx(j))).point);
                A(G.link(G.node(i).links(idx(j))).n2,i)=length(G.link(G.node(i).links(idx(j))).point);
            end;
            if(i==G.link(G.node(i).links(idx(j))).n2)
                A(i,G.link(G.node(i).links(idx(j))).n1)=length(G.link(G.node(i).links(idx(j))).point);
                A(G.link(G.node(i).links(idx(j))).n1,i)=length(G.link(G.node(i).links(idx(j))).point);
            end;
        end;
    end;

end
