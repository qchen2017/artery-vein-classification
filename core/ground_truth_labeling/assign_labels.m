
function [Gout] = assign_labels(Gout, labels)

    % assign labels to each node
    for i = 1 : length(Gout.node)
        
        % retrieve label on the segment pixels
        labels_on_segment = unique(labels(Gout.node(i).idx));
        labels_on_segment(labels_on_segment==0) = [];
        
        % assign the most frequent label in the segment as the ground truth
        % label of the node
        if length(labels_on_segment) > 1
            labels_on_segment = mode(labels_on_segment);
        end
        if labels_on_segment==2
            labels_on_segment = 0;
        end
        Gout.node(i).label = labels_on_segment;
        
    end
    
    % Assign labelling as an image mask
    Gout.labels = labels;

end