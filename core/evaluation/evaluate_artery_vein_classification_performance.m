
function performance = evaluate_artery_vein_classification_performance(classified_graph, ground_truth_graph, evaluation_metric)

    % initialize performance in 0
    performance = 0;

    % depending on the evaluation metric, compute performance
    switch evaluation_metric
        
        % Acc = (TP + TN) / (TP + TN + FP + FN)
        case 'accuracy'
            
            % for each node
            for i = 1 : length(classified_graph.node)
                % sum TP + TN when the label is accurate
                if (ground_truth_graph.node(i).label~=0) && (classified_graph.node(i).label == ground_truth_graph.node(i).label)
                    performance = performance + 1;
                end
            end
            % divide performance by the number of nodes
            %performance = performance;
            
    end
    
end