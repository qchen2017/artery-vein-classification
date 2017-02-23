

function [model, performance_on_validation] = learn_artery_vein_classifier(training_data, validation_data, C_space, evaluation_metric)

    % by default, our evaluation metric is accuracy
    if nargin < 4
        evaluation_metric = 'accuracy';
    end
    
    % number of segments on the training set
    n = 0;
    for i = 1 : length(training_data)
        n = n + length(find(training_data{i}.all_my_labels));
    end
    % and the validation set
    n_val = 0;
    for i = 1 : length(validation_data)
        n_val = n_val + length(find(validation_data{i}.all_my_labels));
    end
    
    % initialize an array of accuracies on the validation set
    accuracies_on_validation = zeros(length(C_space), 1);
    % and a cell array of models
    learned_models = cell(size(accuracies_on_validation));
    
    % for each value of C in the space of C's we want to explote
    for c_idx = 1 : length(C_space)
        
        fprintf('\nTraining with %d\n', C_space(c_idx));
        
        % learn a model
        [learned_models{c_idx}] = learn_artery_vein_classifier_for_a_given_c(training_data, C_space(c_idx), n);
        % apply model on validation data
        results = cell(size(validation_data));
        
        correctly_classified_samples = 0;
        
        for j = 1 : length(validation_data)
            
            % Classify arteries and veins
            results{j} = classify_arteries_and_veins(validation_data{j}, learned_models{c_idx});

            correctly_classified_samples = correctly_classified_samples + evaluate_artery_vein_classification_performance(results{j}, validation_data{j}, evaluation_metric);
            
%             if (j==1)
%                 figure
%                 
%                 subplot(1, 2, 1);
%                 [II] = generate_image_from_classified_graph(results{j});
%                 imshow(II)
%                 title('our algorithm')
%                 
%                 subplot(1, 2, 2);
%                 [II] = generate_image_from_classified_graph(validation_data{j});
%                 imshow(II)
%                 title('ground truth')
%             end
            
        end
        % evaluate this model on the validation set
        accuracies_on_validation(c_idx) = correctly_classified_samples / n_val; 
        
        %accuracies_on_validation(c_idx) = ...
        %    evaluate_artery_vein_classification_performance(all_classified_labels, all_ground_truth_labels, evaluation_metric);
        
        
        fprintf('\nAccuracy on validation: %d\n\n', accuracies_on_validation(c_idx));
        
    end
    
    % Get the maximum performance on the validation set
    [performance_on_validation, idx] = max(accuracies_on_validation);
    
    % plot the evolution of C values
    figure, plot(log10(C_space), accuracies_on_validation, 'LineWidth', 2);
    ylim([0 1]);
    xlabel('$i = \log_{10}(C)$ values', 'Interpreter', 'LaTex');
    ylabel('Accuracy', 'Interpreter', 'LaTex');
    box on;
    grid on;
    title('Accuracy evolution on the validation set per each value of $C$', 'Interpreter', 'LaTex');
    
    % And now return the best model
    model = learned_models{idx};

end


function [model] = learn_artery_vein_classifier_for_a_given_c(graphs, C, n)

    % Set up the configuration structure
    % -----------------------------------------------
    
    % Regularization parameter
    config.C = C;% / n;
    % Size of the weight vector:
    %
    % || psi || = || w_u || + || w_p ||
    %
    % where || w_u || is 2 x the number of unary features + 2 (for the bias
    %                 terms)
    %       || w_p || is the number of pairwise features
    
    config.size_w_u = 2 * (graphs{1}.properties.unary_dim + 1);
    config.size_w_p = graphs{1}.properties.pairwise_dim;
    
    %config.size_w_u = 2 * (1 + 1);
    %config.size_w_p = graphs{1}.properties.pairwise_dim * 4;
    
    config.sizePsi = config.size_w_u + config.size_w_p;
    % No positivity constraints needed
    %config.posindx = [];
    config.posindx = config.size_w_u + 1:1:config.sizePsi;
    % Add n to the config struct
    config.n = n;
    
    % -----------------------------------------------
    
    % Learn the CRF using SOSVM
    [model, config] = sosvmCallback(config, graphs);
    
    % Assign the configuration to the model
    model.config = config;

end