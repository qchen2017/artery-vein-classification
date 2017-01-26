
% SCRIPT_TRAIN_ARTERY_VEIN_CLASSIFIER
% -------------------------------------------------------------------------
% This script will train a CRF using a SOSVM for artery/vein
% classification.
% -------------------------------------------------------------------------

config_train_artery_vein_classifier;

%% Setup folders

% Complete data folder
data_folder = fullfile(data_folder, dataset_name);
% Complete results folder
results_folder = fullfile(results_folder, dataset_name);

%% Prepare the training data

% initialize the training data filename
training_data_filename = fullfile(data_folder, 'training_set','training_data.mat'); 

% if the file exists, load it
if exist(training_data_filename, 'file')==0
    
    fprintf('Computing training and validation data...\n');
    
    % Initialize the folder
    mkdir(fullfile(data_folder, 'training_set'));
    % Extract training data
    [data] = prepare_graphs(data_folder, true);
    % Use 70% of the data as training set and the remaining 30% for
    % validation
    N_training = floor(length(data) * 0.7);
    training_data = data(1:N_training);
    validation_data = data(N_training + 1:end);
    % Save the file
    save(training_data_filename, 'training_data', 'validation_data');
    
else
    
    fprintf('Loading training and validation data...\n');
    
    % Load the training data file
    load(training_data_filename);
end

%% Train a CRF using a SOSVM

% learn an artery/vein classifier
[model, performance_on_validation] = learn_artery_vein_classifier(training_data, validation_data, c_values);
