
% SCRIPT_COMPUTE_TIME_STATISTICS_FOR_DATA_PREPARATION
% -------------------------------------------------------------------------
% This script will generate unary and pairwise features to evaluate time
% statistics.
% -------------------------------------------------------------------------

config_compute_time_statistics_for_data_preparation;

%% Setup folders

% Complete data folder
data_folder = fullfile(data_folder, dataset_name);
% Complete results folder
results_folder = fullfile(results_folder, dataset_name);

%% Evaluate time in data
   
fprintf('Analyzing times on data...\n');

% Initialize the folder
mkdir(fullfile(data_folder, 'training_set'));
% Extract training data
[data, computational_times] = prepare_graphs(data_folder, true);


  