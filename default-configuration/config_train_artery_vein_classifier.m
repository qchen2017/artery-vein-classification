
% CONFIG_TRAIN_ARTERY_VEIN_CLASSIFIER
% -------------------------------------------------------------------------
% This script is called by script_train_artery_vein_classifier for setting
% up parameters before training an artery/vein classifier.
% -------------------------------------------------------------------------

% Folder where the data is
data_folder = './data';

% Data set name
dataset_name = 'RITE-training';

% Results folder
results_folder = './results';

% C values
c_values = 10.^[-8:1:8];
%c_values = 10;