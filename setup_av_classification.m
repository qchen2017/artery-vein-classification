
% SETUP_AV_CLASSIFICATION
% -------------------------------------------------------------------------
% This code add folders to Matlab environment
% -------------------------------------------------------------------------

% get current root position
my_root_position = pwd;

% add configuration folders to path
addpath(genpath(fullfile(my_root_position, 'default-configuration')));

% add main folders to path
addpath(genpath(fullfile(my_root_position, 'data-organization'))) ;

clear
clc