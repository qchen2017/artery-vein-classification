
% SCRIPT_SETUP_RITE_DATA
% -------------------------------------------------------------------------
% This script download RITE (DRIVE) data set and prepares all the data to 
% be used for experiments.
% -------------------------------------------------------------------------

%% set up variables

% set up main variables
config_setup_rite_data;

% create folder
mkdir(fullfile(output_folder, 'tmp'));

%% prepare the input data set

% prepare ZIP filename
zip_filename = fullfile(input_folder, 'AV_groundTruth.zip');

% check if the file exists
if exist(zip_filename, 'file') ~= 0
    % unzip the file
    fprintf('Unzipping AV_groundTruth file...\n');
    % unzip on output_folder/tmp
    unzip(zip_filename, fullfile(output_folder, 'tmp'));
else
    throw('File not found');
end

%% generate the input data set

subsets = {'training', 'test'};

for subs_ = 1 : length(subsets)
    
    % get current set
    current_set = subsets{subs_};

    % the input folder will be different now
    input_folder_for_images = fullfile(output_folder, 'tmp', 'AV_groundTruth', current_set, 'images');
    input_folder_for_labels = fullfile(output_folder, 'tmp', 'AV_groundTruth', current_set, 'av');

    % and the output folder will be...
    rite_dataset_folder = fullfile(output_folder, strcat('RITE-', current_set));
    mkdir(rite_dataset_folder);
    mkdir(fullfile(rite_dataset_folder, 'images'));
    mkdir(fullfile(rite_dataset_folder, 'labels'));
    mkdir(fullfile(rite_dataset_folder, 'vessel-segmentations'));
    mkdir(fullfile(rite_dataset_folder, 'veins'));
    mkdir(fullfile(rite_dataset_folder, 'arteries'));

    % retrieve image names
    image_filenames = dir(fullfile(input_folder_for_images, '*.tif'));
    image_filenames = {image_filenames.name};
    % retrieve labels names
    labels_filenames = dir(fullfile(input_folder_for_labels, '*.png'));
    labels_filenames = {labels_filenames.name};

    % for each graph file
    for i = 1 : length(image_filenames)

        % copy current image
        copyfile(fullfile(input_folder_for_images, image_filenames{i}), fullfile(rite_dataset_folder, 'images', image_filenames{i}));
        % copy all the annotations
        copyfile(fullfile(input_folder_for_labels, labels_filenames{i}), fullfile(rite_dataset_folder, 'labels', labels_filenames{i}));
        % open label
        labels = imread(fullfile(input_folder_for_labels, labels_filenames{i}));
        % identify veins and save
        imwrite(labels(:,:,3)==255 .* labels(:,:,1)==0, fullfile(rite_dataset_folder, 'veins', labels_filenames{i}));
        % identify arteries and save
        imwrite(labels(:,:,1)==255 .* labels(:,:,3)==0, fullfile(rite_dataset_folder, 'arteries', labels_filenames{i}));
        % identify vessels and save
        imwrite(sum(labels,3) > 0, fullfile(rite_dataset_folder, 'vessel-segmentations', labels_filenames{i}));

    end

end

%% delete tmp folder
rmdir(fullfile(output_folder, 'tmp'), 's');