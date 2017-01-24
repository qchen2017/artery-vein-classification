
clear, clc, close all

%% Open the image, the segmentations and the ground truth labels
fprintf('Loading an image, its vessel segmentation and its ground truth labelling\n');

% Load the image
I = imread('/Users/ignaciorlando/Documents/MATLAB/artery-vein-classification/data/RITE-training/images/21_training.tif');
% Load the segmentation
segm = imread('/Users/ignaciorlando/Documents/MATLAB/artery-vein-classification/data/RITE-training/vessel-segmentations/21_training.png');
% GT labels
load('/Users/ignaciorlando/Documents/MATLAB/artery-vein-classification/data/RITE-training/labels/21_training.png.mat');

%% Generate a skeletonization
fprintf('Skeletonizating...\n');

% Load the skeletonization
skel = bwmorph(segm, 'skel', Inf);

%% Extract the graph

% Generate graph of crossings
fprintf('Extracting graph of vessel crossings and branching points...\n');
[G_crossings] = initializeGraphFromSkel_new(skel);

% And now transform the graph of crossings into a graph of segments
fprintf('Extracting graph of segments...\n');
[Gout] = generateGraphOfSegments(G_crossings);

%% Extract unary and pairwise features

% Extract unary features
fprintf('Extracting unary features...\n');
[Gout] = compute_unary_features(Gout, I);
fprintf('Extracting pairwise features...\n');
% Extract pairwise features
[Gout] = compute_pairwise_features(Gout, segm);

%% Assign labels to each segment
fprintf('Assigning ground truth labels to each segment...\n');

% Assign labels to the graph
[Gout] = assign_labels(Gout, labels);

%% Output a labeled graph
fprintf('Printing the graph but with ground truth labels...\n');

% Retrieve labeled graph
[labeled_map] = generate_image_from_classified_graph(Gout);
figure, imshow(labeled_map);

%% Let's give a try on using the CRF unsupervisedly
fprintf('Applying a CRF...\n');

% Get energies
[A,T] = get_crf_energies_from_graph(Gout, ones(12,2), ones(1,2));
% Solve using maxflow
[flow, labels] = maxflow(A,T);
