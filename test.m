
clear, clc, close all

%% Open the image, the segmentations and the ground truth labels

% Load the image
I = imread('/Users/ignaciorlando/Documents/MATLAB/artery-vein-classification/data/RITE-training/images/21_training.tif');
% Load the segmentation
segm = imread('/Users/ignaciorlando/Documents/MATLAB/artery-vein-classification/data/RITE-training/vessel-segmentations/21_training.png');
% GT labels
load('/Users/ignaciorlando/Documents/MATLAB/artery-vein-classification/data/RITE-training/labels/21_training.png.mat');

%% Generate a skeletonization

% Load the skeletonization
skel = bwmorph(segm, 'skel', Inf);

%% Extract the graph

% Generate graph of crossings
[G_crossings] = initializeGraphFromSkel_new(skel);

% And now transform the graph of crossings into a graph of segments
[Gout] = generateGraphOfSegments(G_crossings);

%% Extract unary and pairwise features

% Extract unary features
[Gout] = compute_unary_features(Gout, I);
% Extract pairwise features
[Gout] = compute_pairwise_features(Gout, segm);

%% Assign labels to each segment

% Assign labels to the graph
[Gout] = assign_labels(Gout, labels);