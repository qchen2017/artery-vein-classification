
clear, clc, close all

%% Open the image and the segmentations

% Load the image
I = imread('/Users/ignaciorlando/Documents/MATLAB/artery-vein-classification/data/RITE-training/images/21_training.tif');
% Load the segmentation
segm = imread('/Users/ignaciorlando/Documents/MATLAB/artery-vein-classification/data/RITE-training/vessel-segmentations/21_training.png');

%% Generate a skeletonization

% Load the skeletonization
skel = bwmorph(segm, 'skel', Inf);

%% Extract the graph

% Generate graph of crossings
[G_crossings] = initializeGraphFromSkel_new(skel);

% And now transform the graph of crossings into a graph of segments
[Gout] = generateGraphOfSegments(G_crossings);

%% Extract unary and pairwise features

[Gout] = compute_pairwise_features(Gout, segm);