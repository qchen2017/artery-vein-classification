
clear, clc, close all

% Load the segmentation
segm = imread('/Users/ignaciorlando/Documents/MATLAB/artery-vein-classification/data/RITE-training/vessel-segmentations/21_training.png');

% Load the skeletonization
skel = bwmorph(segm,'skel',Inf);

% Generate graph of crossings
[G_crossings] = initializeGraphFromSkel_new(skel);

% And now transform the graph of crossings into a graph of segments
[Gout] = generateGraphOfSegments(G_crossings);