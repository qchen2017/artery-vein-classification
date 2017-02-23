
function [delta] = lossCB(param, y, tildey)
% lossCB Compute the loss
% [delta] = lossCB(param, y, tildey)
% OUTPUT: delta: loss
% INPUT: param: parameters
%        y: ground truth labelling
%        tildey: estimated labelling

    % Retrieve delta
    delta = sum(y ~= tildey);

    %       incorrect veins          +     incorrect arteries
    %delta = sum((y ~= -1) ~= (tildey ~= -1)) + sum((y == 1) ~= (tildey == 1));
    
end