
function [new_segm, onh_mask] = removeDisconnectedStructures(segm, onh_mask)
% REMOVEDISCONNECTEDSTRUCTURES Takes only the major connected component of
% the binary segmentation.
%   Detailed explanation goes here

    % remove error in the central vessel reflex
    segm = imclose(segm, strel('disk', 2, 8));
    % resize the od masks in case the image was previously rescaled
    if (size(onh_mask, 2) > size(segm,2))
        onh_mask = imresize(onh_mask, size(segm));
    end
    % dilate the od mask to remove undisered false positives
    onh_mask = imdilate(onh_mask, strel('disk',5,8));
    % identify od pixels
    od_idx = find(onh_mask);
    % find connected components of the skeletonization
    CC = bwconncomp(segm);
    % preserve only the structures that are connected to the onh
    new_segm = false(size(segm));
    for i = 1 : CC.NumObjects
        if ~isempty(intersect(CC.PixelIdxList{i}, od_idx))
            new_segm(CC.PixelIdxList{i}) = true;
        end
    end
    % remove vessels inside the onh
    new_segm = new_segm .* imcomplement(onh_mask);

end