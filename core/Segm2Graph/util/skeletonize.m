
function [skel, segm, onh_perimeter, onh_mask] = skeletonize(segm, onh_mask)

    % take the perimeter of the onh
    onh_perimeter = imdilate(bwperim(onh_mask, 8), strel('disk',2,8));
    % remove it from the onh mask
    onh_mask(onh_perimeter==1) = 0;
    % get the skeleton
    
    
end