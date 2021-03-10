function [box_out] = ChangeBoxSize(box_in, scaling)
% CHANGEBOXSIZE: Change the size of the box by scaling it with a factor
% 'scaling'
%
% Copyright (c) 2021 Silke Bergeler
%
% Input: 
% - box_in: parameters of box to be changed in size
% - scaling: scaling factor
% Output:
% - box_out: parameters of rescaled box
%
% NOTE: the created box might be outside of the image range, but the
% overhanging regions are cropped in SpotDetection2DArea (however, it might
% be better to implement this here)

topx = box_in(1);
topy = box_in(2);
width = box_in(3);
height = box_in(4);

width_out = round(width*sqrt(scaling));
height_out = round(height*sqrt(scaling));
topx_out = round(topx - (width_out-width)/2);
topy_out = round(topy - (height_out-height)/2);

box_out = [topx_out, topy_out, width_out, height_out];
end

