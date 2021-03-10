function [new_meshData] = deleteFollowingCells(cell_ind, nframe, new_meshData, ntframes)
% DELETEFOLLOWINGCELLS: deletes the entries for the cell with index 'cell_ind' in all
% frames following frame 'nframe'
%
% Copyright (c) 2021 Silke Bergeler
%
% Input:
% - cell_ind: index of the cell to be deleted
% - nframe: entries for this cell are deleted for all frames following this
% frame
% - new_meshData: mesh data 
% - ntframes: total number of frames
% Output:
% - new_meshData: updated mesh data set

for frame = (nframe+1):ntframes
    if length(new_meshData{frame}) >= cell_ind
        new_meshData{frame}(cell_ind) = {[]};
    end
end
end