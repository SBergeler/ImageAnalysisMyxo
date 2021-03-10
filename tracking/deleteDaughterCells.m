function new_meshData = deleteDaughterCells(cell_ind, new_meshData, ntframes)
% DELETEDAUGHTERCELLS: delete the assignment of the daughter cell(s) to the
% mother cell, if selected by visual inspection
%
% Copyright (c) 2021 Silke Bergeler
%
% Input: 
% - cell_ind: cell index of the cell for which the assignment with the
% daughter cells should be deleted
% - new_meshData: mesh data
% - ntframes: total number of frames
% Output:
% - new_meshData: updated mesh data

for frame = 1:ntframes
    for ncell = 1:length(new_meshData{frame})
        if ~isempty(new_meshData{frame}{ncell}) && ~isempty(new_meshData{frame}{ncell}.ancestors)
            % if the ancestor is the cell with index cell_ind, delete this entry
            if new_meshData{frame}{ncell}.ancestors == cell_ind
                new_meshData{frame}{ncell}.ancestors = [];
            end
        end
    end
end
end

