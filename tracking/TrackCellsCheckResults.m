function new_meshData = TrackCellsCheckResults(new_meshData, PH_paths)
% TRACKCELLSCHECKRESULTS: visually check the results from cell tracking
% and delete frames / cells that are not correctly tracked
%
% Copyright (c) 2021 Silke Bergeler
%
% Input:
% - new_meshData: mesh data obtained after tracking
% - PH_paths: paths to PH images
% Output:
% - new_meshData: updated mesh data

% maximal number of cells considered (it might be that tracking is lost and found again, then these cells are considered as two cells)
ntcells = max(cellfun(@(x) length(x), new_meshData));
% total number of frames
ntframes = length(new_meshData); 

% determine the box size used for plotting the images 
boxes = zeros(ntcells,4);
for cell_ind = 1:ntcells
    box_x_min = [];
    box_x_max = [];
    box_y_min = [];
    box_y_max = [];
    for nframe = 1:ntframes 
        if cell_ind <= length(new_meshData{nframe}) && ~isempty(new_meshData{nframe}{cell_ind})
            box_x_min = [box_x_min new_meshData{nframe}{cell_ind}.box(1)];
            box_x_max = [box_x_max new_meshData{nframe}{cell_ind}.box(1) + new_meshData{nframe}{cell_ind}.box(3)];
            box_y_min = [box_y_min new_meshData{nframe}{cell_ind}.box(2)];
            box_y_max = [box_y_max new_meshData{nframe}{cell_ind}.box(2) + new_meshData{nframe}{cell_ind}.box(4)];
        end       
    end
    boxes(cell_ind,1) = min(box_x_min);
    boxes(cell_ind,2) = max(box_x_max);
    boxes(cell_ind,3) = min(box_y_min);
    boxes(cell_ind,4) = max(box_y_max);
end

% save the frame where tracking should be stopped (because something goes wrong)
stop_frames = zeros(ntcells,1,'int8'); 
delete_daughters = false(ntcells,1);
for cell_ind = 1:ntcells
    disp(['cell ' num2str(cell_ind) ' out of ' num2str(ntcells) ' cells' ])
    [stop_frames(cell_ind),delete_daughters(cell_ind)] = gui_cell_tracking({PH_paths, new_meshData, cell_ind, boxes});
end

%% delete the data in 'new_meshData' that should be excluded from the data analysis (as obtained from the visual inspections)
for cell_ind = 1:ntcells
    nframe = stop_frames(cell_ind);
    new_meshData = deleteFollowingCells(cell_ind, nframe, new_meshData, ntframes);
    if delete_daughters(cell_ind)
        new_meshData = deleteDaughterCells(cell_ind, new_meshData, ntframes);
    end
end
end

