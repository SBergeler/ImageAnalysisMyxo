function [] = CheckSpotDetectionTimeLapse(fluor_paths, mat_path, xshift, yshift)
% CHECKSPOTDETECTIONTIMELAPSE: visually check the results from spot
% detection for time lapse images
%
% Copyright (c) 2021 Silke Bergeler
%
% Input: 
% - fluor_paths: paths to fluorescence images
% - mat_paths: path to mat file
% - xshift: how many pixels the meshes have to be shifted in x direction to
% match the fluorescence signals
% - yshift: how many pixels the meshes have to be shifted in y direction to
% match the fluorescence signals

masks = load(mat_path);

cell_ids = horzcat(masks.cellList.cellId{:});
cell_ids = unique(cell_ids);
ntcells = length(cell_ids);
ntframes = length(masks.cellList.meshData); % total number of frames

% determine the box size used for plotting the images 
boxes = zeros(ntcells,4);
for cell_ind = 1:ntcells
    box_x_min = [];
    box_x_max = [];
    box_y_min = [];
    box_y_max = [];
    for nframe = 1:ntframes 
        if any(masks.cellList.cellId{nframe} == cell_ids(cell_ind)) % if the cell_id is found for this frame
            ind = (masks.cellList.cellId{nframe} == cell_ids(cell_ind));
            box_x_min = [box_x_min masks.cellList.meshData{nframe}{ind}.box(1)];
            box_x_max = [box_x_max masks.cellList.meshData{nframe}{ind}.box(1) + masks.cellList.meshData{nframe}{ind}.box(3)];
            box_y_min = [box_y_min masks.cellList.meshData{nframe}{ind}.box(2)];
            box_y_max = [box_y_max masks.cellList.meshData{nframe}{ind}.box(2) + masks.cellList.meshData{nframe}{ind}.box(4)];
        end       
    end
    boxes(cell_ind,1) = min(box_x_min);
    boxes(cell_ind,2) = max(box_x_max);
    boxes(cell_ind,3) = min(box_y_min);
    boxes(cell_ind,4) = max(box_y_max);
end

for cell_ind = 1:ntcells
    disp(['cell ' num2str(cell_ind) ' out of ' num2str(ntcells) ' cells' ])
    stop = gui_spot_detection_timelapse({fluor_paths, mat_path, xshift, yshift, cell_ids(cell_ind), boxes(cell_ind,:)});
    if(stop)
        break;
    end
end
end

