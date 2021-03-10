function [PH_paths, fluor_paths] = CheckCells(mat_path, PH_paths, fluor_paths, division, timelapse)
% CHECKCELLS: check the cells identified by Oufti and delete entries in the
% mat file with an empty mesh or with a mesh that contains infs;
% furthermore check if there are frames with no cells detected at the end/beginning of the series
% and delete them if this is the case
%
% Copyright (c) 2021 Silke Bergeler
%
% INPUT:
% - mat_path: path to mat file
% - PH_paths: paths to PH images
% - fluor_paths: paths to fluorescence images
% - division: boolean, if true, cells divide
% - timelapse: boolean, if true, time lapse images are considered
% OUTPUT:
% - PH_paths: updated paths to PH images
% - fluor_paths: updated paths to fluorescence images

count  = 0;

masks = load(mat_path); % load the mat file 
nframes = length(masks.cellList.meshData);
for frame = 1:nframes
   ncells = length(masks.cellList.meshData{frame});
   for cell = 1:ncells
       mesh = masks.cellList.meshData{frame}{cell}.mesh;
       if size(mesh,2) ~= 4 || any(any(isinf(mesh)))
           count = count + 1;
           masks.cellList.meshData{frame}{cell} = [];
           masks.cellList.cellId{frame}(cell) = NaN;
       end
   end
end

% delete the NaNs from cellId and the [] from meshData
for frame = 1:nframes
    if(~isempty(masks.cellList.meshData{frame}))
        masks.cellList.meshData{frame} = masks.cellList.meshData{frame}(~cellfun('isempty',masks.cellList.meshData{frame}));
        masks.cellList.cellId{frame} = masks.cellList.cellId{frame}(~isnan(masks.cellList.cellId{frame}));
    end
end

if (timelapse && division) || ~timelapse
    % check if there are frames with no cells detected at the end of the series
    % and delete them if this is the case
    vecempty = false(1,nframes);
    for frame = 1:nframes
        vecempty(frame) = isempty(masks.cellList.meshData{frame});
    end
    ind = find(~vecempty,1,'last'); % last index of the frame that has at least one cell detected
    masks.cellList.meshData(ind+1:end) = []; % delete frames with no cells
    masks.cellList.cellId(ind+1:end) = []; % delete frames with no cells
    
    PH_paths = PH_paths(1:ind);
    fluor_paths = fluor_paths(1:ind);
    
    % check if there are frames with no cells detected in the beginning of the series
    % and delete them if this is the case
    ind = find(~vecempty,1,'first'); % first index of the frame that has at least one cell detected
    masks.cellList.meshData(1:ind-1) = []; % delete frames with no cells
    masks.cellList.cellId(1:ind-1) = []; % delete frames with no cells
    
    PH_paths = PH_paths(ind:end);
    fluor_paths = fluor_paths(ind:end);
end
           
save(mat_path,'-struct','masks') % overwrite the mat file
disp(['CheckCells: ' num2str(count) ' cells are deleted because of invalid mesh'])
end
