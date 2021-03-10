function [] = DeleteSignal0(mat_path)
% DELETESIGNAL0: delete the field signal0 for all cells in meshData, if
% existent (because some cells have this info and some don't, errors occur)
%
% Copyright (c) 2021 Silke Bergeler
%
% INPUT:
% - mat_path: path to mat file

masks = load(mat_path);
nframes = length(masks.cellList.meshData);

for frame = 1:nframes
    for ncell = 1:length(masks.cellList.meshData{frame}) 
        if any(strcmp(fieldnames(masks.cellList.meshData{frame}{ncell}),'signal0')) % test if signal0 is a field
            masks.cellList.meshData{frame}{ncell} = rmfield(masks.cellList.meshData{frame}{ncell},'signal0');
        end
    end
end

save(mat_path,'-struct','masks') % overwrite the mat file

end