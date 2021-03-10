function [] = WriteDataSpotSplittingToFile(mat_paths, save_path, pixelsize)
% WRITEDATASPOTSPLITTINGTOFILE: Write data for cluster splitting to a file
%
% Copyright (c) 2021 Silke Bergeler
% 
% Input: 
% - mat_paths: path to mat files
% - save_path: path to save Excel file
% - pixelsize: pixel size of images

dataAll = cell(1,length(mat_paths));
for j = 1:length(mat_paths)
    masks = load(mat_paths{j});
    nframes = length(masks.cellList.meshData);
    
    % get the infos of all daughter cells in the frame upon division (that's
    % when there is an entry 'ancestors' in the meshData)
    % only take the cells with one spot detected
    
    meshes_daughters = [];
    for frame = 1:nframes
        for ncell = 1:length(masks.cellList.meshData{frame})
            if ~isempty(masks.cellList.meshData{frame}{ncell}.ancestors) && ...
                    length(masks.cellList.meshData{frame}{ncell}.spots_matlab) == 1 && ...
                    ~all(structfun(@(x) isempty(x), masks.cellList.meshData{frame}{ncell}.spots_matlab))
                mesh = masks.cellList.meshData{frame}{ncell};
                mesh.birthframe = frame; % update the info from Oufti here
                meshes_daughters = [meshes_daughters, mesh];
            end
        end
    end
    
    % take only the data where there are two daughters of the same ancestor
    ancestors = arrayfun(@(x) meshes_daughters(x).ancestors, 1:length(meshes_daughters));
    [count,val] = hist(ancestors,unique(ancestors));
    selected = val(count == 2);
    if isempty(selected)
        break
    end
    selected_ind = ismember(ancestors, selected);
    meshes_daughters = meshes_daughters(1,selected_ind);
    
    % create a table with the area and intensity of the cluster in each daughter cell, the
    % min and the max value and the asymmetry value, as well as the area and intensity of the
    % cluster in the mother cell before division
    data = [];
    count = 1;
    for i = 1:length(selected) % all ancestor values that appear exactly twice
        data_area_tmp = [];
        data_int_tmp = [];
        for l = 1:length(meshes_daughters)
            if meshes_daughters(l).ancestors == selected(i)
                data_area_tmp = [data_area_tmp, meshes_daughters(l).spots_matlab.spot_area*pixelsize^2];
                data_int_tmp = [data_int_tmp, meshes_daughters(l).spots_matlab.signal_int];
                birthframe = meshes_daughters(l).birthframe;
            end
        end
        spot_area_mother = [masks.cellList.meshData{birthframe-1}{masks.cellList.cellId{birthframe-1} == selected(i)}.spots_matlab.spot_area];
        spot_int_mother = [masks.cellList.meshData{birthframe-1}{masks.cellList.cellId{birthframe-1} == selected(i)}.spots_matlab.signal_int];
        if(length(spot_area_mother) == 1) % only consider cases where there is only one cluster in the mother cell (it might happen that the cluster divides before cell division is detected)
            spot_area_mother = spot_area_mother*pixelsize^2;
            data(count,:) = [selected(i),data_area_tmp(1),data_area_tmp(2),spot_area_mother, data_int_tmp(1),data_int_tmp(2),spot_int_mother];
            count = count + 1;
        end
    end  
    dataAll{j} = data;
end

% concatenate data
data = vertcat(dataAll{:});

if(isempty(data))
    disp('no data to write to excel file')
    return
end

data = array2table(data, 'VariableNames',{'ancestor','area1','area2','area_mother','int1','int2','int_mother'});

description = table({'ancestor' 'id of the cell from which the daughters descend'; ...
    'area1' 'area of the cluster of daughter cell 1 [um^2]'; ...
    'area2' 'area of the cluster of daughter cell 2 [um^2]'; ...
    'area_mother' 'area of the cluster of the mother cell [um^2]'; ...
    'int1' 'absolute signal intensity in the spot of daughter cell 1'; ...
    'int2' 'absolute signal intensity in the spot of daughter cell 2'; ...
    'int_mother' 'absolute signal intensity in the spot of the mother cell'});

% write table and description of variables to comma-separated txt file
writetable(data,[save_path 'data_cluster_splitting.xls'],'Sheet',1); 
writetable(description,[save_path 'data_cluster_splitting.xls'],'Sheet',2,'WriteVariableNames',false); 

end

