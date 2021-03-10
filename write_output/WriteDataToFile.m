function celldataAll = WriteDataToFile(mat_paths, save_path, pixelsize, timelapse)
% WRITEDATATOFILE: generate table with the informations about the cells and
% the spots and write it to an excel file
%
% Copyright (c) 2021 Silke Bergeler
%
% Input:
% - mat_paths: path to mat files
% - save_path: path to save Excel file
% - pixelsize: pixel size of images
% - timelapse: Boolean, if true, time lapse images are considered
%
% Output:
% - celldataAll: cell array with data about cells and spots

disp('write data to file ...')

celldataAll = cell(1,length(mat_paths));

for j = 1:length(mat_paths)
    load(mat_paths{j},'cellList');
    
    % create array with the frame numbers
    nframe = cellfun(@(x) length(x), cellList.cellId); % get number of cells per frame
    frame_no = arrayfun(@(x) repmat(x,1,nframe(x)),1:length(nframe),'UniformOutput',false);
    frame_no = cell2mat(frame_no)';
    
    % create array with the cell ids in all frames
    cellid = cell2mat(cellList.cellId)';
    
    % create a table with the cell (+ spot infos if existent)
    celldata = horzcat(cellList.meshData{:}); % cell array with cells of all frames
    
    % check if all cells in celldata have same fields and otherwise delete this field
    fields = fieldnames(celldata{1});
    fields_diff = [];
    for c = 1:length(celldata)
        fields_tmp = fieldnames(celldata{c});
        diff1 = setdiff(fields, fields_tmp);
        diff1 = cellfun(@(x) string(x), diff1)';
        diff2 = setdiff(fields_tmp,fields);
        diff2 = cellfun(@(x) string(x), diff2)';
        if(~isempty(diff1) || ~isempty(diff2))
            fields_diff = [fields_diff, diff1, diff2];
        end
    end
    
    fields_diff = unique(fields_diff);
    for i = 1:length(fields_diff)
        disp(['WARNING: field ' char(fields_diff(i)) ' is deleted because not all cells have this field!'])
    end
    
    % delete all fields in fields_diff
    for f = 1:length(fields_diff)
        for c = 1:length(celldata)
            if any(strcmp(fieldnames(celldata{c}),char(fields_diff(f)))) % test if the entry in fields_diff is one of the fields in celldata{c}
                celldata{c} = rmfield(celldata{c},char(fields_diff(f)));
            end
        end
    end
    
    celldata = cat(1,celldata{:}); % create one struct with all cell infos (only works if the fields for each cell are the same!)
    celldata = struct2table(celldata); % convert to a table
    
    if ismember('spots_matlab', celldata.Properties.VariableNames)
        % if there is exactly one spot per cell for all detected cells or only one cell with spots
        % detected, celldata.spots_matlab is a struct, otherwise a cell array
        % with structs of possibly different sizes as entries; the two cases need to be
        % considered separately
        if isstruct(celldata.spots_matlab)
            spot_data = celldata.spots_matlab;
        else
            spot_data = cellfun(@(x) ConvertSpotData(x),celldata.spots_matlab);
        end
        spot_data = struct2table(spot_data);
        spot_data = spot_data(:,{'l','d','x','y','spot_major_axis',...
            'spot_minor_axis','spot_ecc','spot_area','signal_int','rel_signal_int','rel_orient'}); % select only the columns of interest
    else
        spot_data = table();
    end
    
    celldata_variables = celldata.Properties.VariableNames;
    if (timelapse)
        variables_sel = {'ancestors', 'length', 'area', 'width'};
    else
        variables_sel = {'length', 'area', 'width'};
    end
    if(any(strcmp(celldata_variables,'threshold_spot')))
        variables_sel{end+1} = 'threshold_spot';
    end
    
    if(any(strcmp(celldata_variables,'background')))
        variables_sel{end+1} = 'background';
    end
    
    if(any(strcmp(celldata_variables,'cell_signal_int')))
        variables_sel{end+1} = 'cell_signal_int';
    end
    
    celldata = celldata(:,variables_sel); % select only the columns of interest
    data_set = repmat(j,length(frame_no),1);
    celldata = horzcat(table(data_set),table(frame_no),table(cellid), celldata, spot_data);
    
    % change the units from pixel to um
    celldata(:,{'length','width'}) = ...
        array2table(table2array(celldata(:,{'length','width'})).*pixelsize);
    celldata(:,{'area'}) = ...
        array2table(table2array(celldata(:,{'area'})).*pixelsize^2);
    if ismember('spot_area', celldata.Properties.VariableNames)
        data = table2array(celldata(:,{'spot_area'}));
        if iscell(data)
            celldata(:,{'spot_area'}) = ...
                array2table(arrayfun(@(x) cellfun(@(y) y.*pixelsize^2,x,'UniformOutput',false),...
                table2array(celldata(:,{'spot_area'}))));
            celldata(:,{'l','d','spot_major_axis','spot_minor_axis'}) = ...
                array2table(arrayfun(@(x) cellfun(@(y) y.*pixelsize,x,'UniformOutput',false),...
                table2array(celldata(:,{'l','d','spot_major_axis','spot_minor_axis'}))));
        else
            celldata(:,{'spot_area'}) = ...
                array2table(arrayfun(@(x) x.*pixelsize^2, table2array(celldata(:,{'spot_area'}))));
            celldata(:,{'l','d','spot_major_axis','spot_minor_axis'}) = ...
                array2table(arrayfun(@(x) x.*pixelsize, table2array(celldata(:,{'l','d','spot_major_axis','spot_minor_axis'}))));
        end
        % add column for the diffuse signal in the cell
        if iscell(celldata.signal_int)
            celldata.diffuse_signal_int = celldata.cell_signal_int - cellfun(@sum, celldata.signal_int);
        else
            celldata.diffuse_signal_int = celldata.cell_signal_int - celldata.signal_int;
        end
    end
    celldataAll{j} = celldata;
end

% concatenate celldata tables to one
celldataAll = vertcat(celldataAll{:});

description = table({'data_set' 'data set (important if several data sets are combined)'; ...
    'frame_no' 'number of the frame (starting from 1, which might be different to the names of the images)'; ...
    'cellid' 'unique id of the cell'; ...
    'ancestors' 'id of ancestor cell when the cell just divided'; ...
    'length' 'cell length [um]'; ...
    'area' 'cell area [um^2]';...
    'width' 'cell width [um] (average cell width using widths between 30% and 70% cell length)'; ...
    'threshold_spot' 'threshold intensity for spots (without background correction)'; ...
    'background' 'average intensity outside the cell in a surrounding box'; ...
    'cell_signal_int' 'total intensity inside the cell (background corrected)'; ...
    'l' 'position of spot along center line [um] (if there are subindices, i.e. l_1, ... then more than one spot is found in a cell and the values refer to each spot)'; ...
    'd' 'position of spot perpendicular to center line [um]'; ...
    'x' 'position of spot in x direction on fluorescence image [pixel]';...
    'y' 'position of spot in y direction on fluorescence image [pixel]';...
    'spot_major_axis' 'major axis of ellipse fitted to spot'; ...
    'spot_minor_axis' 'minor axis of ellipse fitted to spot'; ...
    'spot_ecc' 'eccentricity of ellipse fitted to spot'; ...
    'spot_area' 'area of the spot [um^2]'; ...
    'signal_int' 'absolute signal intensity in the spot (background corrected)'; ...
    'rel_signal_int' 'signal intensity in the spot relative to the intensity in the cell (using background corrected values)'
    'rel_orient' 'orientation of spot relative to orientation of cell center line [degree]';
    'diffuse_signal_int' 'total signal intensity in the cell (excluding the spots if spots are detected; background corrected)'});

% write table and description of variables to comma-separated txt file
% delete existing file because writetable does not overwrite file
if exist([save_path 'cell_data.xls'],'file')
    delete([save_path 'cell_data.xls'])
end
writetable(celldataAll,[save_path 'cell_data.xls'],'Sheet',1); 
writetable(description,[save_path 'cell_data.xls'],'Sheet',2,'WriteVariableNames',false); 

disp('done')
end
