function [] = PlotSpotObservablesOverTime(mat_paths, save_path, pixelsize, deltat)
% PLOTSPOTOBSERVABLESOVERTIME: Plot the spot intensity and the spot size
% over time for each cell
%
% Copyright (c) 2021 Silke Bergeler
%
% Input: 
% - mat_paths: paths to mat files
% - save_path: path to save figure and Excel file
% - pixelsize: pixel size of images
% - deltat: time step between frames (in min)

dataAll = cell(1,length(mat_paths));
for j = 1:length(mat_paths)
    masks = load(mat_paths{j});
    nframes = length(masks.cellList.meshData);
    data = []; % save dataset, cellid, frame, and spot observables
    
    for frame = 1:nframes
        ncells = length(masks.cellList.meshData{frame});
        for c = 1:ncells
            this_cell = masks.cellList.meshData{frame}{c};
            this_cell_id = masks.cellList.cellId{frame}(c);
            cell_length = this_cell.length;
            if isfield(this_cell,'spots_matlab')
                if length(this_cell.spots_matlab) == 1 && ~isempty(this_cell.spots_matlab.l) % only consider cells with one spot
                    this_cell_spots = this_cell.spots_matlab;
                    spot_major_axis = this_cell_spots.spot_major_axis;
                    spot_minor_axis = this_cell_spots.spot_minor_axis;
                    spot_ecc = this_cell_spots.spot_ecc;
                    spot_area = this_cell_spots.spot_area;
                    rel_orient = this_cell_spots.rel_orient;
                    rel_signal_int = this_cell_spots.rel_signal_int;
                    signal_int = this_cell_spots.signal_int;
                    cluster_position = this_cell_spots.l/this_cell.length;
                    data = [data; j this_cell_id frame spot_major_axis spot_minor_axis spot_ecc ...
                        spot_area rel_orient rel_signal_int signal_int cell_length cluster_position];
                else
                    data = [data; j this_cell_id frame 0 0 0 0 0 0 0 cell_length NaN];
                end
            else
                data = [data; j this_cell_id frame 0 0 0 0 0 0 0 cell_length NaN];
            end
        end
    end
    dataAll{j} = data;
end

% concatenate data
data = vertcat(dataAll{:});

data = array2table(data, 'VariableNames',{'data_set','cellid','frame_no','spot_major_axis' ...
    'spot_minor_axis' 'spot_ecc' 'spot_area' 'rel_orient' 'rel_signal_int' 'signal_int' 'cell_length' 'cluster_position'});

% convert units from pixel to um
data.spot_major_axis = data.spot_major_axis.*pixelsize;
data.spot_minor_axis = data.spot_minor_axis.*pixelsize;
data.spot_area = data.spot_area.*(pixelsize^2);
data.cell_length = data.cell_length.*pixelsize;

% map the cluster positions (relative to cell length to the interval
% [0,0.5]
data.cluster_position = arrayfun(@Mapping,data.cluster_position);

%%
figure;
set(gcf,'Unit','Inches','position',[0,0,7,7],'PaperUnits', 'Inches', 'PaperSize', [7, 7])
subplot(2,2,1);

for d=1:length(unique(data.data_set))
    data_sel = data(data.data_set == d,:);
    cellids = unique(data_sel.cellid);
    for i=1:length(cellids)
        celllengths = data_sel.cell_length(data_sel.cellid == cellids(i));
        frames = data_sel.frame_no(data_sel.cellid == cellids(i));
        times = (frames - min(frames))*deltat; % set all trajectories to the same starting point
        plot(times, celllengths)
        hold on
    end
end
xlabel('time [min]');
ylabel('cell length [\mum]');

subplot(2,2,2);
for d=1:length(unique(data.data_set))
    data_sel = data(data.data_set == d,:);
    cellids = unique(data_sel.cellid);
    for i=1:length(cellids)
        signalints = data_sel.signal_int(data_sel.cellid == cellids(i));
        frames = data_sel.frame_no(data_sel.cellid == cellids(i));
        times = (frames - min(frames))*deltat; % set all trajectories to the same starting point
        plot(times, signalints)
        hold on
    end
end
xlabel('time [min]');
ylabel('cluster signal intensity');

subplot(2,2,3);
for d=1:length(unique(data.data_set))
    data_sel = data(data.data_set == d,:);
    cellids = unique(data_sel.cellid);
    for i=1:length(cellids)
        areas = data_sel.spot_area(data_sel.cellid == cellids(i));
        frames = data_sel.frame_no(data_sel.cellid == cellids(i));
        times = (frames - min(frames))*deltat; % set all trajectories to the same starting point
        plot(times, areas)
        hold on
    end
end
xlabel('time [min]');
ylabel('spot area [\mum^2]');

% save the figure as pdf
saveas(gcf,[save_path 'spot_observables_over_time.pdf'])
close(gcf)

% save the data to excel file
description = table({'data_set' 'data set (important if several data sets are combined)';...
    'cellid' 'unique id of the cell'; ... 
    'frame_no' 'number of the frame (starting from 1, which might be different to the names of the images)'; ...
    'spot_major_axis' 'major axis of ellipse fitted to spot'; ...
    'spot_minor_axis' 'minor axis of ellipse fitted to spot'; ...
    'spot_ecc' 'eccentricity of ellipse fitted to spot'; ...
    'spot_area' 'area of the spot [um^2]'; ...
    'rel_orient' 'orientation of spot relative to orientation of cell center line [degree]'; ...
    'rel_signal_int' 'signal intensity in the spot relative to the intensity in the cell'; ...
    'signal_int' 'absolute signal intensity in the spot'; ...
    'cell_length' 'cell length [um]'; ...
    'cluster_position' 'position of cluster relative to cell length'
    });

% write table and description of variables to comma-separated txt file
writetable(data,[save_path 'data_spot_observables_and_cell_length.xls'],'Sheet',1); 
writetable(description,[save_path 'data_spot_observables_and_cell_length.xls'],'Sheet',2,'WriteVariableNames',false); 

end

