%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use this script to plot the results of the image analysis script (it
% is possible to combine several data sets)
%
% Copyright (c) 2021 Silke Bergeler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear workspace
clear(); close all;

%% folders with data to visualize (NEED TO BE ADJUSTED)
paths = {'PATH_TO_FOLDER1','PATH_TO_FOLDER2'}; % can be a cell array with paths to different folders

%% get data
n_datasets = length(paths); % number of data sets to combine

% create folder to save the plots and excel files (in the folder of the
% first entry in paths)
save_path = [paths{1} 'image_analysis_combined/'];
if ~exist(save_path, 'dir') 
    mkdir(save_path)
end

% load parameters that were used for image analysis
info_params = cell(1,n_datasets);
for i = 1:n_datasets
    info_params_tmp = load([paths{i} 'info_params.mat'], 'parameters_image_analysis');
    info_params{i} = info_params_tmp.parameters_image_analysis;
end
info_params = cell2mat(info_params);

% check if all relevant image analysis parameters are the same for all data sets 
fnames = fieldnames(info_params);
fnames = fnames(1:6);

info_params_selected = struct([]);
for i = 1:n_datasets
    for j = 1:length(fnames)
        info_params_selected(i).(fnames{j}) = info_params(i).(fnames{j});
    end
end

check = 0;
for i = 1:(n_datasets-1)
    check = check + ~isequal(info_params_selected(i),info_params_selected(i+1));
end
if check > 0
    disp('The data sets that should be combined do not have the same image analysis parameters! --> stop execution')
    return
end

% define parameters based on info mat files
pixelsize = info_params(1).pixelsize;
deltat = info_params(1).deltat;
timelapse = info_params(1).timelapse;
fluorescence_images = info_params(1).fluorescence_images;
species = info_params(1).species;
division = info_params(1).division;
xshift = arrayfun(@(i) info_params(i).xshift,1:length(info_params),'UniformOutput',false);
yshift = arrayfun(@(i) info_params(i).yshift,1:length(info_params),'UniformOutput',false);

mat_paths = cell(1,n_datasets);
for i = 1:n_datasets
    full_path = info_params(i).mat_path;
    [pathstr,name,ext] = fileparts(full_path);
    mat_paths{i} = [paths{i} name ext];
end

fluor_paths = cell(1,n_datasets);
if (fluorescence_images)
    for i = 1:n_datasets
        full_paths = info_params(i).fluor_paths;
        fluor_paths_tmp = cell(length(full_paths),1);
        for j = 1:length(full_paths)
            full_path = full_paths{j};
            [pathstr,name,ext] = fileparts(full_path);
            fluor_paths_tmp{j} = [paths{i} 'Fluor/' name ext];
        end
        fluor_paths{i} = fluor_paths_tmp;
    end
end

%% get data to visualize
% Write data xls file (does not contain the infos about spot tracking for division = false --> this info is saved in a separate file)
celldataAll = WriteDataToFile(mat_paths, save_path, pixelsize, timelapse);

%% Visualization - Cell statistics and growth
if ~timelapse
    PlotCellStatistics(celldataAll, save_path)
end

if(~fluorescence_images) 
    return
end

%% Visualization (if fluorescence images exist)
if timelapse && division
    % plot spot sizes after splitting
    PlotSpotSplittingResults(mat_paths, save_path, pixelsize)
    WriteDataSpotSplittingToFile(mat_paths, save_path, pixelsize)
    % plot spot observables
    PlotSpotObservablesOverTime(mat_paths, save_path, pixelsize, deltat)
    if ismember(species, ['X','Y','Z'])
        str = input('Shall movies be generated? [yes/no]','s');
        if strcmp(str, 'yes')
            % plot cluster trajectories, changes of size and signal intensity over time
            PlotClusterTrajectoryAndProperties(mat_paths, fluor_paths, save_path, pixelsize, xshift, yshift);
        end
    end
elseif ~timelapse % SNAPSHOTS 
    % histogram of number of spots per cell
    CountSpots(mat_paths, save_path)
    % only if at least 2 spots exist in the data set, plot the spot
    % statistics / density plots
    if ismember('l',celldataAll.Properties.VariableNames) 
        if iscell(celldataAll.l) && sum(cellfun(@(x) ~isempty(x),celldataAll.l)) > 2 || ...
                ~iscell(celldataAll.l) && length(celldataAll.l) > 2
            % plot histogram of spot observables
            PlotSpotStatistics(mat_paths, save_path, pixelsize)
            % plot spot density in the cell
            if ~ismember(species, ['X','Y','Z'])
                % spot positions as a heat map (arbitrary number of spots)
                DensityPlotSpots(mat_paths, save_path, pixelsize);
            end
        end
    end
    % plot demograph
    % the second to last parameter is the maximal number of cells used for the plot
    % the last parameter is the percentage of pixels in the demograph which are
    % plotted with highest intensity value
    PlotDemograph(mat_paths, fluor_paths, save_path, xshift, yshift, pixelsize, 500, 3)
end
