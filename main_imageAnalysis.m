%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use this script to perform the image analysis (including data processing, 
% cell tracking and spot detection, if applicable)
%
% Copyright (c) 2021 Silke Bergeler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear workspace
clear(); close all;

% add folders to mat path
addpath(genpath(pwd))

%% Settings (NEED TO BE ADJUSTED)
% path to folder with tif images (unstacked) and mat file from Oufti
path = 'PATH_TO_FOLDER';

% size of one pixel (in um)
pixelsize = 0.065;

% time step between frames (in min)
deltat = 30;

%-- OPTIONS --:
% whether time-lapse images are considered or not (cell segmentation in
% Oufti does not have to be done in the time-lapse mode in order to choose this option)
timelapse = false; 
% whether fluorescence images are also considered or only phase contrast
% images
fluorescence_images = true; 
% if fluorescence images exist, which proteins are considered (either 'X' =
% PomX, 'Y' = PomY or 'Z' = PomZ, 'B' = ParB) to allow for differences in
% the spot detection parameters (see below)
species = 'Y';
% whether cells divide or not
% if division is set to false: for time lapse images of cells that do not 
% divide (cephalexin-treated cells), the cluster size is recorded over time 
% (only the cell outlines in the first frame are needed from Oufti)
% Note: the case timelapse = true && division = false is still work in progress!
division = false;

% If fluorescence images exist, set the parameters for spot detection:
% - threshold value: X for method 1 or 2 
% - method: method for setting threshold
%   1: threshold = mean(cytoplasm) + X*std(cytoplasm)
%   2: threshold = X*mean(cytoplasm));
% - npixel: number of connected pixels to be called a spot

if species == 'X'
    threshold_value = 5;
    method = 1;
    npixel = 25;
elseif species == 'Y'
    threshold_value = 4;
    method = 1;
    npixel = 10;
elseif species == 'Z'
    threshold_value = 1.5;
    method = 1;
    npixel = 20;
else
    threshold_value = 1.5;
    method = 1;
    npixel = 20;
end

%% Data preprocessing
% get the name of the .mat file (if there is only one .mat file in the
% folder); if a mat file with the ending '_tracked.mat' already exists delete the file

disp('data preprocessing ...')

% create folders for images and image stacks as required for further analysis
MakeFileStructure(path)

mat_files = dir([path,'*.mat']);

% check if a file 'info_params.mat' exist and delete it from the list of
% mat files
indices = arrayfun(@(i) ~strcmp('info_params.mat',mat_files(i).name),1:length(mat_files));
mat_files = mat_files(indices);

if isempty(mat_files)
    disp('Error: No mat file in the directory!')
    return
elseif length(mat_files) == 1 % if there is only one mat file, consider this one
    mat_file = mat_files;
else
    mat_names = {mat_files.name};
    % find the common part in the names
    mat_name = mat_names{1}(all(~diff(char(mat_names(:))),1));
    if ~isempty(dir([path,'*_tracked.mat']))
        delete([path, '*_tracked.mat'])
    end
    if ~isempty(dir([path,'*_original.mat']))
        delete([path mat_name '.mat'])
        movefile([path mat_name '_original.mat'], [path mat_name '.mat'])
    end
    mat_file = dir([path mat_name '.mat']);
end

mat_path = [path,mat_file.name];

% delete the fields 'signal0' in MeshData (from Oufti) because this can lead to errors
DeleteSignal0(mat_path);

%%
global IMAGEDIMX 
global IMAGEDIMY 
if fluorescence_images
    fluor_dir = dir([path,'Fluor/*.tif']);
    if contains(fluor_dir(1).name,"_RAW_ch0")
        fluor_names = arrayfun(@(x) x.name, fluor_dir, 'UniformOutput',false);
        fluor_paths = cellfun(@(x) [path,'Fluor/',x], fluor_names, 'UniformOutput',false);
    else
        for i = 1:length(fluor_dir)
            % split the name (name + format)
            name_splitted = strsplit(fluor_dir(i).name,'.');
            % select the name
            name_splitted = name_splitted{1};
            % test if the last character is a number
            name_digits = isstrprop(name_splitted,'digit');
            if ~name_digits(end)
                fluor_dir_sep = fluor_dir(i); % save i-th entry separately
                fluor_dir(i) = []; % delete i-th entry
                continue; % there should be at most one entry without a number at the end
            end
        end
        fluor_names = sort_nat(arrayfun(@(x) x.name, fluor_dir, 'UniformOutput',false)); % sort names of files in natural order
        if exist('fluor_dir_sep','var')
            % add the name without a number as first entry
            C = {fluor_dir_sep.name, fluor_names};
            C = vertcat(C{:});
            fluor_paths = cellfun(@(x) [path,'Fluor/',x], ...
                C, 'UniformOutput',false);
        else
            fluor_paths = cellfun(@(x) [path,'Fluor/',x], fluor_names, 'UniformOutput',false);
        end
    end
    % define image dimensions as global variable (such that we can use it
    % later)
    image = imread(fluor_paths{1});
    dims = size(image);
    IMAGEDIMX = dims(2);
    IMAGEDIMY = dims(1);
end

PH_dir = dir([path,'PH/*.tif']);
if contains(PH_dir(1).name,"_RAW_ch0")
    PH_names = arrayfun(@(x) x.name, PH_dir, 'UniformOutput',false);
    PH_paths = cellfun(@(x) [path,'PH/',x], PH_names, 'UniformOutput',false);
else
    for i = 1:length(PH_dir)
        % split the name (name + format)
        name_splitted = strsplit(PH_dir(i).name,'.');
        % select the name
        name_splitted = name_splitted{1};
        % test if the last character is a number
        name_digits = isstrprop(name_splitted,'digit');
        if ~name_digits(end)
            PH_dir_sep = PH_dir(i); % save i-th entry separately
            PH_dir(i) = []; % delete i-th entry
            continue; % there should be at most one entry without a number at the end
        end
    end
    PH_names = sort_nat(arrayfun(@(x) x.name, PH_dir, 'UniformOutput',false)); % sort names of files in natural order
    if exist('PH_dir_sep','var')
        % add the name without a number as first entry
        C = {PH_dir_sep.name, PH_names};
        C = vertcat(C{:});
        PH_paths = cellfun(@(x) [path,'PH/',x], ...
            C, 'UniformOutput',false);
    else
        PH_paths = cellfun(@(x) [path,'PH/',x], PH_names, 'UniformOutput',false);
    end
end

% rename the original mat_file to keep it
mat_name = strsplit(mat_file.name,'.mat');
mat_path_new = [path mat_name{1} '_original.mat'];
copyfile(mat_path,mat_path_new)

% delete cells that do not have a mesh
if exist('fluor_paths','var')
    [PH_paths, fluor_paths] = CheckCells(mat_path, PH_paths, fluor_paths, division, timelapse);
else
    PH_paths = CheckCells2(mat_path, PH_paths, division, timelapse);
end
    
% add length, width and area of bacterial cells to mat file
AddLengthArea(mat_path);

%% Cell tracking for time-lapse images 
% Since cell outline detection did not always work well in Oufti when we considered growing M. xanthus cells in
% time-lapse images, we implemented the option to segment cells using the individual mode in Oufti and then track the
% cells using a modified version of the tracking script from Oufti / MicrobeTracker.
% Note: this also works for cells that do not divide, but for which Oufti
% detected the cell outlines
if timelapse && division 
    disp('cell tracking ...')
    masks = load(mat_path);
    % use trackstack_modified to create new cellList with tracked cells
    new_meshData = trackstack_modified(masks.cellList.meshData, [1 1 1 2], 1); 
    % default values for trackstack: weights = [0.2 1 1 1], cutoff = 1
    % weigths are: 
    % 1) distance perpendicular to the main direction
    % 2) distance along the main direction
    % 3) angle
    % 4) ratio of the areas
    
    % check the results and modify the tracked data accordingly
    prompt = 'Accept tracking of all cells? Y\\N: ';
    checked = input(prompt,'s');
    if(checked == 'N')
        new_meshData = TrackCellsCheckResults(new_meshData, PH_paths);
    end
    
    % convert mat file with tracked cells to new format used in Oufti (with cellIds)
    masks.cellList = oufti_makeNewCellListFromOld(new_meshData);
    
    % save data for tracked cells to new mat file
    mat_name = strsplit(mat_file.name,'.mat');
    save([path mat_name{1} '_tracked.mat'],'-struct','masks');
    mat_path = [path mat_name{1} '_tracked.mat']; % in the following use the new mat file
    disp('done')
end

%% the following functions only apply if fluorescence images exist
if(~fluorescence_images)
    parameters_image_analysis = struct('pixelsize',pixelsize,'deltat',deltat,'timelapse',...
        timelapse,'fluorescence_images',fluorescence_images,'species',species,'division',...
        division,'mat_path',mat_path,'xshift',NaN,'yshift',NaN);
    
    save([path 'info_params'],'parameters_image_analysis')
    return
end

%% Spot detection
% get the right shift of the fluorescence images relative to the phase
% contrast images (using the first frame)
masks = load(mat_path);

% get shifts in x and y direction (x: distance to left edge of image, y:
% distance to top of image)
% Note: only if the data is not shifted, the positions of the spots
% detected with Oufti agree with the segments of the mesh
[xshift,yshift] = gui_shift({fluor_paths{1}, masks.cellList.meshData{1}});

disp('spot detection ...')

% spot detection inside cells detected by Oufti (if ~division, the first
% frame needs to have cell outlines and for those the spots are identified)
AddSpotsToMatFile(mat_path, fluor_paths, xshift, yshift, threshold_value, method, npixel, timelapse, division)

% calculate orientation of spot relative to cell orientation
GetOrientationRelativeToCell(mat_path, xshift, yshift)

% check spot detection results
if timelapse && division
    CheckSpotDetectionTimeLapse(fluor_paths, mat_path, xshift, yshift);
elseif ~timelapse % for snapshots
    gui_spot_detection_indep_frames({fluor_paths, mat_path, xshift, yshift});
end
% NOTE: I do not check spot detection for timelapse && ~division here; this
% is done in the next section!

disp('done')

%% Spot tracking (for cells that do not divide)

if timelapse && ~division
    [align_shifts,fluor_paths_aligned] = AlignImages(5, path, PH_paths, fluor_paths, true);
    [data_spots, boxes] = TrackSpots(mat_path, fluor_paths_aligned, align_shifts, threshold_value, method, npixel, [1,1]);
    data_spots = TrackSpotsCheckResults(fluor_paths_aligned, data_spots, boxes);
end

%% Save mat file to the data directory to store information on parameters used for image analysis

parameters_image_analysis = struct('pixelsize',pixelsize,'deltat',deltat,'timelapse',...
    timelapse,'fluorescence_images',fluorescence_images,'species',species,'division',...
    division,'mat_path',mat_path,'fluor_paths',{fluor_paths},'xshift',xshift,'yshift',yshift,...
    'threshold_value',threshold_value,'method',method,'npixel',npixel);

save([path 'info_params'],'parameters_image_analysis')
