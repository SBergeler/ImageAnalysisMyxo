function [ ] = AddSpotsToMatFile(mat_path, fluor_paths, xshift, yshift, threshold_value, method, npixel, timelapse, division)
% ADDSPOTSTOMATFILE: This function determines spots inside cells using
% SpotDetection2D and adds the spots as spots_matlab to the mat file
%
% Copyright (c) 2021 Silke Bergeler
%
% Input:
% - mat_path: path to mat file
% - fluor_paths: paths to fluorescence images
% - xshift: how many pixels the meshes have to be shifted in x direction to
% match the fluorescence signals
% - yshift: how many pixels the meshes have to be shifted in y direction to
% match the fluorescence signals
% - threshold_value: X for method 1 or 2 
% - method: method for setting threshold
%   1: threshold = mean(cytoplasm) + X*std(cytoplasm)
%   2: threshold = X*mean(cytoplasm));
% - npixel: number of connected pixels to be called a spot
% - timelapse: Boolean, true if time lapse images are considered
% - division: Boolean, true if cells divide

%% Spot Detection
masks = load(mat_path);
nfluoim = size(fluor_paths,1); % total number of fluorescence images
ntframes = size(masks.cellList.meshData,2); % total number of frames analyzed by Oufti
% Note: typically the number of fluorescence images and the number of
% frames analyzed by Oufti are the same. But they might differ, if we are
% only interested in a subset of cells

image_fluor_cell_array = cell(nfluoim,1); % cell array storing the fluorescence image for each frame
binary_image_cell_array = cell(nfluoim,1); % cell array storing the binary images for each frame
mesh_cell_array = cell(nfluoim,1); % cell array storing the mesh informations for all cells in each frame
data_spots = cell(nfluoim,1); % cell array storing the spot position data for each frame

show_debug_plots = true;
count = 0;

for nframe = 1:min(nfluoim,ntframes)
    image_fluor = imread(fluor_paths{nframe});
    % get an image with the intensities inside the cells set to zero
    background_fluor = getBackgroundImage(image_fluor,masks.cellList.meshData{nframe}, xshift, yshift);
    % create binary image of the same size as the fluorescence image
    binary_image = logical(false(size(image_fluor,1),size(image_fluor,2)));
    data = struct('ncell',{},'spot_positions',{},'spot_major_axis',{}, ...
        'spot_minor_axis',{},'spot_ecc',{},'spot_area',{},'spot_orient',{},...
        'spot_boundary',{},'signal_int',{},'rel_signal_int',{},'threshold_spot',{},'background',{},'cell_signal_int',{});
    ntcells = length(masks.cellList.meshData{1,nframe}); % total number of cells detected in this frame
    meshs = cell(ntcells,1); % store the meshes of the cells in this frame (shift corrected)
    for ncell = 1:ntcells
        count = count + 1;
        if count == 1
            show_debug_plots = false;
        end
        this_cell = masks.cellList.meshData{1,nframe}{1,ncell};
        % change the meshes according to the shift 
        this_cell.mesh(:,1) = this_cell.mesh(:,1)+xshift;
        this_cell.mesh(:,2) = this_cell.mesh(:,2)+yshift;
        this_cell.mesh(:,3) = this_cell.mesh(:,3)+xshift;
        this_cell.mesh(:,4) = this_cell.mesh(:,4)+yshift;
        % use the function SpotDetection2D to detect the spots in this cell
        [binary_image_spots, spot_positions, spot_major_axis, ...
            spot_minor_axis, spot_ecc, spot_area, spot_orient, spot_boundary, PixelIdxList, ...
            signal_int, rel_signal_int, threshold_spot, background, cell_signal_int] ...
            = SpotDetection2D(image_fluor,background_fluor,this_cell,threshold_value,method,npixel,show_debug_plots);
        if length(spot_major_axis) > 1 % if more than one spot is detected
            cell_id = masks.cellList.cellId{nframe}(ncell);
            if timelapse && division && nframe > 1
                ind = masks.cellList.cellId{nframe-1} == cell_id;
                % if the cell was not identified in the previous frame
                if ~any(ind)
                    PlotCellAndSpots(this_cell,image_fluor,PixelIdxList,nframe,cell_id,spot_positions,spot_area,signal_int)
                else
                    previous_cell = masks.cellList.meshData{nframe-1}{ind};
                    % change the meshes according to the shift
                    previous_cell.mesh(:,1) = previous_cell.mesh(:,1)+xshift;
                    previous_cell.mesh(:,2) = previous_cell.mesh(:,2)+yshift;
                    previous_cell.mesh(:,3) = previous_cell.mesh(:,3)+xshift;
                    previous_cell.mesh(:,4) = previous_cell.mesh(:,4)+yshift;
                    previous_image_fluor = imread(fluor_paths{nframe-1});
                    previous_spot_positions = data_spots{nframe-1}(ind).spot_positions;
                    previous_spot_area = data_spots{nframe-1}(ind).spot_area;
                    previous_signal_int = data_spots{nframe-1}(ind).signal_int;
                    previous_PixelIdxList = data_spots{nframe-1}(ind).PixelIdxList;
                    PlotCellAndSpotsTL(this_cell,previous_cell,image_fluor,PixelIdxList,previous_image_fluor,previous_PixelIdxList,nframe,cell_id,...
                        spot_positions,spot_area,signal_int,previous_spot_positions,previous_spot_area,previous_signal_int)
                end
            else
                PlotCellAndSpots(this_cell,image_fluor,PixelIdxList,nframe,cell_id,spot_positions,spot_area,signal_int)
            end
            prompt = 'Which spot ID to use? [1,2,...; press return for all; press 0 for no spots] '; % Note: at the moment it is not possible to choose more than 1 spot, but not all
            spot_id = input(prompt);
            if ~isempty(spot_id) && (spot_id > length(spot_major_axis) || spot_id < 0)
                disp('invalid input')
                figure(gcf)
                prompt = 'Which spot ID to use? [1,2,...; press return for all; press 0 for no spots] '; % Note: at the moment it is not possible to choose more than 1 spot, but not all
                spot_id = input(prompt);
            end
                
            close(gcf)
            if ~isempty(spot_id) 
                if spot_id ~= 0
                    spot_positions = spot_positions(spot_id,:);
                    spot_major_axis = spot_major_axis(spot_id);
                    spot_minor_axis = spot_minor_axis(spot_id);
                    spot_ecc = spot_ecc(spot_id);
                    spot_area = spot_area(spot_id);
                    spot_orient = spot_orient(spot_id);
                    spot_boundary = spot_boundary{spot_id};
                    PixelIdxList = PixelIdxList{spot_id};
                    signal_int = signal_int(spot_id);
                    rel_signal_int = rel_signal_int(spot_id);
                else 
                    spot_positions = double.empty(0,2);
                    spot_major_axis = double.empty(0,1);
                    spot_minor_axis = double.empty(0,1);
                    spot_ecc = double.empty(0,1);
                    spot_area = double.empty(0,1);
                    spot_orient = double.empty(0,1);
                    spot_boundary = [];
                    PixelIdxList = [];
                    signal_int = [];
                    rel_signal_int = [];
                end
            end
        end
        data(ncell).ncell = ncell;
        data(ncell).spot_positions = spot_positions;
        data(ncell).spot_major_axis = spot_major_axis;
        data(ncell).spot_minor_axis = spot_minor_axis;
        data(ncell).spot_ecc = spot_ecc;
        data(ncell).spot_area = spot_area;
        data(ncell).spot_orient = spot_orient;
        data(ncell).spot_boundary = spot_boundary;
        if ~iscell(PixelIdxList) % sometimes the PixelIdxList is a double array
            PixelIdxList_tmp = cell(1,1);
            PixelIdxList_tmp{1} = PixelIdxList;
            PixelIdxList = PixelIdxList_tmp;
        end
        data(ncell).PixelIdxList = PixelIdxList;
        data(ncell).signal_int = signal_int;
        data(ncell).rel_signal_int = rel_signal_int;
        data(ncell).threshold_spot = threshold_spot;
        data(ncell).background = background;
        data(ncell).cell_signal_int = cell_signal_int;
        binary_image = binary_image + binary_image_spots; % add the cell mask to the binary image
        meshs{ncell,1} = this_cell.mesh; % save the mesh (shift corrected) of this cell
    end
    image_fluor_cell_array{nframe,1} = image_fluor;
    binary_image_cell_array{nframe,1} = binary_image;
    mesh_cell_array{nframe,1} = meshs; 
    data_spots{nframe,1} = data;
end
                      
%% Add the spot data to the masks as a struct
masks2 = load(mat_path); % load the data in the mat file again and add the spot positions from Matlab (this is because the mesh data should not be changed)
for nframe = 1:size(masks2.cellList.meshData,2) % consider all frames analyzed by Oufti
    if nframe <= length(data_spots)
        data_nframe = data_spots{nframe,1};
    else
        data_nframe = zeros(0,3); % if we didn't search for spots in all frames analyzed by Oufti
    end
    for ncell = 1:size(masks2.cellList.meshData{1,nframe},2)
        data_spot_in_cell = data_nframe(ncell); 
        % save the threshold value used to detect spots, the background fluorescence and the total cell intensity to the mesh data
        % for the cell
        masks2.cellList.meshData{nframe}{ncell}.threshold_spot = data_spot_in_cell.threshold_spot;
        masks2.cellList.meshData{nframe}{ncell}.background = data_spot_in_cell.background;
        masks2.cellList.meshData{nframe}{ncell}.cell_signal_int = data_spot_in_cell.cell_signal_int;
        % save the data for the spots to spots_matlab
        spot_positions = data_spot_in_cell.spot_positions;
        spot_major_axis = num2cell(data_spot_in_cell.spot_major_axis');
        spot_minor_axis = num2cell(data_spot_in_cell.spot_minor_axis');
        spot_ecc = num2cell(data_spot_in_cell.spot_ecc');
        spot_area = num2cell(data_spot_in_cell.spot_area');
        spot_orient = num2cell(data_spot_in_cell.spot_orient');
        spot_boundary = num2cell(data_spot_in_cell.spot_boundary');
        signal_int = num2cell(data_spot_in_cell.signal_int');
        rel_signal_int = num2cell(data_spot_in_cell.rel_signal_int');
        % calculate the coordinates relative to the cell for each spot
        ntspots = size(rel_signal_int, 2); % total number of spots found in this cell
        if ntspots == 0
            masks2.cellList.meshData{nframe}{ncell}.spots_matlab = ...
                struct('l',[],'d',[],'x',[],'y',[],'positions',[],...
                'spot_major_axis',[],'spot_minor_axis',[], ...
                'spot_ecc',[],'spot_area',[],'spot_orient',[],...
                'spot_boundary',[],'signal_int',[],'rel_signal_int',[]);
        else
            l = zeros(1,ntspots);
            d = zeros(1,ntspots);
            x = zeros(1,ntspots);
            y = zeros(1,ntspots);
            positions = zeros(1,ntspots);
            ind = 1;
            err = false; % boolean which indicates if something goes wrong in ChangeCoordinates and if yes, the spots are not stored
            for nspot = 1:ntspots
                % Use the function 'ChangeCoordinates' to get the
                % coordinates of the spot relative to the cell; if a spot is
                % detected at the boundary of the cell, the spot might not
                % be inside any of the mesh segments, which leads to l = d
                % = position = -1
                [lnew,dnew,position] = ChangeCoordinates(masks2.cellList.meshData{nframe}{ncell}, ...
                    [spot_positions(nspot,1),spot_positions(nspot,2)], xshift, yshift);
                if position ~= -1
                    l(ind) = lnew;
                    d(ind) = dnew;
                    x(ind) = spot_positions(nspot,1);
                    y(ind) = spot_positions(nspot,2);
                    positions(ind) = position;
                    ind = ind + 1;
                else
                    err = true;
                end
            end
            if err % Note: we neglect all spots also if an error occured only for one spot!
                masks2.cellList.meshData{nframe}{ncell}.spots_matlab = ...
                struct('l',[],'d',[],'x',[],'y',[],'positions',[],...
                'spot_major_axis',[],'spot_minor_axis',[], ...
                'spot_ecc',[],'spot_area',[],'spot_orient',[],...
                'spot_boundary',[],'signal_int',[],'rel_signal_int',[]);
            else              
                l = num2cell(l);
                d = num2cell(d);
                x = num2cell(x);
                y = num2cell(y);
                positions = num2cell(positions);
                masks2.cellList.meshData{nframe}{ncell}.spots_matlab = ...
                    struct('l',l,'d',d,'x',x,'y',y,'positions',positions,...
                    'spot_major_axis',spot_major_axis,'spot_minor_axis',spot_minor_axis, ...
                    'spot_ecc',spot_ecc,'spot_area',spot_area,'spot_orient',spot_orient,...
                    'spot_boundary',spot_boundary,'signal_int',signal_int,...
                    'rel_signal_int',rel_signal_int);
            end
        end
    end
end
% add the positions of the spots determined with the Matlab script SpotDetection2D to the mat
% file (the mesh is not shifted!)
save(mat_path,'-struct','masks2');

end
