function [data_spots, boxes] = TrackSpots(mat_path, fluor_paths_aligned, align_shifts, threshold_value, method, npixel, M)
% TRACKSPOTS: track spots in time lapse images with cell outlines only for
% the first frame
% Note: Cells that do not divide become very long and curvy quickly such that
%   we could not get satisfactory cell outlines using Oufti --> in this case we use
%   the cell outlines in the first frames and track the spots identified in
%   these cells without using the cell outlines in the following frames
%
% Copyright (c) 2021 Silke Bergeler
%
% Input:
% - mat_path: path to mat file
% - fluor_paths_aligned: paths to aligned fluorescence images
% - align_shifts: array with the shifts of the fluorescence images to
%     align them
% - threshold_value: threshold value chosen for spot detection
% - method: method chosen for spot detection
% - npixel: number of pixels for spot detection
% - M: weights used for the distance measure between spots (used to assign spots in different frames) 
%
% Output:
% - data_spots: cell array storing the spot data
% - boxes: box sizes

masks = load(mat_path);
nfluoim = size(fluor_paths_aligned,1); % total number of fluorescence images
ntframes = size(masks.cellList.meshData,2); % total number of frames analyzed by Oufti
if ntframes == 0
    print('No cell outlines are detected with Oufti in the first frame!')
end

%% Write spot detection results for first frame to a cell array
data_spots = cell(nfluoim,1); % cell array storing the spot data for each frame
meshes = masks.cellList.meshData{1,1}; % only consider spots that are detected in the first frame

% write the data obtained from spot detection of the first frame to
% data_spots
data_spots{1} = struct('spot_ind',{},'box',{},'background',{},'x',{},'y',{},'spot_major_axis',{}, ...
    'spot_minor_axis',{},'spot_ecc',{},'spot_area',{},'spot_orient',{},...
    'spot_boundary',{},'signal_int',{},'threshold_spot',{});

% index of the spots that are considered 
spot_ind = 1;
% width and height of the box considered for spot detection (the position of the box is based on the
% spot detected in the previous frame)
box_width = 100;
box_height = 100;
image_fluor_first_frame = imread(fluor_paths_aligned{1});
% for all cells detected in the first frame
for ncell = 1:length(meshes) 
    % number of spots detected in the cell
    no_spots = length(meshes{ncell}.spots_matlab);
    if no_spots == 1 && ~any(structfun(@isempty, meshes{ncell}.spots_matlab)) % only consider cells with exactly one spot detected
        data_spots{1}(spot_ind).spot_ind = spot_ind;
        % positions of the spot in the fluorescence image (not aligned) as
        % obtained with SpotDetection2D
        x_old = meshes{ncell}.spots_matlab.x; 
        y_old = meshes{ncell}.spots_matlab.y;
        box = [round(x_old-box_width/2), round(y_old-box_height/2), box_width, box_height] + [align_shifts(1,1), align_shifts(1,2),0,0];
        background = meshes{ncell}.background;
        data_spots{1}(spot_ind).box = box;
        data_spots{1}(spot_ind).background = background;
        [~, spot_positions, spot_major_axis, ...
            spot_minor_axis, spot_ecc, spot_area, spot_orient, spot_boundary, signal_int, threshold_spot] ...
            = SpotDetection2DArea(image_fluor_first_frame,box,background,threshold_value,method,npixel);
        % Since we now consider a larger region for spot detection, it
        % might be that more than one spot is detected. We chose the one
        % that is most similar to the one detected with SpotDetection2D.
        % However, since we consider cells with one spot detected, at least one spot
        % (the one detected previously) should be detected again
        if length(spot_minor_axis) ~= 1 
            set1 = [x_old + align_shifts(1,1); y_old + align_shifts(1,2); meshes{ncell}.spots_matlab.spot_area];
            set2 = [spot_positions(:,1)'; spot_positions(:,2)'; spot_area'];
            if(isempty(set2))
                disp('No spot is detected in the first frame (should not be the case!)')
            end
            A = getweightsSpot(set1,set2,M);
            [~,ind] = min(A);
        else
            ind = 1;
        end
        data_spots{1}(spot_ind).x = spot_positions(ind,1);
        data_spots{1}(spot_ind).y = spot_positions(ind,2);
        data_spots{1}(spot_ind).spot_major_axis = spot_major_axis(ind);
        data_spots{1}(spot_ind).spot_minor_axis = spot_minor_axis(ind);
        data_spots{1}(spot_ind).spot_ecc = spot_ecc(ind);
        data_spots{1}(spot_ind).spot_area = spot_area(ind);
        data_spots{1}(spot_ind).spot_orient = spot_orient(ind);
        data_spots{1}(spot_ind).spot_boundary = spot_boundary(ind);
        data_spots{1}(spot_ind).signal_int = signal_int(ind);
        data_spots{1}(spot_ind).threshold_spot = threshold_spot;
        % Alternative: use the values obtained from spot detection inside
        % the cell only
        %data_spots{1}(cell_ind).x = InRangeX(meshes{1,ncell}.spots_matlab.x + align_shifts(1,1));
        %data_spots{1}(cell_ind).y = InRangeY(meshes{1,ncell}.spots_matlab.y + align_shifts(1,2));
        %data_spots{1}(cell_ind).spot_major_axis = meshes{1,ncell}.spots_matlab.spot_major_axis;
        %data_spots{1}(cell_ind).spot_minor_axis = meshes{1,ncell}.spots_matlab.spot_minor_axis;
        %data_spots{1}(cell_ind).spot_ecc = meshes{1,ncell}.spots_matlab.spot_ecc;
        %data_spots{1}(cell_ind).spot_area = meshes{1,ncell}.spots_matlab.spot_area;
        %data_spots{1}(cell_ind).spot_orient = meshes{1,ncell}.spots_matlab.spot_orient;
        %data_spots{1}(cell_ind).spot_boundary = meshes{1,ncell}.spots_matlab.spot_boundary; % spot boundary is not needed, or?
        %data_spots{1}(cell_ind).threshold_spot = meshes{1,ncell}.threshold_spot;
        spot_ind = spot_ind + 1;
    end
end
ntspots = length(data_spots{1}); % number of spots in first frame considered further

%% for the following frames use 'SpotDetection2DArea' to detect spots
% search for spots in the defined regions 
disp('tracking the spots ...')
% cell array storing the binary images for each frame
binary_image_cell_array = cell(nfluoim,1); % first frame is missing!
% binary vector that stores if the spot is still tracked or not (if it is lost
% once, tracking is not continued)
cell_tracking = true(ntspots, 1); 
for nframe = 2:nfluoim
    disp(['frame: ', num2str(nframe)])
    image_fluor = imread(fluor_paths_aligned{nframe});
    binary_image = logical(false(size(image_fluor,1),size(image_fluor,2)));
    data = struct('spot_ind',{},'box',{},'background',{},'x',{},'y',{},'spot_major_axis',{}, ...
        'spot_minor_axis',{},'spot_ecc',{},'spot_area',{},'spot_orient',{},...
        'spot_boundary',{},'signal_int',{},'threshold_spot',{});
    % row index in data_spots
    row_ind = 1;
    for spot_ind = 1:ntspots
        if cell_tracking(spot_ind)
            % use the function SpotDetection2DArea to detect the spots in the
            % area around the spot (including the spot)
            % index in data_spots of the cell considered 
            index = arrayfun(@(x) x.spot_ind == spot_ind, data_spots{nframe-1});
            box = [round(data_spots{nframe-1}(index).x-box_width/2), round(data_spots{nframe-1}(index).y-box_height/2), box_width, box_height];
            background = data_spots{1}(index).background;
            [binary_image_spots, spot_positions, spot_major_axis, ...
                spot_minor_axis, spot_ecc, spot_area, spot_orient, spot_boundary, signal_int, threshold_spot] ...
                = SpotDetection2DArea(image_fluor,box,background,threshold_value,method,npixel);
            % chose the spot that is the most similar one to the spot in the
            % previous frame
            if length(spot_area) ~= 1
                if(isempty(spot_area))
                    % no spot is detected --> stop considering this cell in the following frames
                    cell_tracking(spot_ind) = false;
                    continue;
                end
                set1 = [data_spots{nframe-1}(index).x; data_spots{nframe-1}(index).y; data_spots{nframe-1}(index).spot_area];
                set2 = [spot_positions(:,1)'; spot_positions(:,2)'; spot_area'];
                A = getweightsSpot(set1,set2,M);
                [~,ind] = min(A); % TODO: use a threshold to avoid bad assignment!
            else
                ind = 1;
            end
            data(row_ind).box = box;
            data(row_ind).background = background;
            data(row_ind).spot_ind = spot_ind; 
            data(row_ind).x = spot_positions(ind,1);
            data(row_ind).y = spot_positions(ind,2);
            data(row_ind).spot_major_axis = spot_major_axis(ind);
            data(row_ind).spot_minor_axis = spot_minor_axis(ind);
            data(row_ind).spot_ecc = spot_ecc(ind);
            data(row_ind).spot_area = spot_area(ind);
            data(row_ind).spot_orient = spot_orient(ind);
            data(row_ind).spot_boundary = spot_boundary(ind);
            data(row_ind).signal_int = signal_int(ind);
            data(row_ind).threshold_spot = threshold_spot;
            binary_image = binary_image + binary_image_spots; % add the cell mask to the binary image
            row_ind = row_ind + 1;
        end
    end
    data_spots{nframe} = data;
    binary_image_cell_array{nframe} = binary_image;
end
disp('done')

%% determine the box size used for plotting the images (such that all boxes for spot detection are included)
boxes = zeros(ntspots,4);
for spot_ind = 1:ntspots
    box_x_min = [];
    box_x_max = [];
    box_y_min = [];
    box_y_max = [];
    for nframe = 1:nfluoim 
        if ismember(spot_ind, arrayfun(@(x) x.spot_ind, data_spots{nframe}))
            index = arrayfun(@(x) x.spot_ind == spot_ind, data_spots{nframe});
            box_x_min = [box_x_min data_spots{nframe}(index).box(1)];
            box_x_max = [box_x_max data_spots{nframe}(index).box(1) + data_spots{nframe}(index).box(3)];
            box_y_min = [box_y_min data_spots{nframe}(index).box(2)];
            box_y_max = [box_y_max data_spots{nframe}(index).box(2) + data_spots{nframe}(index).box(4)];
        else
            break;
        end       
    end
    boxes(spot_ind,1) = min(box_x_min);
    boxes(spot_ind,2) = max(box_x_max);
    boxes(spot_ind,3) = min(box_y_min);
    boxes(spot_ind,4) = max(box_y_max);
end

end