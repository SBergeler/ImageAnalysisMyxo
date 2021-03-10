function data_spots = TrackSpotsCheckResults(fluor_paths_aligned, data_spots, boxes)
% TRACKSPOTSCHECKRESULTS: visually check the results from spot tracking
% (for cells that do not divide) and delete frames / cells that are not
% correctly tracked
%
% Copyright (c) 2021 Silke Bergeler
%
% Input:
% - fluor_paths_aligned: paths to the aligned fluorescence images
% - data_spots: cell array storing the spot data
% - boxes: box sizes
% Output:
% - data_spots: updated cell array with the spot data

ntspots = length(data_spots{1});
ntframes = length(data_spots);

% plot the fluorescence images and the detected spot to visually check the results
stop_frames = zeros(ntspots,1,'int8'); % variable that stores when the spot detection should be stopped (because
% e.g. a different spot is tracked than the desired one)
for spot_ind = 1:ntspots
    disp(['spot ' num2str(spot_ind) ' out of ' num2str(ntspots) ' spots' ])
    stop_frames(spot_ind) = gui_spot_tracking({fluor_paths_aligned, data_spots, spot_ind, boxes});
end

%% delete the data in 'data_spots' that should be excluded from the data analysis (as obtained from the visual inspections)
for spot_ind = 1:ntspots
    nframe = stop_frames(spot_ind);
    data_spots = deleteFollowing(spot_ind, nframe, data_spots, ntframes);
end

end