function [data_spots] = deleteFollowing(spot_ind, nframe, data_spots, nfluoim)
% DELETEFOLLOWING: deletes the entries with index 'spot_ind' in all
% frames following frame 'nframe'
%
% Copyright (c) 2021 Silke Bergeler
%
% Input:
% - spot_ind: index of the spot to be deleted
% - nframe: entries for this spot are deleted for all frames following this
% frame
% - data_spots: data set with the spot infos
% - nfluoim: total number of fluorescence images
% Output:
% - data_spots: updated data set with spot infos

for frame = (nframe+1):nfluoim
    index = arrayfun(@(x) x.spot_ind == spot_ind, data_spots{frame});
    if ~all(index == 0)
        data_spots{frame}(index) = [];
    end
end
end