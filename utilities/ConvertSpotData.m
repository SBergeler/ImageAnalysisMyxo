function [new_data] = ConvertSpotData(data)
% CONVERTSPOTDATA: this function converts the spot data in the form of a
% struct with one entry per spot to a struct with arrays as entries
%
% Copyright (c) 2021 Silke Bergeler
%
% Input:
% - data: spot data to be converted
% Output:
% - new_data: spot data in the new format

fields = fieldnames(data);
new_data = struct;
for i = 1:numel(fields)
    new_data.(fields{i}) = [data.(fields{i})];
end
end