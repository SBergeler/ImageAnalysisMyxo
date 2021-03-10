function y = Mapping(x)
% MAPPING: maps values in [0,1] to [0,0.5] by mirroring the values at 0.5
%
% Copyright (c) 2021 Silke Bergeler
%
% Input: 
% - x: value to be transformed, in [0,1]
% Output:
% - y: transformed value, in [0,0.5]

if x <= 0.5
    y = x;
else
    y = 1-x;
end
end