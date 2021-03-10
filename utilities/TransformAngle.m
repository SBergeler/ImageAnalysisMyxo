function angle_out = TransformAngle(angle_in)
% TRANSFORMANGLE: This function checks if angle_in is in the range [-90,90]
% and if not calculates the corresponding angle in [-90,90]
%
% Copyright (c) 2021 Silke Bergeler
%
% Input:
% - angle_in: input angle, angle in [-180,180] degree range
% Output:
% - angle_out: angle in [-90,90] degree range obtained from adding /
% subtracting 180 degrees

if angle_in > 90
    angle_out = angle_in - 180;
elseif angle_in < -90
    angle_out = angle_in + 180;
else
    angle_out = angle_in;
end
end