function inbox = InBox(box,x,y)
% INBOX: Determine whether the point [x,y] is in the box [topx, topy, width, height]
%
% Copyright (c) 2021 Silke Bergeler
%
% Input:
% - box: box of interest
% - x / y: coordinates of the point to test if they are in the box
% Output: 
% - inbox: Boolean, true if the point is in the box

topx = box(1);
topy = box(2);
width = box(3);
height = box(4);

if topx <= x && x <= topx + width && topy <= y && y <= topy + height
    inbox = true;
else 
    inbox = false;
end
end

