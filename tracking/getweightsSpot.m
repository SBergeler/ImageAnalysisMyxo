% This function is a modified version of the function 'getweights' in the
% script 'trackstack.m' from the Oufti software package from the
% Jacobs-Wagner lab, Yale University. The package is published under the
% GNU General Public License as published by the Free Software Foundation
% (version 3).
%
% Copyright (c) 2015, the Jacobs-Wagner lab, Yale University.


function A = getweightsSpot(set1,set2,M) % set: [x, y, area]' for each spot --> 3xn Matrix
n1 = size(set1,2); % n1 should be one since we start with one spot
n2 = size(set2,2);

x1 = repmat(set1(1,:)',1,n2);
x2 = repmat(set2(1,:) ,n1,1);
y1 = repmat(set1(2,:)',1,n2);
y2 = repmat(set2(2,:) ,n1,1);
dx2 = (x2-x1).^2 + (y2-y1).^2;

area1 = repmat(set1(3,:)',1,n2);
area2 = repmat(set2(3,:) ,n1,1);

A = M(1)^2*dx2 + (M(2)*log2(area1./area2)).^2;
end