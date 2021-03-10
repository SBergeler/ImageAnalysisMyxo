function [ l, d, position ] = ChangeCoordinates(cell, spotdata, xshift, yshift)
% CHANGECOORDINATES: Calculate the l and d coordinates (along center line and
% vertical to it) from the x and y positions
%
% Copyright (c) 2021 Silke Bergeler
%
% Input: 
% - cell: data for one cell from Oufti
% - spotdata: spot positions to be transformed in the new coordinates [x y]
% - xshift: shift of the mesh in x direction w.r.t the fluorescence image
% - yshift: shift of the mesh in y direction w.r.t the fluorescence image
% OUTPUT: 
% - l: position of spot in cell coordinates (along center line)
% - d: position of spot in cell coordinates (transversal to center line)
% - position: segment number of the mesh in which the spot is located

%% get the segment of the mesh where the spot is located
% initialize l,d and position
l = -1;
d = -1;
position = -1;

mesh = cell.mesh;
% change the meshes according to the shifts
mesh(:,[1 3]) = mesh(:,[1 3]) + xshift;
mesh(:,[2 4]) = mesh(:,[2 4]) + yshift;
% find the polygon element where the spot is located
for p = 1:(size(mesh,1)-1)
    polyx = [mesh(p,1) mesh(p+1,1) mesh(p+1,3) mesh(p,3) mesh(p,1)];
    polyy = [mesh(p,2) mesh(p+1,2) mesh(p+1,4) mesh(p,4) mesh(p,2)];
    [in,on] = inpolygon(spotdata(1),spotdata(2),polyx,polyy);
    if(in || on)
        position = p;
        break;
    end
end

% if position did not change, the spot is not inside any of the mesh
% segments --> end this function
if position == -1
    return
end

% calculate the distance d and l:
% the centerline is given as:
centerx = mean([mesh(:,1) mesh(:,3)],2);
centery = mean([mesh(:,2) mesh(:,4)],2);
% vector from the centerline point (left) to the spot
diffvec = [spotdata(1)-centerx(position) spotdata(2)-centery(position)];
% vector on the centerline of the segment in which the spot
% is located
centerlinevec = [centerx(position+1)-centerx(position) centery(position+1)-centery(position)];
% vector more or less perpendicular to the center line (left edge of
% segment)
dvec = [mesh(position,1)-mesh(position,3) mesh(position,2)-mesh(position,4)];
% dvec*centerlinevec'; % the scalar product of the 2 vectors is typically nonzero -->
% not exactly perpendicular
% add up all length of the center line elements before the segment in which
% the spot is located
l = 0;
for j = 1:(position-1)
    l = l + norm([centerx(j+1)-centerx(j) centery(j+1)-centery(j)]);
end
% Method 1: project the vector diffvec to the centerlinevec and dvec
% l1 = l + diffvec*centerlinevec'/norm(centerlinevec);
% d1 = diffvec*dvec'/norm(dvec);
% Method 2 (preferred): determine alpha and beta such that
% alpha*centerlinevec+beta*dvec = diffvec
coeffmat = [centerlinevec' dvec']\diffvec'; % coeffmat = [alpha beta]'; A\b is a faster version of inv(A)*b
l2 = l + norm(coeffmat(1)*centerlinevec);
d2 = sign(coeffmat(2))*norm(coeffmat(2)*dvec); 
% test if method2 is correct:
% coeffmat(1)*centerlinevec + coeffmat(2)*dvec == diffvec
% return l and d
l = l2;
d = d2;

end