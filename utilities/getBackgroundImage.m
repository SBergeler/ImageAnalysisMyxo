function background_fluor = getBackgroundImage(image_fluor,cells, xshift, yshift)
% GETBACKGROUNDIMAGE: set the intensities inside the cells to zero (output
% is used for background correction)
%
% Copyright (c) 2021 Silke Bergeler
%
% Input:
% - image_fluor: fluorescence image
% - cells: info about cells in this frame from Oufti (including mesh info)
% - xshift: how many pixels the meshes have to be shifted in x direction to
% match the fluorescence signals
% - yshift: how many pixels the meshes have to be shifted in y direction to
% match the fluorescence signals
% Output:
% - background_fluor: fluorescence image with the intensities inside the
% detected cells set to zero

background_fluor = image_fluor;

for i = 1:length(cells)
    mesh = cells{i}.mesh;
    
    % change the meshes according to the shifts
    mesh(:,1) = mesh(:,1) + xshift;
    mesh(:,2) = mesh(:,2) + yshift;
    mesh(:,3) = mesh(:,3) + xshift;
    mesh(:,4) = mesh(:,4) + yshift;
    
    % Create binary mask (for cell)
    x = double([mesh(:,1);mesh(:,3)])';
    y = double([mesh(:,2);mesh(:,4)])';
    
    % Convert region-of-interest polygon to mask (pixels inside polygon are 1)
    binary_mask = poly2mask(x,y,size(image_fluor,1),size(image_fluor,2)); 
    
    % set the fluorescence values in the cell meshes to zero
    background_fluor = background_fluor.*uint16(imcomplement(binary_mask));
end