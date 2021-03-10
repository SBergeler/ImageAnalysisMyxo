function [binary_image_spots, spot_positions, spot_major_axis, spot_minor_axis, ...
    spot_ecc, spot_area, spot_orient, spot_boundary, PixelIdxList, signal_int, rel_signal_int, threshold_spot, background, cell_signal_int] = ...
    SpotDetection2D(image_fluor,background_fluor,this_cell,threshold_value,method,npixel,show_debug_plots)
% SpotDetection2D: Detect spots in cells, which are defined by two
% characteristics: 1) the signal intensity has to be above a certain
% threshold (given by 'method' and 'threshold_value') and 2) spots have to
% consist of at least 'npixel' pixels
%
% Copyright (c) 2021 Silke Bergeler
% We thank Manon Wigbers for her contribution to this script.
%
% INPUT:
% - image_fluor: fluorescence image
% - background_fluor: same image as image_fluor but with 0's inside the
% meshes found by Oufti (note: not all cells might be detected in Oufti!)
% - this_cell: mask from Oufti analysis
% - threshold_value: X for method 1 or 2 
% - method: method for setting threshold (either 1 or 2)
%    1: threshold = mean(cytoplasm) + X*std(cytoplasm)
%    2: threshold = X*mean(cytoplasm)
% - npixel: number of connected pixels to be called a spot
% - show_debug_plots: Boolean (if true, show plots created in this script
% for debugging)

% OUTPUT:
% - binary_image_spots: binary image with ones when the fluorescence signal
% is larger than the threshold value
% - spot_positions: positions of the spots detected 
% - spot_major_axis: major axis of the ellipse that has the same normalized 
% second central moments as the spot region
% - spot_minor_axis: minor axis of the ellipse that has the same normalized 
% second central moments as the spot region
% - spot_ecc: eccentricity of the ellipse that has the same normalized 
% second central moments as the spot region, formula: 
% sqrt(1-(spot_minor_axis./spot_major_axis).^2);
% - spot_area: number of pixels in the spot region
% - spot_orient: angle between the x-axis and the major axis of the ellipse
% that has the same second-moments as the region
% - spot_boundary: row and column coordinates of boundary pixels
% - PixelIdxList: locations of pixels in the spot regions
% - signal_int: summed fluoresence signal intensity of the spots
% detected
% - rel_signal_int: summed fluoresence signal intensity of the spots
% detected divided by the total intensity in the cell
% - threshold_spot: threshold intensity for spots (without background
% correction)
% - background: average intensity outside the cell in a surrounding box
% - cell_signal_int: total intensity in the cell

%% Read in mesh for cell 
mesh = this_cell.mesh;
mesh_min_x = min(min(mesh(:,1),min(mesh(:,3))));
mesh_max_x = max(max(mesh(:,1),max(mesh(:,3))));
mesh_min_y = min(min(mesh(:,2),min(mesh(:,4))));
mesh_max_y = max(max(mesh(:,2),max(mesh(:,4))));

% check if the mesh coordinates are in the image region (due to shifting of
% the meshes it can happen that the mesh is not fully enclosed in the
% fluorescence image region)
if mesh_min_x < 1
    mesh_min_x = 1;
end

if mesh_min_y < 1
    mesh_min_y = 1;
end

if mesh_max_x > size(image_fluor,2)
    mesh_max_x = size(image_fluor,2);
end

if mesh_max_y > size(image_fluor,1)
    mesh_max_y = size(image_fluor,1);
end

% Create binary mask (for cell)
x = double([mesh(:,1);mesh(:,3)])';
y = double([mesh(:,2);mesh(:,4)])';
binary_mask = poly2mask(x,y,size(image_fluor,1),size(image_fluor,2)); % Convert region-of-interest polygon to mask (pixels inside polygon are 1)

%% Substract background from cell signal

% minimal box that includes cell:
box_xmin = uint16(mesh_min_x);
box_xmax = uint16(mesh_max_x);
box_ymin = uint16(mesh_min_y);
box_ymax = uint16(mesh_max_y);

% get fluorescence image with 0's inside the cells that are detected with
% Oufti
background_fluor_small = background_fluor(box_ymin:box_ymax,box_xmin:box_xmax);

% if number of pixels in the background is smaller than 1000, increase box
% size by a factor of 2
if(sum(sum(background_fluor_small>0)) < 1000)
    box_in = zeros(1,4);
    box_in(1) = box_xmin;
    box_in(2) = box_ymin;
    box_in(3) = box_xmax - box_xmin;
    box_in(4) = box_ymax - box_ymin;
    box_out = ChangeBoxSize(box_in, 2); 
    topx = box_out(1);
    topy = box_out(2);
    width = box_out(3);
    height = box_out(4);
    box_xmin = max(1,topx);
    box_xmax = min(size(image_fluor,2),topx+width);
    box_ymin = max(1,topy);
    box_ymax = min(size(image_fluor,1),topy+height);
    background_fluor_small = background_fluor(box_ymin:box_ymax,box_xmin:box_xmax);
end

% define background as the median non-zero signal outside the cell in the box
background = median(background_fluor_small(background_fluor_small>0));

% subtract background from larger fluor image 
fluor_background_corr = image_fluor - background; 

% total fluorescence signal inside the cell (background corrected)
cell_fluor_signal = fluor_background_corr.*uint16(binary_mask);
cell_signal_int = sum(cell_fluor_signal(:));

%% for debugging only:

if(show_debug_plots)
    % smaller image containing only the cell of interest
    zoomed_fluor = image_fluor(box_ymin:box_ymax,box_xmin:box_xmax);
    
    figure()
    subplot(2,3,1)
    imshow(imadjust(zoomed_fluor));
    hold on
    plot(mesh(:,1)-single(box_xmin),mesh(:,2)-single(box_ymin),'Color','cyan','LineWidth',1.2);
    plot(mesh(:,3)-single(box_xmin),mesh(:,4)-single(box_ymin),'Color','cyan','LineWidth',1.2);
    hold off
    title('not BG corrected')
    
    subplot(2,3,2)
    imshow(imadjust(background_fluor_small))
    title('BG')
    
    subplot(2,3,3)
    imshow(imadjust(fluor_background_corr(box_ymin:box_ymax,box_xmin:box_xmax)))
    title('BG corrected')
    
    subplot(2,3,4:6)
    background_vals = background_fluor_small(background_fluor_small>0);
    histogram(background_vals,'Normalization','pdf')
    [N, ~] = histcounts(background_vals,'Normalization','pdf');
    hold on
    %xline(single(median(background_vals)),'LineWidth',1.5,'color','r');
    %xline(single(mean(background_vals)),'LineWidth',1.5,'color','k');
    line([median(background_vals) median(background_vals)], [0 max(N)*1.1],'Color','red','LineStyle','--');
    line([mean(background_vals) mean(background_vals)], [0 max(N)*1.1],'Color','black','LineStyle','--');
    hold off
    xlabel('pixel intensity (BG)')
    ylabel('pdf')
    ylim([0 max(N)*1.1])
    title('Pixel intensities in the BG (red line: median, black line: mean)')
end

%% Find pixels with high intensity
binary_image_spots = logical(false(size(image_fluor,1),size(image_fluor,2))); 

% sort the fluorescence intensities in the cell in ascending order
cell_fluor_signal_diffuse_notBGC = sort(image_fluor(binary_mask));
% use only fraction of 70% smallest intensities (for calculation of
% threshold)
cell_fluor_signal_diffuse_notBGC_small = cell_fluor_signal_diffuse_notBGC(1:round(length(cell_fluor_signal_diffuse_notBGC)*0.7));

% Set threshold spot detection, based on cytoplasmic signal
mean_diffuse_notBGC = mean(double(cell_fluor_signal_diffuse_notBGC_small));
std_diffuse_notBGC = std(double(cell_fluor_signal_diffuse_notBGC_small));

if method == 1
    threshold_spot = ...
        mean_diffuse_notBGC + ...
        threshold_value*std_diffuse_notBGC;
elseif method == 2
    threshold_spot = ...
        threshold_value*mean_diffuse_notBGC;
else
    threshold_spot = ...
        mean_diffuse_notBGC + ...
        threshold_value*std_diffuse_notBGC;
    disp('error: threshold method selection')
end

% fluorescence in the cell (not BG corrected)
cell_fluor_signal_notBGC = image_fluor.*uint16(binary_mask);
binary_image_spots(cell_fluor_signal_notBGC > threshold_spot) = true;

%% for debugging only

if (show_debug_plots)
    figure;
    subplot(2,2,1)
    imshow(imadjust(image_fluor(box_ymin:box_ymax,box_xmin:box_xmax)))
    title('fluorescence signal in cell (not BG corr.)')
    
    subplot(2,2,2)
    imshow(binary_image_spots(box_ymin:box_ymax,box_xmin:box_xmax))
    title('high intensity pixels')
    
    subplot(2,2,[3,4])
    histogram(cell_fluor_signal_notBGC(binary_mask),100,'Normalization','pdf');
    [N,~] = histcounts(cell_fluor_signal_notBGC(binary_mask),100,'Normalization','pdf');
    hold on
    line([threshold_spot threshold_spot], [0 max(N)*1.1],'Color','red','LineStyle','--');
    hold on
    line([mean_diffuse_notBGC mean_diffuse_notBGC], [0 max(N)*1.1],'Color','black','LineStyle','--');
    title('diffuse intensity inside cell, not BGC (red: threshold, black: mean (ignoring 30% highest ints)')
    ylabel('pdf')
    xlabel('fluorescence intensity')
    ylim([0 max(N)*1.1])
end

%% Find spots by connectedness
% save all connected components with connectivity 4
components_1 = bwconncomp(binary_image_spots,4); 
% get statistics of the spots (largest / smallest extension, area, ...) 
stats = regionprops('table',components_1,'Area','Centroid',...
    'MajorAxisLength','MinorAxisLength','PixelIdxList','Eccentricity','Orientation', ...
    'Image');
if size(stats,2) == 1 % to correct for an error that occurred for one-pixel connected components 
    newnames = {'Area','Centroid','MajorAxisLength','MinorAxisLength','PixelIdxList','Eccentricity','Orientation', 'Image'};
    stats = splitvars(stats,'PixelIdxList','NewVariableNames',newnames);
    stats.Area = cell2mat(stats.Area);
end
stats = stats(stats.Area > npixel,:); % choose only the spots with at least 'npixel' pixels
spot_positions = stats.Centroid;
spot_major_axis = stats.MajorAxisLength;
spot_minor_axis = stats.MinorAxisLength;
spot_area = stats.Area;
spot_ecc = stats.Eccentricity;
spot_orient = stats.Orientation;
PixelIdxList = stats.PixelIdxList;
    
% calculate the signal intensity of the spots (sum of all fluorescence
% intensity values (background corrected) in the spot regions)
if isempty(stats.PixelIdxList) % if no spot is found
    signal_int = [];
    rel_signal_int = [];
    spot_boundary = [];    
else
    signal_int = cellfun(@(x) sum(cell_fluor_signal(x)), stats.PixelIdxList);
    rel_signal_int = cellfun(@(x) sum(cell_fluor_signal(x)), stats.PixelIdxList)./cell_signal_int;
    spot_boundary = cellfun(@(x) bwboundaries(x,'noholes'),stats.Image,'UniformOutput',false); % row and column coordinates 
end
end
