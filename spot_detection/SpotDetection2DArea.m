function [binary_image_spots, spot_positions, spot_major_axis, spot_minor_axis, ...
    spot_ecc, spot_area, spot_orient, spot_boundary, signal_int,threshold_spot] = ...
    SpotDetection2DArea(image_fluor,this_area,background,threshold_value,method,npixel)
% SPOTDETECTION2DAREA: Detect spots in the area 'this_area', which are defined by two
% characteristics: 1) the signal intensity has to be above a certain
% threshold (given by 'method' and 'threshold_value') and 2) spots have to
% consist of at least 'npixel' pixels
%
% Copyright (c) 2021 Silke Bergeler
% We thank Manon Wigbers for her contribution to this script.
%
% INPUT:
% - image_fluor: fluorescence image
% - this_area: area of the fluorescence image searched for spots
% - background: background signal which is subtracted
% - threshold_value: X for method 1 or 2 
% - method: method for setting threshold (either 1 or 2)
%    1: threshold = mean(cytoplasm) + X*std(cytoplasm)
%    2: threshold = X*mean(cytoplasm)
% - npixel: number of connected pixels to be called a spot
% Note: this method only works for a signal that is mainly in a spot
% and not diffuse inside the cell!

% OUTPUT:
% - binary_image_spots: binary image with ones where the fluorescence signal
% is larger than the threshold value
% - spot_positions: positions of the spots detected 
% - spot_major_axis: major axis of the ellipse that has the same normalized 
% second central moments as the spot region
% - spot_minor_axis: minor axis of the ellipse that has the same normalized 
% second central moments as the spot region
% - spot_ecc: eccentricity of the ellipse that has the same normalized 
% second central moments as the spot region, formula: 
% sqrt(1-(spot_minor_axis./spot_major_axis).^2);
% - spot_area: area of the spot
% - spot_orient: orientation of the spot on the image
% - spot_boundary: boundary of the spots
% - signal_int: total fluorescence signal of the spot
% - threshold_spot: threshold intensity for spots

%% Find pixels with high intensity

% background correction
image_fluor = image_fluor - background; % subtract background from fluor image (values are rounded to next int and negative values are set to 0)

binary_image_spots = logical(false(size(image_fluor,1),size(image_fluor,2))); 

% this_area has the format: top x, top y, width, height
topx = this_area(1);
topy = this_area(2);
width = this_area(3);
height = this_area(4);
% due to the shift of the images when aligning, the box might be outside
% the valid area -> crop the box
zoomed_fluor = image_fluor(max(1,topy):min(size(image_fluor,1),topy+height),max(1,topx):min(size(image_fluor,2),topx+width)); 

% binary_mask is 1 where the area of interest is and zero otherwise
% binary_mask = logical(false(size(image_fluor,1),size(image_fluor,2))); 
% binary_mask(max(1,topy):min(size(image_fluor,1),topy+height),max(1,topx):min(size(image_fluor,2),topx+width)) = true;

% % for debugging only
% imshow(imadjust(zoomed_fluor))
% hold on
% mesh = masks.cellList.meshData{1,1}{1,2}.mesh;
% plot(mesh(:,1)-topx,mesh(:,2)-topy,'Color','cyan','LineWidth',1.5)
% hold on
% plot(mesh(:,3)-topx,mesh(:,4)-topy,'Color','cyan','LineWidth',1.5)

% Determine diffuse signal in the area
fluor_signal_diffuse = sort(zoomed_fluor(:)); % sort the intensities in ascending order
 
% To ignore the spots when calculating the mean and std of the diffuse intensity signal in the area, we cut off the 5% highest intensities 
fluor_signal_diffuse_small = fluor_signal_diffuse(1:round(length(fluor_signal_diffuse)*0.95)); 

% Set threshold spot detection, based on diffuse signal in the area
mean_diffuse = mean(double(fluor_signal_diffuse_small));
if method == 1
    threshold_spot = ...
        mean_diffuse + ...
        threshold_value*std(double(fluor_signal_diffuse_small));
elseif method == 2
    threshold_spot = ...
        threshold_value*mean_diffuse;
else
    threshold_spot = ...
        mean_diffuse + ...
        threshold_value*std(double(fluor_signal_diffuse_small));
    disp('error: threshold method selection')
end
 
binary_image_small = logical(false(size(zoomed_fluor,1),size(zoomed_fluor,2)));
binary_image_small(zoomed_fluor > threshold_spot) = true;
binary_image_spots(max(1,topy):min(size(image_fluor,1),topy+height),max(1,topx):min(size(image_fluor,2),topx+width)) = binary_image_small;

% % for debugging only
% % plot the intensity distribution inside the cell with the threshold to
% % detect spots
% figure;
% histogram(zoomed_fluor, 'Normalization','probability');
% hold on 
% line([threshold_spot threshold_spot], [0 0.1],'Color','red','LineStyle','--');
% hold on 
% line([mean_diffuse mean_diffuse], [0 0.1],'Color','black','LineStyle','--');
% figure;
% imshow(binary_image_spots(topy:(topy+height),topx:(topx+width)));

%% Find spots by connectedness
components_1 = bwconncomp(binary_image_spots,4); % save all connected components with connectivity 4
% get statistics of the spots (largest / smallest extension, area, ...) 
stats = regionprops('table',components_1,'Area','Centroid',...
    'MajorAxisLength','MinorAxisLength','PixelIdxList','Eccentricity', 'Orientation', ...
    'Image');
stats = stats(stats.Area > npixel,:); % choose only the spots with at least 'npixel' pixels
spot_positions = stats.Centroid;
spot_major_axis = stats.MajorAxisLength;
spot_minor_axis = stats.MinorAxisLength;
spot_area = stats.Area;
spot_ecc = stats.Eccentricity;
spot_orient = stats.Orientation;
    
if isempty(stats.PixelIdxList) % if no spot is found
    spot_boundary = []; 
    signal_int = [];
else
    spot_boundary = cellfun(@(x) bwboundaries(x,'noholes'),stats.Image,'UniformOutput',false); % row and column coordinates 
    signal_int = cellfun(@(x) sum(image_fluor(x)), stats.PixelIdxList);
end

%% For debugging only
% figure(2)
% subplot(1,3,1)
% imshow(binary_image_spots(topy:(topy+height),topx:(topx+width)))
% title('spots')
% 
% subplot(1,3,2)
% imshow(imadjust(zoomed_fluor))
% title('fluorescence intensity')
% 
% subplot(1,3,3)
% imshow(binary_image_spots(topy:(topy+height),topx:(topx+width)))
% hold on
% plot(spot_positions(:,1)-topx,spot_positions(:,2)-topy, 'xr')
% title('high intensity pxs')
% 
% figure(3)
% int_threshold = fluor_signal_diffuse(round(length(fluor_signal_diffuse)*0.7));
% histogram(fluor_signal_diffuse,50)
% hold on
% line([int_threshold int_threshold], get(gca, 'ylim'),'Color','red');
% title('Fluorescence intensity inside the cell')
end
