function [align_shifts, fluor_paths_aligned] = AlignImages(nframes, path, PH_paths, fluor_paths, create_stack)
% ALIGNIMAGES: manually align the first 'nframes' images to frame nframes +
% 1
%
% Copyright (c) 2021 Silke Bergeler
%
% INPUT:
% - nframes: number of frames that are to be aligned (starting from frame 1), since 
% there is typically an initial shift of the cells
% - path: path to the directory with the data
% - PH_paths: paths to the PH images 
% - fluor_paths: paths to the fluorescence images
% - create_stack: Boolean, if true: an image stack with the
% aligned PH images is created to visually check the results

% OUTPUT:
% - align_shifts: array with the shifts of the fluorescence images
% - fluor_paths_aligned: paths to the aligned fluorescence images

% NOTE: the phase contrast images are used for alignment, because they are
% easier to compare to each other; we approximate the shift of the
% fluorescence images by the shift of the PH images

%% align the PH images

global IMAGEDIMX
global IMAGEDIMY
nfluoim = size(fluor_paths,1); % total number of fluorescence images

% if there is no file align_shifts.txt (i.e. images have been aligned
% already), align the images
if ~isfile([path, 'align_shifts.txt'])
    % array to save the shifts in x and y direction for the frames considered
    align_shifts = zeros(nframes, 2); 
    % the first 'nframes' frames are compared to frame 'nframes' + 1
    for frame = 1:nframes
        fixed = imread(PH_paths{nframes+1});
        moving = imread(PH_paths{frame});
        [xalign,yalign] = gui_align_images({fixed, moving, frame});
        align_shifts(frame,:) = [xalign,yalign];
    end

    % create an image stack of the aligned PH images to visually check the result
    if(create_stack)
        image_aligned = [];
        for frame = 1:nfluoim
            image = imread(PH_paths{frame});
            if frame <= nframes
                xalign = align_shifts(frame,1);
                yalign = align_shifts(frame,2);
                shifted_image = zeros(IMAGEDIMY, IMAGEDIMX,'uint16');
                x1 = max(1 + xalign,1);
                x2 = min(IMAGEDIMX, IMAGEDIMX + xalign);
                y1 = max(1 + yalign,1);
                y2 = min(IMAGEDIMY, IMAGEDIMY + yalign);
                xm1 = max(1 - xalign,1);
                xm2 = min(IMAGEDIMX, IMAGEDIMX - xalign);
                ym1 = max(1 - yalign,1);
                ym2 = min(IMAGEDIMY, IMAGEDIMY - yalign);
                shifted_image(y1:y2,x1:x2) = image(ym1:ym2,xm1:xm2);
                image_aligned = cat(3,image_aligned,shifted_image);
            else
                image_aligned = cat(3,image_aligned,image);
            end
        end
        saveimagestack(image_aligned, [path 'PH_images_aligned.tif'])
    end
    % write align_shifts to a file
    writetable(array2table(align_shifts,'VariableNames',{'x_shift' ,'y_shift'}),[path,'align_shifts.txt'])
else
    align_shifts = table2array(readtable([path, 'align_shifts.txt']));
end

if ~exist([path, 'Fluor_aligned'], 'dir')
    mkdir([path, 'Fluor_aligned'])
    % shift images and save them to folder Fluor_aligned
    for frame = 1:nfluoim
        image = imread(fluor_paths{frame});
        if frame <= nframes
            xalign = align_shifts(frame,1);
            yalign = align_shifts(frame,2);
            shifted_image = zeros(IMAGEDIMY, IMAGEDIMX,'uint16');
            x1 = max(1 + xalign,1);
            x2 = min(IMAGEDIMX, IMAGEDIMX + xalign);
            y1 = max(1 + yalign,1);
            y2 = min(IMAGEDIMY, IMAGEDIMY + yalign);
            xm1 = max(1 - xalign,1);
            xm2 = min(IMAGEDIMX, IMAGEDIMX - xalign);
            ym1 = max(1 - yalign,1);
            ym2 = min(IMAGEDIMY, IMAGEDIMY - yalign);
            shifted_image(y1:y2,x1:x2) = image(ym1:ym2,xm1:xm2);
        else
            shifted_image = image;
        end
        imwrite(shifted_image,[path,'Fluor_aligned/TXR_t',num2str(frame),'.tif']);
    end
end

% the paths to the aligned images (set the entries of fluor_paths_aligned to the paths of
% the aligned fluorescence images as defined before)
fluor_paths_aligned = arrayfun(@(x) [path,'Fluor_aligned/TXR_t',num2str(x),'.tif'], 1:nfluoim, 'UniformOutput', false)';

end