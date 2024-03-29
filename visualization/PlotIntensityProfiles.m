function [] = PlotIntensityProfiles(mat_paths, fluor_paths, save_path, xshift, yshift, pixelsize)
% PLOTINTENSITYPROFILES: Plot normalized intensity profiles along cells in snapshot images
%
% Copyright (c) 2021 Silke Bergeler
%
% Input:
% - mat_paths: paths to mat files
% - fluor_paths: cell array with paths to fluorescence images
% - save_path: path to save figures
% - xshift: how many pixels the meshes have to be shifted in x direction to
% match the fluorescence signals
% - yshift: how many pixels the meshes have to be shifted in y direction to 
% match the fluorescence signals
% - pixelsize: pixel size of images

disp('Plotting normalized intensity profiles ... ');

% save the length of the cell and the integrated signal intensities in the cells (integrated along the short cell axis)
cell_struct = struct('dataset',[],'frame',[],'cell',[],'length',[],'no_segments',[],'summed_signal',[],'conc',[],'diffs',[]); 

i = 1;

for m = 1:length(mat_paths)
    load(mat_paths{m},'cellList');
    nframes = size(cellList.meshData,2);
    
    for frame = 1:nframes
        ncells = size(cellList.meshData{frame},2);
        image_fluor = imread(fluor_paths{m}{frame});
        background_fluor = getBackgroundImage(image_fluor,cellList.meshData{frame},xshift{m},yshift{m});
        for cell = 1:ncells
            this_cell = cellList.meshData{frame}{cell};
            mesh = this_cell.mesh;
            box = this_cell.box;
            % change the box according to the shifts
            box(1) = box(1) + xshift{m};
            box(2) = box(2) + yshift{m};
            % consider only a smaller image with the cell to reduce
            % computation time
            image_dim_x = size(image_fluor,2);
            image_dim_y = size(image_fluor,1);
            topx = box(1);
            topy = box(2);
            width_box = box(3);
            height_box = box(4);
            box_xmin = max(1,topx);
            box_xmax = min(image_dim_x,topx + width_box);
            box_ymin = max(1,topy);
            box_ymax = min(image_dim_y,topy + height_box);
            image_fluor_selected = image_fluor(box_ymin:box_ymax,box_xmin:box_xmax);
            background_fluor_selected = background_fluor(box_ymin:box_ymax,box_xmin:box_xmax);
            % change the meshes according to the shifts
            mesh(:,[1 3]) = mesh(:,[1 3]) - box(1) + xshift{m};
            mesh(:,[2 4]) = mesh(:,[2 4]) - box(2) + yshift{m};
            cell_length = this_cell.length*pixelsize;
            % get lengths of mesh segments (used for plotting later)
            centerline_x = mean([mesh(:,1), mesh(:,3)],2);
            centerline_y = mean([mesh(:,2), mesh(:,4)],2);
            diff_x = diff(centerline_x);
            diff_y = diff(centerline_y);
            diffs = sqrt(diff_x.^2 + diff_y.^2);
            
            % Create binary mask
            mesh_min_x = min(min(mesh(:,1),min(mesh(:,3))));
            mesh_max_x = max(max(mesh(:,1),max(mesh(:,3))));
            mesh_min_y = min(min(mesh(:,2),min(mesh(:,4))));
            mesh_max_y = max(max(mesh(:,2),max(mesh(:,4))));
            if mesh_min_x < 1
                mesh_min_x = 1;
            end
            if mesh_min_y < 1
                mesh_min_y = 1;
            end
            if mesh_max_x > size(image_fluor_selected,2)
                mesh_max_x = size(image_fluor_selected,2);
            end
            if mesh_max_y > size(image_fluor_selected,1)
                mesh_max_y = size(image_fluor_selected,1);
            end
            % x = double([mesh(:,1);mesh(:,3)])';
            % y = double([mesh(:,2);mesh(:,4)])';
            % Convert region-of-interest polygon to mask (pixels inside polygon are 1)
            % binary_mask = poly2mask(x,y,size(image_fluor_selected,1),size(image_fluor_selected,2)); 
            
            % Substract background from cell signal
            % minimal box that includes cell:
            box_xmin = uint16(mesh_min_x);
            box_xmax = uint16(mesh_max_x);
            box_ymin = uint16(mesh_min_y);
            box_ymax = uint16(mesh_max_y);
            
%             % smaller image containing only the cell of interest
%             zoomed_fluor = image_fluor_selected(uint16(mesh_min_y):uint16(mesh_max_y),uint16(mesh_min_x):uint16(mesh_max_x)); 
%             % cut out the cell to analyze the background
%             binary_mask_small = binary_mask(uint16(mesh_min_y):uint16(mesh_max_y),uint16(mesh_min_x):uint16(mesh_max_x));
%             background = median(zoomed_fluor(~binary_mask_small));
%             % subtract background from fluor image (values are rounded to next int and negative values are set to 0)
%             image_fluor_selected = image_fluor_selected - background; 
            
            % get fluorescence image with 0's inside the cells that are detected with
            % Oufti
            background_fluor_small = background_fluor_selected(box_ymin:box_ymax,box_xmin:box_xmax);
            
            % if number of pixels in the background is smaller than 1000, increase box
            % size by a factor of 2 (up to size of box, so might be smaller than 2x)
            if(sum(sum(background_fluor_small>0)) < 1000)
                box_in = zeros(1,4);
                box_in(1) = box_xmin;
                box_in(2) = box_ymin;
                box_in(3) = box_xmax - box_xmin;
                box_in(4) = box_ymax - box_ymin;
                box_out = ChangeBoxSize(box_in, 2);
                topx = box_out(1);
                topy = box_out(2);
                width_box = box_out(3);
                height_box = box_out(4);
                box_xmin = max(1,topx);
                box_xmax = min(size(background_fluor_selected,2),topx+width_box);
                box_ymin = max(1,topy);
                box_ymax = min(size(background_fluor_selected,1),topy+height_box);
                background_fluor_small = background_fluor_selected(box_ymin:box_ymax,box_xmin:box_xmax);
            end
            
            % define background as the median non-zero signal outside the cell in the box
            background = median(background_fluor_small(background_fluor_small>0));
            
            % subtract background from fluor image (values are rounded to next int and negative values are set to 0)
            image_fluor_selected = image_fluor_selected - background;
            
%             % for debugging only:
%             figure()
%             subplot(1,2,1)
%             imshow(imadjust(image_fluor_selected));
%             hold on
%             plot(mesh(:,1),mesh(:,2),'Color','cyan','LineWidth',1.2);
%             plot(mesh(:,3),mesh(:,4),'Color','cyan','LineWidth',1.2);
%             title('BG corrected')
%             subplot(1,2,2)
%             imshow(imadjust(background_fluor_small))
%             title('BG')
            
            % Convert region-of-interest polygon to mask (pixels inside polygon
            % are 1) and sum over signal of these pixels
            signal = arrayfun(@(x) ...
                (sum(sum(image_fluor_selected.*uint16(poly2mask(double([mesh(x,1) mesh(x+1,1) mesh(x+1,3) mesh(x,3)]),...
                double([mesh(x,2) mesh(x+1,2) mesh(x+1,4) mesh(x,4)]),size(image_fluor_selected,1),size(image_fluor_selected,2)))))), ...
                1:(size(mesh,1)-1));
            n_pixel = arrayfun(@(x) sum(sum(uint16(poly2mask(double([mesh(x,1) mesh(x+1,1) mesh(x+1,3) mesh(x,3)]),...
                double([mesh(x,2) mesh(x+1,2) mesh(x+1,4) mesh(x,4)]),size(image_fluor_selected,1),size(image_fluor_selected,2))))), ...
                1:(size(mesh,1)-1));
            % set segments with no pixels to 1 (since the signal is zero
            % there, the concentration becomes zero)
            n_pixel(n_pixel == 0) = 1; 
            conc = signal./n_pixel;
            no_segments = length(signal);
            % get number of pixels per segment to find empty ones
            % no_pixels = arrayfun(@(x) ...
            %    (sum(sum(uint16(poly2mask(double([mesh(x,1) mesh(x+1,1) mesh(x+1,3) mesh(x,3)]),...
            %    double([mesh(x,2) mesh(x+1,2) mesh(x+1,4) mesh(x,4)]),size(image_fluor_selected,1),size(image_fluor_selected,2)))))), ...
            %    1:(size(mesh,1)-1));
            % take only the signal values which were obtained from summing at
            % least one pixel
            % signal = signal(no_pixels > 0);
            % diffs = diffs(no_pixels > 0);
            
            % orient the cells such that the highest intensity value is always on
            % the same side
            [~, ind] = max(signal);
            if(ind > round(no_segments/2))
                signal = flip(signal);
                conc = flip(conc);
            end
            cell_struct(i).dataset = m;
            cell_struct(i).frame = single(frame);
            cell_struct(i).cell = cellList.cellId{frame}(cell);
            cell_struct(i).length = double(cell_length);
            cell_struct(i).no_segments = no_segments;
            cell_struct(i).summed_signal = transpose(signal);
            cell_struct(i).conc = transpose(conc);
            cell_struct(i).diffs = double(diffs);
            i = i + 1;
        end
    end
end

%%
% reshape the data into a table
cell_array = struct2cell(cell_struct);
sz = size(cell_array);
cell_array = reshape(cell_array, sz(1), []);
cell_array = cell_array'; 

cell_table = cell2table(cell_array);
cell_table.Properties.VariableNames = {'dataset' 'frame' 'cell' 'length' 'no_segments' 'summed_signal' 'conc' 'diffs'};

% some cells have a 'weird' mesh in the sense that some segments are a lot
% bigger than usual, this could be due to atypical meshes from Oufti; we
% display cells with atypically large segments
weird_cells = cellfun(@(x) any(x >= 1.5), cell_table.diffs);
frame_weird_cells = cell_table.frame(weird_cells);
cellid_weird_cells = cell_table.cell(weird_cells);
diffs_weird_cells = cellfun(@max, cell_table.diffs(weird_cells));
disp('weird cells')
for i=1:length(frame_weird_cells)
    disp(['frame ' num2str(frame_weird_cells(i)) ' cellid ' num2str(cellid_weird_cells(i)) ' dx ' num2str(diffs_weird_cells(i))])
end

% get maximal and mean difference between two segments
diffs_mat = cell2mat(cell_table.diffs);
max_dx = max(diffs_mat)*pixelsize;
% mean_dx = mean(diffs_mat)*pixelsize;

% add column with number of mesh segments (with at least one pixel)
% cell_table.no_mesh_segs = cellfun(@(x) size(x,2), cell_array(:,4));

%% diffs -> position
for i = 1:height(cell_table)
    pos = cell_table.diffs{i};
    int = cell_table.summed_signal{i}; 
    conc = cell_table.conc{i}; 
    pos = cumsum(pos);
    pos = pos - cell_table.diffs{i}/2;
    pos = pos*pixelsize;
    pos = pos./pos(end); % normalize
    cell_table.pos_x{i} = pos;
    % interpolate profiles such that each profile is evaluated at 100 points
    cell_table.summed_signal_interpolated{i} = interp1(pos,int,linspace(0,1,100),'spline');
    cell_table.conc_interpolated{i} = interp1(pos,conc,linspace(0,1,100),'spline');
end

% calculate average profile
avg_profile = vertcat(cell_table.summed_signal_interpolated{:});
avg_profile = mean(avg_profile,1);

avg_conc_profile = vertcat(cell_table.conc_interpolated{:});
avg_conc_profile = mean(avg_conc_profile,1);

%% Plot normalized intensity plots (summed intensities)

figure;

for i = 1:height(cell_table)
    p = plot(cell_table.pos_x{i},cell_table.summed_signal{i},'k');
    p.Color(4) = 0.2;
    hold on
end
% add average profile
plot(linspace(0,1,100),avg_profile,'red')

ylabel('summed signal intensity')
xlabel('x-position (normalized)')

% save the figure as pdf
set(gcf,'Unit','Inches','position',[0,0,5,3],'PaperUnits', 'Inches', 'PaperSize', [5, 3])
saveas(gcf,[save_path 'summed_signal_intensity_profile.pdf'])

% save as matlab figure
savefig([save_path 'summed_signal_intensity_profile.fig'])
close(gcf)

disp('done')

%% Plot normalized intensity plots (concentrations)

figure;

for i = 1:height(cell_table)
    p = plot(cell_table.pos_x{i},cell_table.conc{i},'k');
    p.Color(4) = 0.2;
    hold on
end
% add average profile
plot(linspace(0,1,100),avg_conc_profile,'red')

ylabel('intensity per pixel')
xlabel('x-position (normalized)')

% save the figure as pdf
set(gcf,'Unit','Inches','position',[0,0,5,3],'PaperUnits', 'Inches', 'PaperSize', [5, 3])
saveas(gcf,[save_path 'signal_per_pixel_intensity_profile.pdf'])

% save as matlab figure
savefig([save_path 'signal_per_pixel_intensity_profile.fig'])
close(gcf)

disp('done')

end