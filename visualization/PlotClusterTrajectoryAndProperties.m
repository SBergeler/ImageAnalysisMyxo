function [ ] = PlotClusterTrajectoryAndProperties(mat_paths, fluor_paths, save_path, pixelsize, xshift, yshift)
% PLOTCLUSTERTRAJECTORYANDPROPERTIES: Plot the cluster position over time
%
% Copyright (c) 2021 Silke Bergeler
%
% Input:
% - mat_paths: paths to mat files
% - fluor_paths: paths to fluorescence images
% - save_path: path to save figures
% - pixelsize: pixel size of images
% - xshift: how many pixels the meshes have to be shifted in x direction to
% match the fluorescence signals
% - yshift: how many pixels the meshes have to be shifted in y direction to 
% match the fluorescence signals

% if folder to save the images already exists, remove it and create new one
newdir = strcat([save_path 'single_cell_stacks/']);
if exist(newdir, 'dir')
    rmdir(newdir,'s')
end
mkdir(newdir);

for j = 1:length(mat_paths)
    masks = load(mat_paths{j}); % load the mat file with the mesh and cluster position information
    
    cellids = masks.cellList.cellId; % the cell ids for each frame
    cellids = cellids(~cellfun('isempty',cellids)); % remove empty cells (if a frame does not contain any cell, which is e.g. the case when cell detection was not run over all frames)
    cellids_mat = cell2mat(cellids); % transform cell array to vector
    cellids_unique = unique(cellids_mat); % the unique cell ids
    cellids_with_counts = table(cellids_unique',histcounts(cellids_mat, [unique(cellids_mat) max(unique(cellids_mat))+1])','VariableNames',{'cell_id','no_frames'}); % matrix which stores the cell ids and the number of frames the cell is analyzed
    
    nframes = length(cellids); % number of frames (which contain cells)
    if nframes == 1
        disp('Only one frame --> no movie is generated')
        return
    end
    
    figure;
    histogram(cellids_with_counts.no_frames,[(min(cellids_with_counts.no_frames):max(cellids_with_counts.no_frames))-0.5 max(cellids_with_counts.no_frames)+0.5]) % histogram which shows the number of cells that exists over a certain number of frames
    xlabel('Number of frames the cell is detected')
    ylabel('Count')
    title('Distribution of time traces')
    
    % save the figure as pdf
    saveas(gcf,[save_path 'Oufti_cell_detection_distribution_time_traces_dataset_' num2str(j) '.png'])
    close(gcf)
    
    % only consider cells that are identified in all frames
    cellids_unique_allframes = cellids_with_counts.cell_id(cellids_with_counts.no_frames == nframes);
    
    % create matrix with the unique ids of the cells that occur in all frames and the corresponding positions they
    % refer to in the frames
    cellids_numbers_mat = zeros(length(cellids_unique_allframes),nframes);
    for i = 1:length(cellids_unique_allframes) % for all unique cell ids
        for k = 1:nframes % search in all frames at which position the data for this cell id can be found
            index = find(cellids{k} == cellids_unique_allframes(i));
            cellids_numbers_mat(i,k) = index;
            if(isempty(index))
                disp('TrackClusterPosition: no index found, should not be possible!');
                break;
            end
        end
    end
    
    cluster_rel_mat = NaN(length(cellids_unique_allframes),nframes,2); % matrix which saves the cluster position (relative to the cell boundary) over time
    cluster_size_mat = NaN(length(cellids_unique_allframes),nframes); % matrix which saves the cluster size over time
    cluster_int_mat = NaN(length(cellids_unique_allframes),nframes); % matrix which saves the relative intensity in the cluster over time
    box_dims = zeros(length(cellids_unique_allframes),nframes,4); % save the positions of a box where the cell is found (from Oufti)
    
    % get the cluster positions and sizes
    celllength_vec = zeros(1,length(cellids_unique_allframes)*nframes); % save all cell lengths to a vector (used for plotting)
    cellwidth_vec = zeros(1,length(cellids_unique_allframes)*nframes); % save all cell widths to a vector (used for plotting)
    for i = 1:length(cellids_unique_allframes)
        for k = 1:nframes
            cellid = cellids_numbers_mat(i,k); % get the id of the cell in frame k
            if(cellid == 0) % if the id is zero, the cell does not exist in frame k (should not be the case)
                disp('TrackClusterPosition: cellid = 0, this should not be the case!');
                break;
            end
            this_cell = masks.cellList.meshData{k}{cellid}; % the cell data set
            this_cell_spots = this_cell.spots_matlab;
            box_dims(i,k,:) = this_cell.box;
            celllength = this_cell.length*pixelsize; % length of the cell (in um)
            celllength_vec((i-1)*nframes + k) = celllength;
            cellwidth = this_cell.width*pixelsize; % width of the cell (in um)
            cellwidth_vec((i-1)*nframes + k) = cellwidth;
            if(length(this_cell_spots)==1 && ~isempty(this_cell_spots(1).l)) % only proceed if there is exactly one spot in the cell
                cluster_rel_mat(i,k,:) = [this_cell_spots.l*pixelsize/celllength this_cell_spots.d*pixelsize/cellwidth];
                cluster_size_mat(i,k,:) = this_cell_spots.spot_area*pixelsize^2;
                cluster_int_mat(i,k,:) = this_cell_spots.rel_signal_int;
            end
        end
    end
    % make a stack for each cell with the detected cluster position
    % chose the box around the cell such that the size is big enough to fit the
    % cell in all time frames
    box_per_cell = [min(box_dims(:,:,1),[],2) min(box_dims(:,:,2),[],2) ...
        max(box_dims(:,:,1)+box_dims(:,:,3),[],2) - min(box_dims(:,:,1),[],2) ...
        max(box_dims(:,:,2)+box_dims(:,:,4),[],2) - min(box_dims(:,:,2),[],2)];
    
    t = linspace(0,2*pi,50);
    
    for i = 1:length(cellids_unique_allframes)
        mkdir([newdir num2str(i) '/'])
        disp(['Creating images of cell ',num2str(i), ' out of ', num2str(length(cellids_unique_allframes)), ' ...']);
        % prepare the new video
        % uncomment the next line to automatically close the previous figures
        % close all; % close all images
        % vidObj = VideoWriter([newdir 'cell_' num2str(i) '.avi']);
        % open(vidObj);
        for k = 1:nframes
            figure('visible','off');
            set(gcf,'Unit','Inches','position',[0,0,4,4],'PaperUnits', 'Inches', 'PaperSize', [4, 4])
            
            subplot(2,3,[1 2])
            fluorescence_image = imread(fluor_paths{k});
            cellid = cellids_numbers_mat(i,k); % get the id of the cell in frame k
            this_cell = masks.cellList.meshData{k}{cellid}; % the cell data set
            this_cell_spots = this_cell.spots_matlab;
            box = box_per_cell(i,:); % xmin, ymin, width, height
            %imwrite(fluor_stack{k},strcat(string{1},'_single_cell_stacks/cell_',num2str(i),'.tiff'),'WriteMode','append')
            imshow(imadjust(fluorescence_image(box(2):(box(2)+box(4)), box(1):(box(1)+box(3)))));
            title(['Cell id = ', num2str(cellids_unique_allframes(i)), ', frame no = ', num2str(k)])
            hold on
            % plot the mesh of the cell
            mesh = this_cell.mesh;
            % change the meshes according to the shifts and adjust the coordinates to the zoomed image
            mesh(:,[1 3]) = mesh(:,[1 3]) - box(1) + xshift;
            mesh(:,[2 4]) = mesh(:,[2 4]) - box(2) + yshift;
            plot(mesh(:,1),mesh(:,2),'Color','cyan','LineWidth',1.2);
            plot(mesh(:,3),mesh(:,4),'Color','cyan','LineWidth',1.2);
            hold on
            % plot center line
            centerx = mean([mesh(:,1) mesh(:,3)],2);
            centery = mean([mesh(:,2) mesh(:,4)],2);
            plot(centerx, centery,'Color','cyan','LineWidth',1.2);
            hold on
            % plot the cluster positions (all identified clusters are shown)
            if ~isempty(this_cell_spots(1).l)
                for n = 1:length(this_cell_spots)
                    % plot ellipses fitted to spots
                    a = this_cell_spots(n).spot_major_axis/2;
                    b = this_cell_spots(n).spot_minor_axis/2;
                    Xc = this_cell_spots(n).x-box(1);
                    Yc = this_cell_spots(n).y-box(2);
                    phi = deg2rad(-this_cell_spots(n).spot_orient);
                    x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
                    y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
                    % plot spot position
                    plot(Xc, Yc, 'xk', 'MarkerSize',10, 'LineWidth', 1.2)
                    hold on
                    % plot ellipse
                    plot(x,y,'magenta','Linewidth',1.2)
                    hold on
                end
            end
            
            subplot(2,3,[4 5])
            plot(cluster_rel_mat(i,1:k,1),cluster_rel_mat(i,1:k,2),'-','color',[0 0 0 0.5],'MarkerSize',3);
            xlim([0,1]);
            ylim([-0.5,0.5]); % d can have negative and positive values
            hold on
            line([0.5 0.5],[-0.5 0.5],'Color','black','LineStyle','--');
            hold on
            line([0 1], [0 0],'Color','black','LineStyle','--');
            pbaspect([mean(celllength_vec) mean(cellwidth_vec) 1]);
            set(gca,'FontSize',6)
            
            subplot(2,3,3)
            plot(1:k,cluster_size_mat(i,1:k),'-','color',[0 0 0 0.5],'MarkerSize',3);
            xlim([1,nframes]);
            ylim([0,max(max(cluster_size_mat(i,:)),0.5)]);
            xlabel('frame no.');
            ylabel('cluster area [um^2]');
            pbaspect([1 1 1])
            xticks((1:5)*nframes/5)
            set(gca,'FontSize',6)
            
            subplot(2,3,6)
            plot(1:k,cluster_int_mat(i,1:k),'-','color',[0 0 0 0.5],'MarkerSize',3);
            xlim([1,nframes]);
            ylim([0,1]);
            xlabel('frame no.');
            ylabel('rel. cluster int.');
            pbaspect([1 1 1])
            xticks((1:5)*nframes/5)
            set(gca,'FontSize',6)
            
            set(findall(gcf,'type','text'),'FontSize',6)
            saveas(gcf, [newdir 'dataset_' num2str(j) '_cellid_' num2str(cellids_unique_allframes(i)) '/frame_' num2str(k) '.tif']);
            close(gcf)
            
            % currFrame = getframe(gcf);
            % writeVideo(vidObj,currFrame);
        end
        % close(vidObj);
    end
end
disp('done');
end
