function [] = DensityPlotSpots(mat_paths, save_path, pixelsize)
% DENSITYPLOTSPOTS: plots the spot positions as X-Y-coordinates in an
% artificial cell as a heat map.
% We use the spots and cell meshes detected by Oufti and plot a heat map
% which shows the cluster positions (in % of cell length and width)
%
% Copyright (c) 2021 Silke Bergeler
% 
% Input:
% - mat_paths: paths to mat files
% - save_path: path to save figure
% - pixelsize: pixel size of images

lvec = []; % vector to save spot position along the centerline of the Oufti mesh
dvec = []; % vector to save spot position distance from centerline
celllength_vec = []; % vector to save the cell lengths of the cells with 1 spot
cellwidth_vec = []; % vector to save the cell widths of the cells with 1 spot
    
for m = 1:length(mat_paths)
    load(mat_paths{m},'cellList') % load the mat file with the mesh and spot information
    
    nframes = length(cellList.meshData); % number of frames
    
    for i = 1:nframes
        if(~isempty(cellList.meshData{i})) % only proceed if the frame contains at least one cell
            ncells = length(cellList.meshData{i}); % number of cells in frame i
            for k = 1:ncells
                disp(strcat('frame number=',num2str(i),', cell number=',num2str(k)))
                celldata = cellList.meshData{i}{k}; % the cell data set
                positions = extractfield(celldata.spots_matlab,'positions');
                if(~isempty(positions) && ~any(isnan(positions)) && size(celldata.mesh,2)==4) % only proceed if there is at least one spot and a proper cell mesh
                    x = mean([celldata.mesh(:,1), celldata.mesh(:,3)], 2);
                    y = mean([celldata.mesh(:,2), celldata.mesh(:,4)], 2);
                    steplength = (diff(x).^2 + diff(y).^2).^(1/2);
                    celllength = sum(steplength)*pixelsize; % length of the cell (along center line)
                    celllength_vec = [celllength_vec celllength]; % save the cell length to the vector
                    cellwidths = sqrt((celldata.mesh(:,1)-celldata.mesh(:,3)).^2+(celldata.mesh(:,2)-celldata.mesh(:,4)).^2)*pixelsize; % width of the cell along the center line
                    for j = 1:length(positions)
                        cellwidth = cellwidths(positions(j));
                        cellwidth_vec = [cellwidth_vec cellwidth]; % save the cell width to the vector
                        lvec = [lvec celldata.spots_matlab(j).l*pixelsize/celllength]; % add the spot position on center line (in % of cell length) to lvec
                        % there are lvec entries which are larger than 1, how can that be?
                        dvec = [dvec celldata.spots_matlab(j).d*pixelsize/cellwidth]; % add the spot position distance from centerline (in % of cell width) to dvec
                    end
                end
            end
        end
    end
end

% plot a density map of the positions
hist = hist3([lvec' dvec'],[50 10]);

figure;
colormap('hot')
imagesc(hist') % same as image, but data is scaled
colorbar
pbaspect([mean(celllength_vec) mean(cellwidth_vec) 1])
% reverse y direction
set(gca,'YDir','normal')

% save the figure as pdf
saveas(gcf,[save_path 'density_spots.pdf'])
close(gcf)

% figure;
% %surf(hist');
% surf(hist','EdgeColor', 'None', 'facecolor', 'interp');
% view(2);
% pbaspect([mean(celllength_vec) mean(cellwidth_vec) 1])
% colormap('hot')

end

