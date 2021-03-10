function [] = GetOrientationRelativeToCell(mat_path, xshift, yshift)
% GETORIENTATIONRELATIVETOCELL: calculate orientation of spot relative to cell orientation
%
% Copyright (c) 2021 Silke Bergeler
%
% Input: 
% - mat_path: path to mat file
% - xshift: shift of the mesh in x direction w.r.t the fluorescence image
% - yshift: shift of the mesh in y direction w.r.t the fluorescence image

masks = load(mat_path);
nframes = size(masks.cellList.meshData,2);

for frame = 1:nframes
    ncells = length(masks.cellList.meshData{frame});
    for cell = 1:ncells
        this_cell = masks.cellList.meshData{frame}{cell};
        if isfield(this_cell,'spots_matlab')
            this_cell_spots = this_cell.spots_matlab;
            mesh = this_cell.mesh;
            % change the meshes according to the shifts
            mesh(:,[1 3]) = mesh(:,[1 3]) + xshift;
            mesh(:,[2 4]) = mesh(:,[2 4]) + yshift;
            centerline = [mean([mesh(:,1) mesh(:,3)],2) mean([mesh(:,2) mesh(:,4)],2)];
            rel_orient = []; % orientation of spot relative to cell orientation
            if isempty(this_cell_spots(1).l)
                masks.cellList.meshData{frame}{cell}.spots_matlab.rel_orient = [];
                continue
            end
            for n = 1:length(this_cell_spots) % for all spots
                % use the centerline of the cell to get orientation of the cell
                % close to the cluster
                spot_position = this_cell_spots(n).positions;
                centerline_selected = centerline(max((spot_position-2),1):min((spot_position+2),length(centerline)),:);
                
                % get line with shortest orthogonal distance to the points
                % using PCA
                [coeff,~,~] = pca(centerline_selected);
                dirVect = coeff(:,1);
                
                % calculate orientation of fitted line to x axis
                cell_orient = -rad2deg(atan(dirVect(2)/dirVect(1)));
                
                % calculate orientation of spot relative to cell
                % orientation
                rel_orient_val = TransformAngle(this_cell_spots(n).spot_orient - cell_orient);
                rel_orient = [rel_orient rel_orient_val];
            end
            if isnan(rel_orient)
                disp('rel_orient is nan!')
            end
            rel_orient = num2cell(rel_orient);
            [masks.cellList.meshData{frame}{cell}.spots_matlab.rel_orient] = rel_orient{:};
        end
    end
end

save(mat_path,'-struct','masks');

end