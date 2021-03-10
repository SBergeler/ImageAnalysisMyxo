function [] = AddLengthArea(mat_path)
% ADDLENGTHAREA: add cell length, width and area to an existing mat file
% from Oufti
% Note: Oufti calculates these quantities, but does not seem to save them in the
% mat file. This function calculates the quantities based on the mesh and 
% adds them to the mat file. 
%
% Copyright (c) 2021 Silke Bergeler
%
% INPUT: 
% - mat_path: path to mat file

masks = load(mat_path);
nframes = length(masks.cellList.meshData);
for frame = 1:nframes
   ncells = length(masks.cellList.meshData{frame});
   for cell = 1:ncells
       mesh = masks.cellList.meshData{frame}{cell}.mesh;
       if length(mesh) >= 4
           centerline = [(mesh(:,1)+mesh(:,3))/2. (mesh(:,2)+mesh(:,4))/2.]; % x and y coordinates of centerline points
           steps_centerline = diff(centerline);
           length_steps_centerline = sqrt(sum(steps_centerline.^2,2));
           cell_length = sum(length_steps_centerline);
           widths = sqrt((mesh(:,1)-mesh(:,3)).^2+(mesh(:,2)-mesh(:,4)).^2);
           avg_widths = mean([widths(1:(length(widths)-1)) widths(2:length(widths))],2);
           area = sum(avg_widths.*length_steps_centerline);
           masks.cellList.meshData{frame}{cell}.length = cell_length;
           masks.cellList.meshData{frame}{cell}.area = area;
           masks.cellList.meshData{frame}{cell}.width = mean(widths(round(0.3*length(widths)):round(0.7*length(widths))));
       end
   end
end
save(mat_path,'-struct','masks') % overwrite the mat file
end