function [] = PlotCellAndSpots(this_cell,image_fluor,PixelIdxList,nframe,cell_id,spot_positions,spot_area,signal_int)
% PLOTCELLANDSPOTS: plot a cell together with its spots
%
% Copyright (c) 2021 Silke Bergeler
%
% Input:
% - this_cell: data for cell of interest
% - image_fluor: fluoresence image (frame 'nframe')
% - PixelIdxList: locations of pixels in the spot regions
% - nframe: id of the frame plotted
% - cell_id: id of the cell plotted
% - spot_positions: positions of the spots detected
% - spot_area: number of pixels in the spot region
% - signal_int: summed fluoresence signal intensity of the spots detected

box = this_cell.box;
topx = box(1);
topy = box(2);
width_box = box(3);
height_box = box(4);

box_xmin = max(1,topx);
box_xmax = min(size(image_fluor,2),topx+width_box);
box_ymin = max(1,topy);
box_ymax = min(size(image_fluor,1),topy+height_box);

colors = [228,26,28;55,126,184;77,175,74;152,78,163;255,127,0;255,255,51;166,86,40;247,129,191;153,153,153]/255;
if size(spot_positions,1) > 9
    disp('PlotCellAndSpots:not enough colors, more than 9 spots detected!')
end

figure('Position',[500,500,600,600])
image_fluor_sel = image_fluor(box_ymin:box_ymax,box_xmin:box_xmax);
imshow(imadjust(image_fluor_sel),'InitialMagnification','fit')
hold on
% plot cell surrounding
plot(this_cell.mesh(:,1)-topx+1,this_cell.mesh(:,2)-topy+1,'Color','cyan','MarkerSize',10);
plot(this_cell.mesh(:,3)-topx+1,this_cell.mesh(:,4)-topy+1,'Color','cyan','MarkerSize',10);
hold on
% plot spots
str = cell(1,size(spot_positions,1));
p = zeros(1,size(spot_positions,1));
for i = 1:size(spot_positions,1)
    p(i) = plot(spot_positions(i,1)-topx+1,spot_positions(i,2)-topy+1,'x','Color',colors(i,:),'MarkerSize',10,'LineWidth',1.5);
    % text(spot_positions(i,1)-topx+5,spot_positions(i,2)-topy,['Spot ID: ' num2str(i)],'Color','white');
    str{i} = ['Spot ID: ' num2str(i) newline 'Spot area: ' num2str(spot_area(i)) newline 'Signal int: ' num2str(signal_int(i)) newline '-----'];
    % plot spot boundaries
    [row,col] = ind2sub(size(image_fluor),PixelIdxList{i});
    bin_tmp = logical(false(size(image_fluor_sel)));
    for j = 1:length(row)
        bin_tmp(row(j)-topy+1,col(j)-topx+1) = true;
    end
    B = bwboundaries(bin_tmp,'noholes');
    ps = plot(B{1}(:,2), B{1}(:,1),'Color',colors(i,:),'LineWidth', 1.5);
    ps.Color(4) = 0.5;
end
title(['Frame: ' num2str(nframe), ', cell id: ' num2str(cell_id)],'FontSize',14)
legend(p,str,'Location','northeastoutside','FontSize',14)
end