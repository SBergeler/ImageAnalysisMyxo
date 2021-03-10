function [] = PlotCellAndSpotsTL(this_cell,previous_cell,image_fluor,PixelIdxList,previous_image_fluor,previous_PixelIdxList,nframe,cell_id,...
    spot_positions,spot_area,signal_int,previous_spot_positions,previous_spot_area,previous_signal_int)
% PLOTCELLANDSPOTSTL: plot a cell in the frame 'nframe' and the same cell
% in the previous frame together with the spots detected within the cell
%
% Copyright (c) 2021 Silke Bergeler
%
% Input:
% - this_cell: data for cell of interest in frame 'nframe'
% - previous_cell: data for cell of interest in frame 'nframe'-1
% - image_fluor: fluoresence image (frame 'nframe')
% - PixelIdxList: locations of pixels in the spot regions
% - previous_image_fluor: fluoresence image (frame 'nframe'-1)
% - previous_PixelIdxList: locations of pixels in the spot regions
% (previous frame)
% - nframe: id of the frame plotted
% - cell_id: id of the cell plotted
% - spot_positions: positions of the spots detected
% - spot_area: number of pixels in the spot region
% - signal_int: summed fluoresence signal intensity of the spots detected
% - previous_spot_positions: positions of the spots detected (previous
% frame)
% - previous_spot_area: number of pixels in the spot region (previous
% frame)
% - previous_signal_int: summed fluoresence signal intensity of the spots
% detected (previous frame)


box = this_cell.box;
topx = box(1);
topy = box(2);
width_box = box(3);
height_box = box(4);

box_xmin = max(1,topx);
box_xmax = min(size(image_fluor,2),topx+width_box);
box_ymin = max(1,topy);
box_ymax = min(size(image_fluor,1),topy+height_box);

box_previous = previous_cell.box;
topx_previous = box_previous(1);
topy_previous = box_previous(2);
width_box_previous = box_previous(3);
height_box_previous = box_previous(4);

box_xmin_previous = max(1,topx_previous);
box_xmax_previous = min(size(image_fluor,2),topx_previous+width_box_previous);
box_ymin_previous = max(1,topy_previous);
box_ymax_previous = min(size(image_fluor,1),topy_previous+height_box_previous);

colors = [228,26,28;55,126,184;77,175,74;152,78,163;255,127,0;255,255,51;166,86,40;247,129,191;153,153,153]/255;
if size(previous_spot_positions,1) > 9
    disp('PlotCellAndSpotsTL:not enough colors, more than 9 spots detected!')
end

figure('Position',[500,500,1000,600])
image_fluor_sel = image_fluor(box_ymin:box_ymax,box_xmin:box_xmax);
previous_image_fluor_sel = previous_image_fluor(box_ymin_previous:box_ymax_previous,box_xmin_previous:box_xmax_previous);

subplot(1,2,1)
imshow(imadjust(previous_image_fluor_sel),'InitialMagnification','fit')
hold on
% plot cell surrounding
plot(previous_cell.mesh(:,1)-topx_previous,previous_cell.mesh(:,2)-topy_previous,'Color','cyan','MarkerSize',10);
plot(previous_cell.mesh(:,3)-topx_previous,previous_cell.mesh(:,4)-topy_previous,'Color','cyan','MarkerSize',10);
hold on
% plot spots
str = cell(1,size(previous_spot_positions,1));

p = zeros(1,size(previous_spot_positions,1));
for i = 1:size(previous_spot_positions,1)
    p(i) = plot(previous_spot_positions(i,1)-topx_previous,previous_spot_positions(i,2)-topy_previous,'x','Color',colors(i,:),'MarkerSize',10,'LineWidth',1.5);
    %text(previous_spot_positions(i,1)-topx_previous+5,previous_spot_positions(i,2)-topy_previous,['Spot ID: ' num2str(i)],'Color','r','FontWeight','Bold','FontSize',14);
    %x = [30 previous_spot_positions(i,1)-topx_previous];
    %y = [40 previous_spot_positions(i,2)-topy_previous];
    %[figx, figy] = axxy2figxy(gca, x, y);
    %annotation('textarrow',figx,figy,'String',['Spot ID: ' num2str(i)])
    str{i} = ['Spot ID: ' num2str(i) newline 'Spot area: ' num2str(previous_spot_area(i)) newline 'Signal int: ' num2str(previous_signal_int(i)) newline '-----'];
    % plot spot boundaries
    [row,col] = ind2sub(size(previous_image_fluor),previous_PixelIdxList{i});
    bin_tmp = logical(false(size(previous_image_fluor_sel)));
    for j = 1:length(row)
        bin_tmp(row(j)-topy_previous,col(j)-topx_previous) = true;
    end
    B = bwboundaries(bin_tmp,'noholes');
    ps = plot(B{1}(:,2), B{1}(:,1),'Color',colors(i,:),'LineWidth', 1.5);
    ps.Color(4) = 0.5;
end
title('Previous frame','FontSize',14)
legend(p,str,'Location','northeastoutside','FontSize',14)

subplot(1,2,2)
imshow(imadjust(image_fluor_sel),'InitialMagnification','fit')
hold on
% plot cell surrounding
plot(this_cell.mesh(:,1)-topx,this_cell.mesh(:,2)-topy,'Color','cyan','MarkerSize',10);
plot(this_cell.mesh(:,3)-topx,this_cell.mesh(:,4)-topy,'Color','cyan','MarkerSize',10);
hold on
% plot spots
str = cell(1,size(spot_positions,1));
p = zeros(1,size(spot_positions,1));
for i = 1:size(spot_positions,1)
    p(i) = plot(spot_positions(i,1)-topx,spot_positions(i,2)-topy,'x','Color',colors(i,:),'MarkerSize',10,'LineWidth',1.5);
    %text(spot_positions(i,1)-topx+5,spot_positions(i,2)-topy,['Spot ID: ' num2str(i)],'Color','r','FontWeight','Bold','FontSize',14);
    str{i} = ['Spot ID: ' num2str(i) newline 'Spot area: ' num2str(spot_area(i)) newline 'Signal int: ' num2str(signal_int(i)) newline '-----'];
    % plot spot boundaries
    [row,col] = ind2sub(size(image_fluor),PixelIdxList{i});
    bin_tmp = logical(false(size(image_fluor_sel)));
    for j = 1:length(row)
        bin_tmp(row(j)-topy,col(j)-topx) = true;
    end
    B = bwboundaries(bin_tmp,'noholes');
    ps = plot(B{1}(:,2), B{1}(:,1),'Color',colors(i,:),'LineWidth', 1.5);
    ps.Color(4) = 0.5;
end
title(['Frame: ' num2str(nframe), ', cell id: ' num2str(cell_id)],'FontSize',14)
legend(p,str,'Location','northeastoutside','FontSize',14)

end