function [] = PlotSpotStatistics(mat_paths, save_path, pixelsize)
% PLOTSPOTSTATISTICS: this function plots statistics about the spots 
% (largest and smallest extension, eccentricity, ...)
% made for stacks of independent images
% only cells with one spot are considered
%
% Copyright (c) 2021 Silke Bergeler
%
% Input:
% - mat_paths: paths to mat files
% - save_path: path to save figures
% - pixelsize: pixel size of images

spot_major_axis = [];
spot_minor_axis = [];
spot_ecc = [];
spot_area = [];
rel_orient = [];
rel_signal_int = [];
signal_int = [];
cell_length = [];
cluster_position = []; % cluster position relative to cell length
    
for j = 1:length(mat_paths)
    load(mat_paths{j},'cellList');
    nframes = length(cellList.meshData);
    
    for frame = 1:nframes
        ncells = length(cellList.meshData{frame});
        for cell = 1:ncells
            this_cell = cellList.meshData{frame}{cell};
            if isfield(this_cell,'spots_matlab')
                if length(this_cell.spots_matlab) == 1 && ~isempty(this_cell.spots_matlab.l) % only consider cells with one spot
                    this_cell_spots = this_cell.spots_matlab;
                    spot_major_axis = [spot_major_axis this_cell_spots.spot_major_axis];
                    spot_minor_axis = [spot_minor_axis this_cell_spots.spot_minor_axis];
                    spot_ecc = [spot_ecc this_cell_spots.spot_ecc];
                    spot_area = [spot_area this_cell_spots.spot_area];
                    rel_orient = [rel_orient this_cell_spots.rel_orient];
                    rel_signal_int = [rel_signal_int this_cell_spots.rel_signal_int];
                    signal_int = [signal_int this_cell_spots.signal_int];
                    cell_length = [cell_length this_cell.length];
                    cluster_position = [cluster_position this_cell_spots.l/this_cell.length];
                end
            end
        end
    end
end
    
% convert units from pixel to um
spot_major_axis = spot_major_axis.*pixelsize;
spot_minor_axis = spot_minor_axis.*pixelsize;
spot_area = spot_area.*(pixelsize^2);
cell_length = cell_length.*pixelsize;

% map the cluster positions (relative to cell length to the interval
% [0,0.5]
cluster_position = arrayfun(@Mapping,cluster_position);

% calculate the mean and sd of the relative signal intensity for different
% cluster positions
nbins = 5; % number of bins to discretize the data

data = table(rel_signal_int', cluster_position', discretize(cluster_position,nbins)');
data.Properties.VariableNames = {'rel_signal_int' 'cluster_position' 'group'};

[G, ID] = findgroups(data.group);
mean_rel_signal_int = splitapply(@mean,data.rel_signal_int,G);
sd_rel_signal_int = splitapply(@std,data.rel_signal_int,G);

[~,edges] = discretize(cluster_position,nbins);
midpoints = mean([edges(1:length(edges)-1);edges(2:length(edges))],1).*100;

data2 = table(signal_int', cluster_position', discretize(cluster_position,nbins)');
data2.Properties.VariableNames = {'signal_int' 'cluster_position' 'group'};

[G2, ID2] = findgroups(data2.group);
mean_signal_int = splitapply(@mean,data2.signal_int,G2);
sd_signal_int = splitapply(@std,data2.signal_int,G2);

figure;
set(gcf,'Unit','Inches','position',[0,0,7,7],'PaperUnits', 'Inches', 'PaperSize', [7, 7])
subplot(2,2,1);
histogram(spot_major_axis, 'Normalization','probability')
hold on 
line([mean(spot_major_axis) mean(spot_major_axis)], [0 0.3],'Color','red','LineStyle','--')
text(0.05,0.85,['avg ' char(177) 'sd = ' newline num2str(round(mean(spot_major_axis),3)) char(177) num2str(round(std(spot_major_axis),3)) '\mum' newline '(N = ' num2str(length(spot_major_axis)) ')'],'Units','normalized')
xlabel('spot major axis [\mum]');
ylabel('probability');

subplot(2,2,2);
histogram(spot_minor_axis, 'Normalization','probability');
hold on 
line([mean(spot_minor_axis) mean(spot_minor_axis)], [0 0.4],'Color','red','LineStyle','--')
text(0.05,0.85,['avg ' char(177) ' sd = ' newline num2str(round(mean(spot_minor_axis),3)) char(177) num2str(round(std(spot_minor_axis),3)) '\mum' newline '(N = ' num2str(length(spot_minor_axis)) ')'],'Units','normalized')
xlabel('spot minor axis [\mum]');
ylabel('probability');

subplot(2,2,3);
histogram(spot_area, 'Normalization','probability');
hold on 
line([mean(spot_area) mean(spot_area)], [0 0.4],'Color','red','LineStyle','--')
text(0.05,0.85,['avg ' char(177) ' sd = ' newline num2str(round(mean(spot_area),3)) char(177) num2str(round(std(spot_area),3)) '\mum^2' newline '(N = ' num2str(length(spot_area)) ')'],'Units','normalized')
xlabel('spot area [\mum^2]');
ylabel('probability');

% save the figure as pdf
saveas(gcf,[save_path 'spot_statistics.pdf'])
close(gcf)

figure;
set(gcf,'Unit','Inches','position',[0,0,7,7],'PaperUnits', 'Inches', 'PaperSize', [7, 7])

subplot(2,2,1);
histogram(rel_signal_int.*100, 'Normalization','probability');
hold on 
line([mean(rel_signal_int)*100 mean(rel_signal_int)*100], [0 0.4],'Color','red','LineStyle','--')
text(0.05,0.85,['avg' char(177) ' sd = ' newline num2str(round(mean(rel_signal_int.*100),3)) char(177) num2str(round(std(rel_signal_int.*100),3)) newline '(N = ' num2str(length(rel_signal_int)) ')'],'Units','normalized')
xlabel('cluster signal intensity [% in cell]');
ylabel('probability');

subplot(2,2,2);
scatter(cell_length, rel_signal_int.*100);
xlabel('cell length [\mum]');
ylabel('cluster signal intensity [% in cell]');
axis([min(cell_length) max(cell_length) ...
    min(rel_signal_int)*100 max(rel_signal_int)*100])

subplot(2,2,3);
scatter(cluster_position.*100, rel_signal_int.*100);
xlabel('cluster position [% of cell length]');
ylabel('cluster signal intensity [% in cell]');
axis([min(cluster_position)*100 max(cluster_position)*100 ...
    min(rel_signal_int)*100 max(rel_signal_int)*100])
hold on 
errorbar(midpoints(ID), mean_rel_signal_int*100, sd_rel_signal_int*100);
xlim([0 50]);

% save the figure as pdf
saveas(gcf,[save_path 'spot_statistics_rel_signal_intensity.pdf'])
close(gcf)

figure;
set(gcf,'Unit','Inches','position',[0,0,7,7],'PaperUnits', 'Inches', 'PaperSize', [7, 7])

subplot(2,2,1);
H = histogram(signal_int, 'Normalization','probability');
hold on 
line([mean(signal_int) mean(signal_int)], [0 max(H.Values)],'Color','red','LineStyle','--')
text(0.05,0.85,['avg' char(177) ' sd = ' newline num2str(round(mean(signal_int),3)) char(177) num2str(round(std(signal_int),3)) newline '(N = ' num2str(length(signal_int)) ')'],'Units','normalized')
xlabel('cluster signal intensity');
ylabel('probability');

subplot(2,2,2);
scatter(cell_length, signal_int);
xlabel('cell length [\mum]');
ylabel('cluster signal intensity');
axis([min(cell_length) max(cell_length) ...
    min(signal_int) max(signal_int)])

subplot(2,2,3);
scatter(cluster_position.*100, signal_int);
xlabel('cluster position [% of cell length]');
ylabel('cluster signal intensity');
axis([min(cluster_position)*100 max(cluster_position)*100 ...
    min(signal_int) max(signal_int)])
hold on 
errorbar(midpoints(ID2), mean_signal_int, sd_signal_int);
xlim([0 50]);

% save the figure as pdf
saveas(gcf,[save_path 'spot_statistics_signal_intensity.pdf'])
close(gcf)

figure
set(gcf,'Unit','Inches','position',[0,0,7,7],'PaperUnits', 'Inches', 'PaperSize', [7, 7])
subplot(2,2,1)
scatter(spot_major_axis, spot_minor_axis)
xlabel('spot major axis [\mum]');
ylabel('spot minor axis [\mum]');
axis([min([spot_major_axis, spot_minor_axis]) max([spot_major_axis, spot_minor_axis]) ...
    min([spot_major_axis, spot_minor_axis]) max([spot_major_axis, spot_minor_axis])])

subplot(2,2,2)
histogram(spot_ecc,'Normalization','probability');
hold on
line([mean(spot_ecc) mean(spot_ecc)], [0 0.5],'Color','red','LineStyle','--')
text(0.05,0.85,['avg ' char(177) ' sd = ' newline num2str(round(mean(spot_ecc),3)) char(177) num2str(round(std(spot_ecc),3)) newline '(N = ' num2str(length(spot_ecc)) ')'],'Units','normalized')
xlabel('eccentricity')
ylabel('probability')

subplot(2,2,3)
histogram(rel_orient,'Normalization','probability')
hold on
line([mean(rel_orient) mean(rel_orient)], [0 0.5],'Color','red','LineStyle','--')
text(0.05,0.85,['avg ' char(177) ' sd = ' newline num2str(round(mean(rel_orient),3)) char(177) num2str(round(std(rel_orient),3)) newline '(N = ' num2str(length(rel_orient)) ')'],'Units','normalized')
xlabel(['orientation [' char(176) ']'])
ylabel('probability')

% save the figure as pdf
saveas(gcf,[save_path 'spot_statistics_eccentricity.pdf'])
close(gcf)
end