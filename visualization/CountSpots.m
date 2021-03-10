function [ ] = CountSpots(mat_paths, save_path)
% COUNTSPOTS: Count the number of spots per cell and plot a histogram with
% the data
%
% Copyright (c) 2021 Silke Bergeler
%
% Input:
% - mat_paths: paths to mat files
% - save_path: path to save figure

no_of_spots_Matlab = [];
count_Matlab = 0; % total number of cells (analyzed with Matlab)

for j = 1:length(mat_paths)
    load(mat_paths{j},'cellList')
    for nframe = 1:size(cellList.meshData,2)
        for ncell = 1:size(cellList.meshData{nframe},2)
            if isfield(cellList.meshData{nframe}{ncell},'spots_matlab')
                count_Matlab = count_Matlab + 1;
                if isempty(cellList.meshData{nframe}{ncell}.spots_matlab(1).l)
                    no_of_spots_Matlab = [no_of_spots_Matlab 0];
                else
                    no_of_spots_Matlab = [no_of_spots_Matlab length(cellList.meshData{nframe}{ncell}.spots_matlab)];
                end
            end
        end
    end
end

figure
set(gcf,'Unit','Inches','position',[0,0,4,3],'PaperUnits', 'Inches', 'PaperSize', [4, 3])
[N, edges] = histcounts(no_of_spots_Matlab);
midpoints = mean([edges(2:end)',edges(1:end-1)'],2);
histogram(no_of_spots_Matlab,'Normalization','probability');
for i=1:length(N)
    text(midpoints(i)-0.2,N(i)/sum(N)+0.05,num2str(round(N(i)/sum(N),2)));
end
xlabel('#spots')
ylabel('probability')
text(0.7,0.8,['N = ' num2str(count_Matlab)],'Units','normalized')

% save the figure as pdf
saveas(gcf,[save_path 'number_of_spots_per_cell.pdf'])
close(gcf)
end