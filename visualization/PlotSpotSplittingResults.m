function [] = PlotSpotSplittingResults(mat_paths, save_path, pixelsize)
% PLOTSPOTSPLITTINGRESULTS: Plot the sizes of the clusters after splitting as
% well as the asymmetry in the size of the splitted clusters
%
% Copyright (c) 2021 Silke Bergeler
%
% Input: 
% - mat_paths: path to mat files
% - save_path: path to save figure
% - pixelsize: pixel size of images

dataAll = cell(1,length(mat_paths));
for j = 1:length(mat_paths)
    masks = load(mat_paths{j});
    nframes = length(masks.cellList.meshData);
    
    % get the infos of all daughter cells in the frame upon division (that's
    % when there is an entry 'ancestors' in the meshData)
    % only take the cells with one spot detected
    meshes_daughters = [];
    for frame = 1:nframes
        for ncell = 1:length(masks.cellList.meshData{frame})
            if ~isempty(masks.cellList.meshData{frame}{ncell}.ancestors) && ...
                    length(masks.cellList.meshData{frame}{ncell}.spots_matlab) == 1 && ...
                    ~all(structfun(@(x) isempty(x), masks.cellList.meshData{frame}{ncell}.spots_matlab))
                meshes_daughters = [meshes_daughters, masks.cellList.meshData{frame}{ncell}];
            end
        end
    end
    
    % take only the data where there are two daughters of the same ancestor
    ancestors = arrayfun(@(x) meshes_daughters(x).ancestors, 1:length(meshes_daughters));
    [count,val] = hist(ancestors,unique(ancestors));
    selected = val(count == 2);
    % only proceed if there is at least one
    if isempty(selected)
        break
    end
    selected_ind = ismember(ancestors, selected);
    meshes_daughters = meshes_daughters(1,selected_ind);
    
    % create a table with the area of the cluster in each daughter cell, the
    % min and the max value and the asymmetry value
    data = [];
    for i=1:length(selected) % all ancestor values that appear exactly twice
        data_tmp = [];
        for l=1:length(meshes_daughters)
            if meshes_daughters(l).ancestors == selected(i)
                data_tmp = [data_tmp, meshes_daughters(l).spots_matlab.spot_area*pixelsize^2];
            end
        end
        data(i,:) = [data_tmp(1),data_tmp(2)];
    end
    dataAll{j} = data;
end

% concatenate data
data = vertcat(dataAll{:});

if(isempty(data))
    disp('no data to plot')
    return
end

data(:,3) = min(data,[],2);
data(:,4) = max(data,[],2);
data(:,5) = arrayfun(@(x) (data(x,1)-data(x,2))/(data(x,1) + data(x,2)),1:size(data,1));

data = array2table([repmat(j,size(data,1),1) data],'VariableNames',{'dataset','area1','area2','min','max','asym'});

%%
figure;
set(gcf,'Unit','Inches','position',[0,0,5,5],'PaperUnits', 'Inches', 'PaperSize', [5, 5])

% plot a histogram of the asymmetry values 
subplot(2,1,1);
histogram(data.asym, 'Normalization','probability')
hold on 
line([mean(data.asym) mean(data.asym)], get(gca, 'ylim'),'Color','red','LineStyle','--')
text(0.05,0.85,['avg ' char(177) ' sd = ' newline num2str(round(mean(data.asym),3)) char(177) num2str(round(std(data.asym),3)) newline '(N = ', num2str(length(data.asym)),')'],'Units','normalized')
xlabel('asymmetry values');
ylabel('probability');

% scatter plot of the min and max size of the cluster
subplot(2,1,2);
scatter(data.min, data.max)
xlabel('area (small) [\mum^2]');
ylabel('area (large) [\mum^2]');

% save the figure as pdf
saveas(gcf,[save_path 'cluster_splitting.pdf'])
close(gcf)

end

