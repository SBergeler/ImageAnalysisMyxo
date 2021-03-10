function [] = PlotCellStatistics(celldataAll, save_path)
% PLOTCELLSTATISTICS: this function plots the cell lengths and widths as a histogram
% made for stacks of independent images and not time lapse images
%
% Copyright (c) 2021 Silke Bergeler
%
% Input:
% - celldataAll: cell array with data about cells and spots
% - save_path: path to save the figure

lengths = table2array(celldataAll(:,'length'));
widths = table2array(celldataAll(:,'width'));

figure;
set(gcf,'Unit','Inches','position',[0,0,5,5],'PaperUnits', 'Inches', 'PaperSize', [5, 5])

subplot(2,1,1);
histogram(lengths, 'Normalization','probability')
hold on 
line([mean(lengths) mean(lengths)], get(gca, 'ylim'),'Color','red','LineStyle','--')
text(0.05,0.85,['avg ' char(177) ' sd = ' newline num2str(round(mean(lengths),3)) char(177) ...
    num2str(round(std(lengths),3)) '\mum' newline '(N = ', num2str(length(lengths)),')'],'Units','normalized')
xlabel('cell length [\mum]');
ylabel('probability');

subplot(2,1,2);
histogram(widths, 'Normalization','probability');
hold on 
line([mean(widths) mean(widths)], get(gca, 'ylim'),'Color','red','LineStyle','--')
text(0.05,0.85,['avg ' char(177) ' sd = ' newline num2str(round(mean(widths),3)) char(177) ...
    num2str(round(std(widths),3)) '\mum' newline '(N = ',num2str(length(widths)),')'],'Units','normalized')
xlabel('cell width [\mum]');
ylabel('probability');

% save the figure as pdf
saveas(gcf,[save_path 'cell_statistics.pdf'])
close(gcf)
end