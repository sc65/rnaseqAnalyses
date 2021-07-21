%%
masterFolder = '/Users/sapnachhabra/Desktop/CellTrackercd/Experiments/rnaSeqExperiments/rnaSeq_datasets_fpkmTable';
cellsFile =[masterFolder filesep  'monkeyPostImplantationMa2019' filesep 'singleCells.mat'];
genesFile = [masterFolder filesep 'humanEmbryoAmnionDataD14/forReport' filesep 'forReport_genes.mat'];
%%
load(cellsFile); load(genesFile);

%% enter information about lineages, days and gene lists to plot
samples.lineages = {'CTB', 'STB', 'EVT', 'AME'};
samples.days = {'D12', 'D14'};
cellIdsToRemove = {'D12A3S10'};
genesToPlot = {'CGA', 'CGB1'};

%% enter information about which type of graph to plot
boxPlot = 0;
heatmapPlot = 0;
scatterPlot.toPlot = 1;
scatterPlot.xAxisGene = 'CGA';
scatterPlot.ameNamesOnPlot = 1;

%% run
%geneCellTable = geneCellTableNew_human_ensemblIds;
plotGeneExpression (geneCellTable, ldm, colorCode, samples, genesToPlot, ...
    boxPlot, heatmapPlot, scatterPlot, cellIdsToRemove)

%%
saveInPath = ['/Users/sapnachhabra/Desktop/CellTrackercd/Experiments/rnaSeqExperiments/rnaSeq_datasets_fpkmTable/'...
    'humanEmbryoAmnionDataD14/forReport/revision1/withNewGeneCellTable/'...
    'ctbEvtStbAme_scatterPlot'];
saveAllOpenFigures(saveInPath);
%%
%% ------------------------------------------------------------------------------------------------------------
%% ------------------------------------------------------------------------------------------------------------
%% ------------------------------------------------------------------------------------------------------------

%%
function plotGeneExpression (geneCellTable, ldm, colorCode, samples, geneList1, ...
    boxPlot, heatmapPlot,  scatterPlot, cellIdsToRemove)

% get required genes
[~, ~, idx2] = intersect(geneList1, table2cell(geneCellTable(:,2)), 'stable');
readCountMatrix = log2(table2array(geneCellTable(idx2,3:end))+1);
cellLabels = geneCellTable.Properties.VariableNames(3:end);

% get the required samples
idx1 = ismember(ldm.lineageLabels, samples.lineages);
idx2 = ismember(ldm.dayLabels, samples.days);
idx3 = ismember(ldm.cellLabels, cellLabels);
idx4 = idx1&idx2&idx3;

samples.lineageLabels = strcat(ldm.lineageLabels(idx4), ':', ldm.dayLabels(idx4)); % sorted by days
samples.fpkm = readCountMatrix(:, idx4);
samples.cellLabels = cellLabels(idx4);

%%
% remove unwanted cells
if ~isempty(cellIdsToRemove)
    [~, idx1,~] = intersect(samples.cellLabels, cellIdsToRemove);
    samples.cellLabels(idx1) = [];
    samples.fpkm(:,idx1) = [];
    samples.lineageLabels(:,idx1) = [];
end

%%
% sort according to lineages;
% lineage1_day1, lineage1_day2; lineage2_day1, lineage2_day2

[~, idx1] = natsortfiles(samples.lineageLabels); % makes sure D10 comes after D6.
samples.lineageLabels = samples.lineageLabels(idx1);
samples.fpkm = samples.fpkm(:,idx1);
samples.cellLabels = samples.cellLabels(idx1);
%%
% shuffle again to get the lineage order same as samples.lineages
sampleOrder1 = cell(1, numel(samples.lineages));
for ii = 1:numel(samples.lineages)
    sampleOrder1{ii} = find(contains(samples.lineageLabels, samples.lineages{ii}));
end
samplesOrder = cat(2, sampleOrder1{:});
samples.lineageLabels = samples.lineageLabels(samplesOrder);
samples.fpkm = samples.fpkm(:,samplesOrder);
samples.cellLabels = samples.cellLabels(:,samplesOrder);
%%
% plot
if boxPlot == 1
    % color boxplots of each lineage separately
    lineageLabels = unique(samples.lineageLabels, 'stable');
    lineageColors = zeros(numel(lineageLabels),3);
    for ii = 1:numel(lineageLabels)
        lIdx = ismember(ldm.lineages, strtok(lineageLabels{ii}, ':'));
        lineageColors(ii,:) = colorCode.lineages(lIdx,:);
    end
    %%
    %figure; nRows = 1; nColumns = 1;
    for ii = 1:size(samples.fpkm,1)
        %subplot(nRows, nColumns, ii); hold on;
        figure;
        h1 = boxplot(samples.fpkm(ii,:), samples.lineageLabels);
        set(h1,{'linew'},{2});
        
        h2 = findobj(gca,'Tag','Box');
        for jj = 1:length(h2)
            patch(get(h2(jj),'XData'),get(h2(jj),'YData'),lineageColors(end-jj+1,:),'FaceAlpha',.3);
        end
        
        ax = gca;
        if sum(contains(ax.XTickLabels, 'Intermediate')) > 0 % shorten the label
            ax.XTickLabels = strrep(ax.XTickLabels, 'Intermediate', 'Inter');
        end
        ax.FontSize = 20; ax.FontWeight = 'bold';
        title([geneList1(ii)]); ylim([0 max(samples.fpkm(ii,:))+0.5]);
    end
end
%%
if scatterPlot.toPlot == 1
    
    %%
    % save marker types and marker colors for each lineage
    lineageLabels = samples.lineages;
    lineageColors = zeros(numel(lineageLabels),3);
    lineagePlotSymbols = repmat({'o'}, 1, numel(lineageLabels));
    lineagePlotSymbols(ismember(lineageLabels, {'AME'})) = {'^'};
    lineagePlotAlpha = repmat(0.6, 1, numel(lineageLabels));
    lineagePlotAlpha(ismember(lineageLabels, {'AME'})) = 1;
    
    for jj = 1:numel(lineageLabels)
        lIdx = ismember(ldm.lineages, strtok(lineageLabels{jj}, ':'));
        lineageColors(jj,:) = colorCode.lineages(lIdx,:);
    end
    %%
    % % sort samples according to fpkm values of marker genes on xAxis
    xAxisGene_idx = ismember(geneList1, scatterPlot.xAxisGene);
    xAxisGene_fpkm = samples.fpkm(xAxisGene_idx,:);
    [~, idx1] = sort(xAxisGene_fpkm);
    fpkm1 = samples.fpkm(:,idx1);
    lineageLabels1 = samples.lineageLabels(idx1);
    cellLabels1 = samples.cellLabels(idx1);
    %%
    % plot
    for jj = 1:numel(geneList1)
        figure; hold on;
        for kk = 1:numel(lineageLabels)
            idx1 = contains(lineageLabels1, lineageLabels{kk});
            s = scatter(fpkm1(xAxisGene_idx,idx1), fpkm1(jj,idx1), 120, 'filled', lineagePlotSymbols{kk}, ...
                'MarkerEdgeColor', lineageColors(kk,:),'MarkerFaceColor', lineageColors(kk,:),...
                'LineWidth',1.5);
            alpha(s,lineagePlotAlpha(kk));
        end
        
        if scatterPlot.ameNamesOnPlot == 1
            idx1 = contains(lineageLabels1, 'AME');
            amnionCells = cellLabels1(idx1);
            text(fpkm1(xAxisGene_idx,idx1)+0.25, fpkm1(jj,idx1)+0.1, strrep(amnionCells, '_', ':'), 'Color', 'k', ...
                'FontSize', 18);
        end
        
        %title([samples.days{1} samples.days{2}]);
        ylabel(geneList1{jj}); xlabel(scatterPlot.xAxisGene); legend(lineageLabels);
        xlim([-0.25 max(fpkm1(xAxisGene_idx,:))+0.5]); ylim([-0.25 max(fpkm1(jj,:))+0.5])
        ax = gca; ax.FontSize = 30; ax.FontWeight = 'bold';
    end
    
    
end


if heatmapPlot == 1
    figure;
    h = heatmap(zscore(samples.fpkm, 0, 2));
    colormap('jet');
    caxis([-2 3]);
    
    h.YDisplayLabels = geneList1;
    %h.XDisplayLabels = nan(size(samples.fpkm,2),1);
    %ax = gca;
    %ax.FontSize = 20;
end




end


