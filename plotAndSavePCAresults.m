%% --------------------------------------------------------------------------------------------------------
clearvars; close all;
%%
cellsFile = 'forReport_singleCells.mat';
genesFile = 'forReport_genes.mat';
%
load(cellsFile); load(genesFile);
%%
savePCAinfo = 0; % change to zero if you only need to plot  previously saved results
%% --------------------------------------------------------------------------------------------------------

%% --------------------------------------------------------------------------------------------------------
%% -> Specify the lineages and days for analyses

% 1. find most variable genes between a given set of lineages.
% 2. retain cellIdx correponding to the two lineages
% 3. retain genes - genes1, expressed (fpkm>1) in at least half of the cells of either cell lineage
% calculate cv of genes1
% perform pca on the entire dataset by selecting a portion of genes1 (based
% on cv)

%%
% samples_lineages1 = {'AME', 'CTB', 'STB', 'EVT'};
% samples_daysIdx = [7:10 12 14];
% samples_days1 = strcat('D', strsplit(int2str(samples_daysIdx), ' '));

%
samples_lineages1 = {'AME', 'CTB', 'STB', 'EVT'};
samples_days1 = {'D12', 'D14'};
cellIdToRemove = {'D12A3S10'}; % cell to remove from analyses

PCAtSNE_struct_field = 'ameCtb'; % a prefix appropriate for samples_lineages1

proteinCodingGenesRestrict = 1; % 1 = restrict analyses to protein coding genes
%% -> (1a)Get fpkm values corresponding to expressed genes
% genes with fpkm >1 in atleast half of the cells of every lineage
%%
plotExpressedGeneTypesPie = 0; %a pie chart of expressed gene types in all lineages
plotExpressedGeneTypesBar = 0; %a bar chart with fraction of protein coding and pseudogenes expressed
geneCellTableNew = retainExpressedGenes(ldm, samples_lineages1, samples_days1,  geneCellTable, ...
    plotExpressedGeneTypesPie, plotExpressedGeneTypesBar, ensemblGenes, geneTypes, colorCodeG);

%% -> (1b)(Optional) gene types expressed in single cells of a given lineage
%%
samples_lineages11 = {'AME', 'Intermediate'}; % >=1 lineages
samples_days11 = 'D12'; % a single day

pieChart = 1;
barChart = 1;

for ii = samples_lineages11
    plotExpressedGenesTypesSingleCells(ldm, ii, samples_days11, geneCellTable, ...
        ensemblGenes, geneTypes, colorCodeG, pieChart, barChart)
end
%% -> (2)keep only protein coding genes
if proteinCodingGenesRestrict == 1
    proteinCodingGenes = unique(ensemblGenes.name(ismember(ensemblGenes.type, 'protein_coding')));
    [~, idx1, ~] = intersect(table2cell(geneCellTableNew(:,2)), proteinCodingGenes);
    geneCellTableNew (setxor(1:size(geneCellTableNew,1), idx1),:) = [];
end
%% -> (3)exclude specific cells
if ~isempty(cellIdToRemove)
    [~,  ~, idx2] = intersect(cellIdToRemove, geneCellTableNew.Properties.VariableNames);
    geneCellTableNew(:,idx2) = [];
end
%%
%% ---------------------------------------------------------------------------------------------------------
%% ---------------------------Skip to section6 if analysing previously saved data---------------------------
%% ---------------------------------------------------------------------------------------------------------

%% -> (4) Make the cv vs mean expression graph
% find cv and mean of log transformed (to decrease range) values

readCounts_fpkm = table2array(geneCellTableNew(:,3:end));

readCounts_fpkmNew = log2(readCounts_fpkm+1);
%fpkmNew = zscore(fpkmNew, 0, 'all');
readCounts_genesMean = mean(readCounts_fpkmNew,2);
readCounts_genesSTD = std(readCounts_fpkmNew, 0,2);
readCounts_genesCV = readCounts_genesSTD./readCounts_genesMean;
%
figure;
subplot(1,2,1); histogram(readCounts_genesMean); xlabel('Mean');
ax = gca; ax.FontSize = 14; ax.FontWeight = 'bold';

subplot(1,2,2); histogram(readCounts_genesCV);xlabel('CV');
ax = gca; ax.FontSize = 14; ax.FontWeight = 'bold';

figure; plot(readCounts_genesMean, readCounts_genesCV, 'k.'); xlabel('Mean'); ylabel('CV');
ax = gca; ax.FontSize = 14; ax.FontWeight = 'bold';
%
% figure; plot(genesMean, genesCV.^2, 'k.'); xlabel('Mean'); ylabel('CV^2');
% ax = gca; ax.FontSize = 14; ax.FontWeight = 'bold';
%%
%%  ------------------------------------------------------------------------------
%%  ------------------------------------------------------------------------------
%%  ------------------------------------------------------------------------------
%% -> (5) Two different options to perform PCA/tSNE
%% (A) Save PCA/tSNE results in a structure for later analyses (GO, Gene PC contribution)
%% (B) Perform quick PCA/tSNE for a given cell, gene list. Plot results real time
%% -> save PCA output in a structure : variableGenesPCA

%%  ------------------------------------------------------------------------------
%%  ------------------------------------------------------------------------------
%% (Option A) -> Save PCA/tSNE results in a structure for later analyses
%%  ------------------------------------------------------------------------------
cvThreshold = [0.5]; % variability cutoff.

savePCAinfo = 1;
savetSNEinfo = 0;
daysPlot = 0;
ameNamesOnPlot = 1;
%%
if savePCAinfo == 1
    PCA_struct = struct;
end

if savetSNEinfo == 1
    tSNE_struct = struct;
end
%%
expressedGenes.geneIds = table2cell(geneCellTableNew(:,1));
expressedGenes.geneNames = table2cell(geneCellTableNew(:,2));

counter = 1;

for ii = cvThreshold
    
    readCounts_fpkmNew1 = readCounts_fpkmNew(readCounts_genesCV>ii,:);
    genes1 = expressedGenes.geneNames(readCounts_genesCV>ii);
    geneIds1 = expressedGenes.geneIds(readCounts_genesCV>ii);
    
    nGenes = numel(genes1);
    numel(genes1) - numel(unique(genes1))
    
    
    if savePCAinfo == 1
        PCAplot = 1; tSNEplot = 0;
        [pc_coeff1, pc_score1, ~, ~, pc_varExplained1, ~] = pca(readCounts_fpkmNew1');
        pc_varExplained1(1:5)
        
        pcTsne.score = pc_score1;
        pcTsne.varExplained =  pc_varExplained1(1:2);
        
        [~, idx1] = sort(abs(pc_coeff1(:,1)), 'desc'); % arranged according to contributions to pc1 (high to low)
        PCA_struct = savePCAstructureFieldNames(PCA_struct, {counter}, PCAtSNE_struct_field, ii, genes1(idx1), geneIds1(idx1), ...
            pc_coeff1(idx1,1), pc_coeff1(idx1,2), pc_varExplained1, pc_score1);
        
        if PCAplot == 1
            plotScore(samples_lineages1, samples_days1, ldm, colorCode, geneCellTableNew, ...
                pcTsne,  nGenes, cvThreshold(counter), PCAplot, tSNEplot, daysPlot, ameNamesOnPlot);
        end
    end
    
    if savetSNEinfo == 1
        tSNEplot = 1; PCAplot = 0;
        Y = tsne(readCounts_fpkmNew1','Algorithm','barneshut','NumPCAComponents',50);
        pcTsne.score = Y;
        
        tSNE_struct = savetSNEstructureFieldNames(tSNE_struct, {counter}, PCAtSNE_struct_field, ii,...
            genes1, geneIds1, pcTsne.score);
        
        if tSNEplot == 1
            plotScore(samples_lineages1, samples_days1, ldm, colorCode, geneCellTableNew, ...
                pcTsne, nGenes, cvThreshold(counter), PCAplot, tSNEplot, daysPlot, ameNamesOnPlot);
        end
    end
    
    
    counter = counter+1;
end

%%
saveAllOpenFigures(['/Users/sapnachhabra/Desktop/CellTrackercd/Experiments/rnaSeqExperiments/rnaSeq_datasets_fpkmTable/'...
    'humanEmbryoAmnionDataD14/forReport/revision1/withNewGeneCellTable/pcatsne/'...
    'proteinCodingGenes_allLineages/ameNamesOnPlot'])

%% -------------------------------------------------------------------------------------------------------
%% -------------------------------------------------------------------------------------------------------
%% -------------------------------------------------------------------------------------------------------
%% -------------------------------------------------------------------------------------------------------
%% (OptionB) -> perform and plot pca for a given gene list
%% -------------------------------------------------------------------------------------------------------
geneList1 = geneLists.humanEmbryoLineage_revised;
samples_lineages1 = ldm.lineages([2 3 7]);
samples_days1 = ldm.days;
ameNamesOnPlot = 0;  cvThreshold = [];
PCAplot = 1; daysPlot = 0; tSNEplot = 0;
%%
[~, ~, idx2] = intersect(geneList1, table2cell(geneCellTableNew(:,2)), 'stable');

readCounts_fpkmNew1 = log2(table2array(geneCellTableNew(idx2,3:end))+1);
[pc_coeff1, pc_score1, ~, ~, pc_varExplained1, ~] = pca(readCounts_fpkmNew1');
pc_varExplained1(1:5)

pcTsne.score = pc_score1; pcTsne.varExplained = pc_varExplained1(1:2);
nGenes = size(readCounts_fpkmNew1,1);

plotScore(samples_lineages1, samples_days1, ldm, colorCode, geneCellTable, ...
    pcTsne,nGenes, cvThreshold, PCAplot, tSNEplot, daysPlot, ameNamesOnPlot);


%% -------------------------------------------------------------------------------------------------------
%% -------------------------------------------------------------------------------------------------------
%% -------------------------------------------------------------------------------------------------------
%% -> 6) Analyse previously saved PCA/tSNE data
%% -------------------------------------------------------------------------------------------------------
%% -> (6a) plot PCA
%% -------------------------------------------------------------------------------------------------------
samples_lineages1 = ldm.lineages;
samples_days2 = {'D12', 'D14'}; % [days to be plotted separately]; [Note: PCA is not computed twice]

cvThreshold = [0.5 1 1.5 2];
cvThreshold_indices = [1:3]; % cvThreshold goes from [0.5:0.5:2].

PCAtSNE_struct = variableProteinCodingGenesPCA;
PCAtSNE_struct_field = 'allLineages';

tSNEplot = 0; PCAplot = 1;%
%[if PCAtSNE_struct corresponds to PCA, set PCAplot = 1; tSNEplot = 0;
% else PCAplot = 0; tSNEplot = 1]

daysPlot = 1;
ameNamesOnPlot = 0;

for ii = 1:numel(cvThreshold_indices)
    
    pc_genes1 = getfield(PCAtSNE_struct(cvThreshold_indices(ii)), PCAtSNE_struct_field, 'genes');
    nGenes = numel(pc_genes1);
    
    if PCAplot == 1
        pc_score1 = getfield(PCAtSNE_struct(cvThreshold_indices(ii)), PCAtSNE_struct_field, 'pca', 'score');
        pc_varExplained1 = getfield(PCAtSNE_struct(cvThreshold_indices(ii)), PCAtSNE_struct_field, 'pca', 'varianceExplained');
        pcTsne.varExplained =  pc_varExplained1(1:2);
    else
        pc_score1 = getfield(PCAtSNE_struct(cvThreshold_indices(ii)), PCAtSNE_struct_field, 'tsne', 'score');
        pcTsne.varExplained =  [];
    end
    
    pcTsne.score = pc_score1;
   
    plotScore(samples_lineages1, samples_days2, ldm, colorCode, geneCellTableNew, ...
        pcTsne,nGenes, cvThreshold(ii), PCAplot, tSNEplot, daysPlot, ameNamesOnPlot); 
end

%%
saveAllOpenFigures(['/Users/sapnachhabra/Desktop/CellTrackercd/Experiments/rnaSeqExperiments/rnaSeq_datasets_fpkmTable/'...
    'humanEmbryoAmnionDataD14/forReport/revision1/withNewGeneCellTable/pcatsne/'...
    'proteinCodingGenes_allLineages/withoutAmeNames'])

%% -------------------------------------------------------------------------------------------------------
%% -> (6b_1) plot gene contributions of various gene categories to pc1 and pc2
%% -------------------------------------------------------------------------------------------------------

cvThreshold = [0.5:0.5:2.0];
cvThreshold_indices = [1:4];
PCAtSNE_struct = variableGenesPCA;
PCAtSNE_struct_field = 'allLineages';

for ii = 1:numel(cvThreshold)
    clearvars pc_coeff1 pc_genes1
    
    pc_coeff1(:,1) = getfield(PCAtSNE_struct(cvThreshold_indices(ii)), PCAtSNE_struct_field, 'pca', 'geneCoeffpc1');
    pc_coeff1(:,2) = getfield(PCAtSNE_struct(cvThreshold_indices(ii)), PCAtSNE_struct_field, 'pca', 'geneCoeffpc2');
    pc_genes1 = getfield(PCAtSNE_struct(cvThreshold_indices(ii)), PCAtSNE_struct_field, 'genes');
    
    [pcGeneTypes{ii}, pcCoeffCum{ii}] = makePcBiplotGeneTypes(cvThreshold(ii), pc_genes1, pc_coeff1, ...
        ensemblGenes, geneTypes, colorCodeG);
end

%% -> (6b_2) plot gene contributions to PCA as bar plots

toKeep = {'pseudogene', 'protein_coding'}; %(combine all pseudogene categories)
pcCoeffCumAll = zeros(2,2*length(pcGeneTypes));
%pattern = [[cv1_pc1Coeff_pseudoGenes, cv1_pc2coeff_pseudoGenes, cv2_pc1Coeff_pseudoGenes, cv2_pc2coeff_pseudoGenes;
%cv1_pc1Coeff_proteinCoding, cv1_pc2coeff_proteinCoding, cv2_pc1Coeff_proteinCoding, cv2_pc2coeff_proteinCoding=]
% repeats in

for ii = 1:numel(cvThreshold)
    pcCoeffCum{ii} = pcCoeffCum{ii}./sum(pcCoeffCum{ii});
    for jj = 1:numel(toKeep)
        idx1 = contains(pcGeneTypes{ii}, toKeep{jj});
        pcCoeffCumAll(jj,(ii-1)*2+1:(ii-1)*2+2) = sum(pcCoeffCum{ii}(idx1,:),1);
        
    end
end
%
pcCoeff_pc1 = pcCoeffCumAll(:,[1 3 5]); % only plot the first three cv's [0.5 1 1.5].
pcCoeff_pc2 = pcCoeffCumAll(:,[2 4 6]);

toPlot = {pcCoeff_pc1, pcCoeff_pc2};
for ii = 1:2
    figure;
    h = bar(toPlot{ii}');
    ax = gca; legend(strrep(toKeep, '_', ' '));
    ax.XTickLabel = strsplit(num2str([0.5 1 1.5]), ' ');
    ax.FontSize = 24; ax.FontWeight = 'bold';
    
end
%%
saveAllOpenFigures(['/Users/sapnachhabra/Desktop/CellTrackercd/Experiments/rnaSeqExperiments/rnaSeq_datasets_fpkmTable/'...
    'humanEmbryoAmnionDataD14/forReport/revision1/withNewGeneCellTable/pcatsne/'...
    'allGenes_allLineages/pcaContribution'])


%% -------------------------------------------------------------------------------------------------------
%% -> (6c_1) plot gene contribution to PC1
%% -------------------------------------------------------------------------------------------------------

PCAtSNE_struct = variableGenesPCA;
PCAtSNE_struct_field = 'allLineages';

cvThreshold = [0.5:0.5:2.0];
cvThreshold_indices = [1:4];

figure;
for ii = 1:numel(cvThreshold)
    subplot(2,2, ii);
    histogram([getfield(PCAtSNE_struct(cvThreshold_indices(ii)), PCAtSNE_struct_field, 'pca', 'geneCoeffpc1')], 'BinWidth', 0.01);
    title(['cv' num2str(cvThreshold(ii))]);
end

%% -------------------------------------------------------------------------------------------------------
%% -> (6c_2) save genelists for GO analyses
%% -------------------------------------------------------------------------------------------------------
% save  genes that contribute to max variability under
% different PCA conditions

saveInPath = ['/Users/sapnachhabra/Desktop/CellTrackercd/Experiments/rnaSeqExperiments/rnaSeq_datasets_fpkmTable/'...
    'humanEmbryoAmnionDataD14/forReport/revision1/withNewGeneCellTable/pcatsne/'...
    'allGenes_allLineages/pcaGeneLists'];
mkdir(saveInPath);

coeffThreshold = [0.8, 0.5, 0.2, 0];
PCAtSNE_struct = variableGenesPCA;
PCAtSNE_struct_field = 'allLineages';

cvThreshold = [0.5:0.5:2.0];
cvThreshold_indices = [1:4];

for ii = 1:numel(coeffThreshold)
    
    saveInPath1 = [saveInPath filesep 'variableGenes_names' filesep strrep(num2str(coeffThreshold(ii)), '.', '_')];
    mkdir(saveInPath1);
    
    saveInPath2 = [saveInPath filesep 'variableGenes_ids' filesep strrep(num2str(coeffThreshold(ii)), '.', '_')];
    mkdir(saveInPath2);
    
    
    for jj = 1:numel(cvThreshold)
        
        pc.cv1 = getfield(PCAtSNE_struct(cvThreshold_indices(jj)), PCAtSNE_struct_field, 'cvThreshold');
        pc.coeff1 = getfield(PCAtSNE_struct(cvThreshold_indices(jj)), PCAtSNE_struct_field, 'pca', 'geneCoeffpc1');
        pc.geneNames1 = getfield(PCAtSNE_struct(cvThreshold_indices(jj)), PCAtSNE_struct_field, 'genes');
        pc.geneIds1 = getfield(PCAtSNE_struct(cvThreshold_indices(jj)), PCAtSNE_struct_field, 'genesID');
        
        pc.genesNew = table(sort(pc.geneNames1(pc.coeff1>coeffThreshold(ii)*max(pc.coeff1))));
        pc.geneIdsNew = table(sort(pc.geneIds1(pc.coeff1>coeffThreshold(ii)*max(pc.coeff1))));
        
        fileName1 = [saveInPath1 filesep 'variableGenesCv' strrep(num2str(pc.cv1), '.', '_') '.txt'];
        fileName2 = [saveInPath2 filesep 'variableGenesCv' strrep(num2str(pc.cv1), '.', '_') '.txt'];
        
        writetable(pc.genesNew, fileName1, 'WriteVariableNames', 0);
        writetable(pc.geneIdsNew, fileName2, 'WriteVariableNames', 0);
    end
    
end

%% -----------------------------------------------------------------------------------------------------------------
%% -----------------------------------------------------------------------------------------------------------------
%% -----------------------------------------------------------------------------------------------------------------
%% -----------------------------------------------------------------------------------------------------------------
%%
function struct1 = savePCAstructureFieldNames(struct1, idx1, field1, cv1, genes1, genesId1, ...
    pc1_coeff, pc2_coeff, pc_varExplained, pc_score)

struct1 = setfield(struct1, idx1, field1, 'cvThreshold', cv1);
struct1 = setfield(struct1, idx1, field1, 'genes', genes1);
struct1 = setfield(struct1, idx1, field1, 'genesID', genesId1);
struct1 = setfield(struct1, idx1, field1, 'pca', 'geneCoeffpc1', pc1_coeff);
struct1 = setfield(struct1, idx1, field1, 'pca', 'geneCoeffpc2', pc2_coeff);
struct1 = setfield(struct1, idx1, field1, 'pca', 'varianceExplained', pc_varExplained);
struct1 = setfield(struct1, idx1, field1, 'pca', 'score', pc_score);

end
%% -----------------------------------------------------------------------------------------------------------------

function struct1 = savetSNEstructureFieldNames(struct1, idx1, field1, cv1, genes1, genesId1,tSNE_score)

struct1 = setfield(struct1, idx1, field1, 'cvThreshold', cv1);
struct1 = setfield(struct1, idx1, field1, 'genes', genes1);
struct1 = setfield(struct1, idx1, field1, 'genesID', genesId1);
struct1 = setfield(struct1, idx1, field1, 'tsne', 'score', tSNE_score);

end
%%
%% -----------------------------------------------------------------------------------------------------------------
%% -----------------------------------------------------------------------------------------------------------------
%%
function plotScore(lineages1, days2, ldm, colorCode, geneTableNew, ...
    pcTsne, nGenes, cv, PCAplot, tSNEplot, daysPlot, ameNamesOnPlot)
%% plot cell types in pc space.

if ~exist('daysPlot', 'var')
    daysPlot = 0;
end

if ~exist('ameNamesOnPlot', 'var')
    ameNamesOnPlot = 0;
end

if PCAplot == 1
    axisLabels = {['PC1 (' num2str(pcTsne.varExplained(1)) '%)'],...
        ['PC2 (' num2str(pcTsne.varExplained(2)) '%)']}; % [x,y];
elseif tSNEplot == 1
    axisLabels = {'tSNE1', 'tSNE2'};
end

if isempty(cv)
    title_text = ['genes' num2str(nGenes)];
else
    title_text = ['cv:' num2str(cv) 'genes' num2str(nGenes)];
end

%% ---------------------------------------------------------------------------
[~, lIdx1, ~] = intersect(ldm.lineages, lineages1);
xLimits = [round(min(pcTsne.score(:,1))-5) round(max(pcTsne.score(:,1))+5)];
yLimits = [round(min(pcTsne.score(:,2))-5) round(max(pcTsne.score(:,2))+5)];

daysPlot_rows = round(numel(days2)/2);
daysPlot_columns = round(numel(days2)/daysPlot_rows);

if daysPlot_rows*daysPlot_columns<numel(days2)
    daysPlot_columns = daysPlot_columns+1;
end

if daysPlot == 1
    figure; hold on;
     for jj = 1:numel(days2)
        subplot(daysPlot_rows,daysPlot_columns,jj); hold on;
        title([days2{jj} title_text]);
        
        [~, dIdx,~] = intersect(ldm.days, days2{jj});
        lIdx2 = find(any(ldm.cellCounts(:,dIdx),2));
        lIdx = intersect(lIdx1, lIdx2);
        if sum(ismember(lIdx,1))>0
            lineageOrder = setxor(lIdx',1);
            amnionPlot = 1;
        else
            lineageOrder = lIdx';
            amnionPlot = 0;
            
        end
        
        for kk = lineageOrder % plot amnion at the end
            cIdx = ldm.cellIdx{kk, dIdx};
            [~, idx1, ~] = intersect(geneTableNew.Properties.VariableNames(3:end), cIdx);
            plot(pcTsne.score(idx1,1), pcTsne.score(idx1,2), '.', 'Color', colorCode.lineages(kk,:), 'MarkerSize', 20);
            xlim(xLimits); ylim(yLimits);
        end
        
        if (amnionPlot) % amnion
            cIdx = ldm.cellIdx{1, dIdx};
            [~, idx1, ~] = intersect(geneTableNew.Properties.VariableNames(3:end), cIdx);
            plot(pcTsne.score(idx1,1), pcTsne.score(idx1,2), '^', 'Color', colorCode.lineages(1,:), ...
                'MarkerFaceColor', colorCode.lineages(1,:), 'MarkerSize', 18);
            lineageOrder = [lineageOrder,1];
            xlim(xLimits); ylim(yLimits);
        end
        %%-------legendIdx------%%%%%%%%%%%%%
        %legend(ldm.lineages(lineageOrder));
        xlabel(axisLabels{1}); ylabel(axisLabels{2});
        ax = gca; ax.FontSize = 24; ax.FontWeight = 'bold';
        
    end
end

%%
if sum(ismember(lIdx1,1))>0
    lineageOrder = setxor(lIdx1',1);
    amnionPlot = 1;
else
    lineageOrder = lIdx1';
    amnionPlot = 0;
end

figure; hold on;
title(title_text);
for jj = lineageOrder
    cIdx = cat(1, ldm.cellIdx{jj,:});
    [~, idx1, ~] = intersect(geneTableNew.Properties.VariableNames(3:end), cIdx);
    plot(pcTsne.score(idx1,1), pcTsne.score(idx1,2), '.', 'Color', colorCode.lineages(jj,:),...
        'MarkerSize', 25);
    xlim(xLimits); ylim(yLimits);
end

if amnionPlot == 1
    cIdx = cat(1, ldm.cellIdx{1,:});
    [amnionCells, idx1, ~] = intersect(geneTableNew.Properties.VariableNames(3:end), cIdx);
    plot(pcTsne.score(idx1,1), pcTsne.score(idx1,2), '^', 'Color', colorCode.lineages(1,:), ...
        'MarkerFaceColor', colorCode.lineages(1,:), 'MarkerSize', 15);
    lineageOrder = [lineageOrder,1];
    
    if ameNamesOnPlot == 1
        text(pcTsne.score(idx1,1)+3, pcTsne.score(idx1,2)+1, amnionCells, ...
            'Color', colorCode.lineages(1,:), ...
            'FontSize', 12) %'FontWeight', 'bold');
    end
    
end
%legend(ldm.lineages(lineageOrder));  xlim(xLimits); ylim(yLimits);

xlabel(axisLabels{1}); ylabel(axisLabels{2});
ax = gca; ax.FontSize = 30; ax.FontWeight = 'bold';

end
%%
%% -------------------------------------------------------------------------------------------------------
%% -------------------------------------------------------------------------------------------------------

function [pc_geneTypesUnique, pc_coeff2_cumulative] = makePcBiplotGeneTypes(cv1, pc_genes1, pc_coeff1, ensemblGenes, variableGenesTypes, colorCode)
% a biplot with contributions of different gene types (categories) to PC1
% and PC2.

% 1) find gene categories represented in the genes
% 2) find cumulative coefficient (pc1 and pc2) for all represented gene categories
% 3) plot

% 1)
[~, idx1, idx2] = intersect(ensemblGenes.names, pc_genes1);
pc_geneTypes = ensemblGenes.type(idx1);
pc_geneTypesUnique = unique(pc_geneTypes);

pc_coeff2 = abs(pc_coeff1(idx2,:)); %(get absolute values of pc coefficients and maintain same order as pc_geneTypes)

% 2)
pc_coeff2_cumulative = zeros(numel(pc_geneTypesUnique), 2);
for ii = 1:numel(pc_geneTypesUnique)
    idx1 = ismember(pc_geneTypes, pc_geneTypesUnique{ii});
    pc_coeff2_cumulative(ii,1:2) = sum(pc_coeff2(idx1,1:2),1);
end

% 3)
figure; hold on;
[~, colorIdx, ~] = intersect(variableGenesTypes, pc_geneTypesUnique);

h = biplot(pc_coeff2_cumulative(:,1:2), 'VarLabels', strrep(pc_geneTypesUnique, '_', ':'));
for ii = 1:numel(pc_geneTypesUnique)
    h(ii).Color = colorCode.geneTypes(colorIdx(ii),:);
    h(ii).LineWidth = 4;
    
    h(numel(pc_geneTypesUnique)+ii).Color = colorCode.geneTypes(colorIdx(ii),:);
    h(numel(pc_geneTypesUnique)*2+ii).String = ' ';
end

legend(strrep(pc_geneTypesUnique, '_', ':'), 'Box', 'off'); title(['CV=' num2str(cv1)]);
xlabel('Cumulative PC1 Coefficient');
ylabel('Cumulative PC2 Coefficient');
ax = gca; ax.FontSize = 14; ax.FontWeight = 'bold';

end
%% -------------------------------------------------------------------------------------------------------
%% -------------------------------------------------------------------------------------------------------



function lineageTableNew = retainExpressedGenes(ldm, lineages1, days1, lineageTable, ...
    plotExpressedGeneTypesPie, plotExpressedGeneTypesBar, ensemblGenes, variableGenesTypes, colorCodeG, ...
    fpkmThreshold, cellFractionThreshold, pseudoGenes)
%% find expresssed genes - genes with fpkm >1 in atleast half of the cells of every lineage

lineages = ldm.lineages;
days = ldm.days;
cellLineageDayIdx = ldm.cellIdx;

%%
if ~exist('fpkmThreshold', 'var')
    fpkmThreshold = 1;
end

if ~exist('cellFractionThreshold', 'var')
    cellFractionThreshold = 0.5;
end

%% select cells in desired lineages and days
[~, idx1,~] = intersect(lineages, lineages1);
cellLineageIdx = cat(1, cellLineageDayIdx{idx1,:}); % all cellIds corresponding to lineages1

[~, idx1,~] = intersect(days, days1);
cellDayIdx = cat(1, cellLineageDayIdx{:,idx1});

cellIdx1 = intersect(cellLineageIdx, cellDayIdx);

[~, idx1, ~] = intersect(lineageTable.Properties.VariableNames, cellIdx1);
lineageTable = lineageTable(:, [1 2 idx1']);

%% retain expresssed genes

expressedGenes = false(size(lineageTable,1),numel(lineages1));

for ii = 1:numel(lineages1)
    [~, idx1,~] = intersect(lineages, lineages1{ii});
    cellLineageIdx = intersect(cat(1, cellLineageDayIdx{idx1,:}), cellIdx1);
    
    [~, idx1, ~] = intersect(lineageTable.Properties.VariableNames, cellLineageIdx);
    fpkm1 = table2array(lineageTable(:, idx1));
    expressedGenes(:,ii) = sum(fpkm1>fpkmThreshold,2)>(cellFractionThreshold*size(fpkm1,2));
end
%
lineageTableNew = lineageTable(any(expressedGenes,2),:);
%%
% remove pseudoGenes
if exist('pseudoGenes', 'var')
    [~, idx1, ~] = intersect(table2cell(lineageTableNew(:,2)), pseudoGenes);
    idx2 = setxor([1:size(lineageTableNew,1)], idx1);
    lineageTableNew = lineageTableNew(idx2,:);
end
%%
%% a pie chart of expressed gene types in all lineages
if plotExpressedGeneTypesPie == 1
    figure; hold on;
    genes1 = table2cell(lineageTable(:,2));
    for ii = 1:numel(ldm.lineages)
        genes2 = genes1(expressedGenes(:,ii));
        [~, idx1, ~] = intersect(ensemblGenes.name, genes2);
        geneTypes1 = ensemblGenes.type(idx1);
        
        geneTypes1_unique = unique(geneTypes1);
        [geneTypes1_unique_ensembl, lIdx, ~] = intersect(variableGenesTypes, geneTypes1_unique);
        
        labelsCount=cellfun(@(x) sum(ismember(geneTypes1,x)),geneTypes1_unique_ensembl);
        [~,mIdx] = max(labelsCount); % position corresponding to most represented category
        explode = zeros(1, numel(geneTypes1_unique_ensembl)); % exploding it => separating biggest pie
        explode(mIdx) = 1;
        
        ax = subplot(3,4,ii); % => highlighting percentage corresponding to biggest pie
        p1 = pie(labelsCount, explode);
        p1(mIdx*2).FontSize = 14;
        p1(mIdx*2).FontWeight = 'bold';
        
        colormap(ax, colorCodeG.geneTypes(lIdx,:)); title([ldm.lineages{ii} ' nGenes =' num2str(numel(idx1))]);
        legend(strrep(geneTypes1_unique_ensembl, '_', ':'));
        ax.FontSize = 12; ax.FontWeight = 'bold';
        
    end
    
end

if plotExpressedGeneTypesBar == 1
    %% Two bar charts for all lineages
    %(a) fraction of pseudogenes, (b) fraction of protein coding genes
    
    %%
    % 1) replace gene names with gene types
    [~, idx1, idx2] = intersect(ensemblGenes.name, table2cell(lineageTable(:,2)));
    expressedGenesNew = expressedGenes(idx2,:);
    expressedGenesNew_geneTypes = ensemblGenes.type(idx1);
    %%
    % 2) save fraction of protein coding and pseudogenes
    % number of expressed genes corresponding to
    % [protein_coding, pseudogenes]
    
    proteinCodingGenes = zeros(1, numel(lineages1));
    pseudoGenes = proteinCodingGenes;
    for ii = 1:numel(lineages1)
        expressedGenesNew_geneTypes_lineage1 = expressedGenesNew_geneTypes(expressedGenesNew(:,ii));
        proteinCodingGenes(ii) = sum(contains(expressedGenesNew_geneTypes_lineage1, 'protein_coding'))./(numel(expressedGenesNew_geneTypes_lineage1));
        pseudoGenes(ii) = sum(contains(expressedGenesNew_geneTypes_lineage1, 'pseudogene'))./(numel(expressedGenesNew_geneTypes_lineage1));
    end
    %%
    toPlot = {proteinCodingGenes, pseudoGenes};
    if numel(lineages1) >= 8
        columnOrder = [1 6 9 2:5 7 8];
    else
        columnOrder = [1:numel(lineages1)];
    end
    
    for ii = 1:numel(toPlot)
        figure;
        values = toPlot{ii}(columnOrder);
        idx1 = values>0; % remove lineages not expressed
        values = values(idx1);
        b = bar(values);
        b.FaceColor = [0.5 0.5 0.5];
        ax = gca;
        ax.XTickLabel = ldm.lineages(columnOrder(idx1));
        ax.FontSize = 20;
        ax.FontWeight = 'bold';
    end
end


end
%% --------------------------------------------------------------------------------------------------------
%% --------------------------------------------------------------------------------------------------------
%% --------------------------------------------------------------------------------------------------------

function plotExpressedGenesTypesSingleCells(ldm, lineage1, day1, lineageTable, ...
    ensemblGenes, variableGenesTypes, colorCodeG,  pieChart, barChart, fpkmThreshold)

if ~ exist('fpkmThreshold', 'var')
    fpkmThreshold = 1;
end


data1_genes = table2cell(lineageTable(:,2));
data1_fpkm = table2array(lineageTable(:,3:end));


[~, ~, idx_lineage] = intersect(lineage1, ldm.lineages);
[~, ~, idx_day] = intersect(day1, ldm.days);
names_cells = ldm.cellIdx{idx_lineage, idx_day};
[names_cells, idx_cells, ~] = intersect(lineageTable.Properties.VariableNames, names_cells);

% if pie chart
if (pieChart)
    
    figure_nColumns = 4;
    figure_nRows = round(numel(names_cells)/figure_nColumns);
    if figure_nRows<(numel(names_cells)/figure_nColumns)
        figure_nRows = figure_nRows+1;
    end
    
    
    figure; hold on; % pie chart
    
    for ii = 1:numel(names_cells)
        expressedGenes_names = data1_genes(data1_fpkm(:,idx_cells(ii)-2) >= fpkmThreshold);
        [~, ~, idx2] = intersect(expressedGenes_names, ensemblGenes.name);
        expressedGenes_type = ensemblGenes.type(idx2);
        expressedGenes_type_unique = unique(expressedGenes_type);
        
        [expressedGenes_type_ensembl, colors_idx1, ~] = intersect(variableGenesTypes, expressedGenes_type_unique);
        labelsCount=cellfun(@(x) sum(ismember(expressedGenes_type,x)),expressedGenes_type_ensembl);
        
        [~,mIdx] = max(labelsCount); % position corresponding to most represented category
        explode = zeros(1, numel(expressedGenes_type_ensembl)); % exploding it => separating biggest pie
        explode(mIdx) = 1;
        
        ax = subplot(figure_nRows,figure_nColumns,ii); % => highlighting percentage corresponding to biggest pie
        p1 = pie(labelsCount, explode);
        p1(mIdx*2).FontSize = 14;
        p1(mIdx*2).FontWeight = 'bold';
        
        colormap(ax, colorCodeG.geneTypes(colors_idx1,:));
        title([ldm.lineages{idx_lineage} names_cells{ii} ' nGenes =' num2str(numel(expressedGenes_type))]);
        %legend(strrep(expressedGenes_type_ensembl, '_', ':'));
        ax.FontSize = 12; ax.FontWeight = 'bold';
    end
end

if (barChart)
    % else
    % bar chart 1*2; [fraction of gene types, number of expressedGenes]
    
    genes_types_req = {'lncRNA', 'processed_pseudogene', 'protein_coding'};
    [genes_types_req, colors_idx1, ~] = intersect(variableGenesTypes, genes_types_req);
    genes_types_req = cat(1, genes_types_req, {'rest'});
    %%
    genes_types_fraction = zeros(4, numel(names_cells));
    genes_count = zeros(1, numel(names_cells));
    
    for ii = 1:numel(names_cells)
        expressedGenes_names = data1_genes(data1_fpkm(:,idx_cells(ii)-2) >= fpkmThreshold);
        [~, ~, idx2] = intersect(expressedGenes_names, ensemblGenes.name);
        expressedGenes_type = ensemblGenes.type(idx2);
        
        expressedGenes_type(~ismember(expressedGenes_type, genes_types_req)) = {'rest'};
        
        labelsCount=cellfun(@(x) sum(ismember(expressedGenes_type,x)),genes_types_req);
        labelsCount = labelsCount./(sum(labelsCount));
        genes_types_fraction(:,ii) = labelsCount;
        genes_count(ii) = numel(expressedGenes_type);
    end
    
    %%
    figure; hold on;
    subplot(1,2,1);
    b = bar((1:numel(names_cells)), genes_types_fraction, 'stacked');
    for ii = 1:3
        b(ii).FaceColor = colorCodeG.geneTypes(colors_idx1(ii),:);
        b(ii).BarWidth = 0.8;
    end
    title([lineage1 ' ' day1]); legend(strrep(genes_types_req, '_', ':'));
    ax = gca; ax.FontSize = 20; ax.FontWeight = 'bold';
    
    subplot(1,2,2);
    b = bar((1:numel(names_cells)), genes_count);
    b.FaceColor = [0.5 0.5 0.5];
    b.BarWidth = 0.8;
    title([lineage1 ' ' day1]);
    ax = gca; ax.FontSize = 20; ax.FontWeight = 'bold';
    
end

end
%%





































