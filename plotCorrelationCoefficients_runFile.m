

%%
datasetsFile = '/Users/sapnachhabra/Desktop/CellTrackercd/Experiments/rnaSeqExperiments/rnaSeq_datasets_fpkmTable/datasets.mat';
load(datasetsFile);
%%
caxisLimits = [0 1];
datasetIdx1 = 12; datasetIdx2 = 10; 
samplesIdx1 = [11 7 8 6 10 9 4 2 1 5 3]; samplesIdx2 =[1:43];
rawReads = 0;
%%
counter = 1;
for ii = [2 4 6]
    geneList1{counter} =  variableGenesPCA_proteinCoding(ii).allLineages.genes;
    cv(counter) = variableGenesPCA_proteinCoding(ii).allLineages.cvThreshold;
    counter = counter+1;
end
%%
for ii = 3
%for ii = 1:numel(geneList1)
    plotCorrelationCoefficients(datasets(datasetIdx1), datasets(datasetIdx2),...
        samplesIdx1, samplesIdx2, geneList1{ii}, rawReads, caxisLimits);
    title(['cv' num2str(cv(ii)) 'nGenes' num2str(numel(geneList1{ii}))]);
end

%%
caxisLimits = [0 0.6];
rawReads = 0;
geneList1 = table2cell(readtable('plosbioPaperData.xlsx', 'Sheet', 2));

datasetIdx1 = 10; 
samplesIdx1 = [1:43]; %[1:4 26:29 38:41]%26:31 38:43 5:8 21:25 32:37 10:15 9];

datasetIdx2 = 4;
samplesIdx2 = [1:5]%[1:4 16:20 26:31 38:43 5:8 21:25 32:37 10:15 9];%[3:15 23:29] ;

plotCorrelationCoefficients(datasets(datasetIdx1), datasets(datasetIdx2),...
        samplesIdx1, samplesIdx2, geneList1, rawReads, caxisLimits);

%% ---------------------------------------------------------------------------------------------------------
%% ------------------------- with ebseqFiles ---------------------------------------------------------------
%% ---------------------------------------------------------------------------------------------------------

lineageGenes_ebseq = lineageGenes_ebseq_epiCtbPrED79;

caxisLimits = {[0 1]; [0 0.6]}; rawReads = 0;
dataset1_indices = [8, 4];
samples1_indices =  {[1:35]; [1:4]};  %{[1:4 26:29 38:41]; [1:4]};

datasetIdx2 = 8;
samplesIdx2 = [1:35];% human;  %[1:4 26:29 38:41];monkey



stats.probDe = cell2mat({lineageGenes_ebseq.probDe});
stats.fpkm = cell2mat({lineageGenes_ebseq.fpkm});
stats.foldChange = cell2mat({lineageGenes_ebseq.foldChange});

stats.probDe_unique = unique(stats.probDe);
stats.fpkm_unique = unique(stats.fpkm);
stats.foldChange_unique = unique(stats.foldChange);
stats
%%
geneList_De_hM = cell(2,1);
counter = 1;
%geneCellTableNew = geneCellTable;

for ii = stats.probDe_unique(2)
    for jj = stats.fpkm_unique(2)
        for kk = stats.foldChange_unique([1:5])
            idx = stats.probDe == ii & stats.fpkm == jj & stats.foldChange == kk;
            geneList.ensemblId = lineageGenes_ebseq(idx).ensemblId;
            
            [~, idx1, ~] = intersect(table2cell(geneCellTableNew(:,1)), geneList.ensemblId);
            geneList1 = table2cell(geneCellTableNew(idx1,2));
            
            for ll = 1:numel(dataset1_indices)
                datasetIdx1 = dataset1_indices(ll);
                samplesIdx1 = samples1_indices{ll};
                geneList2 = plotCorrelationCoefficients(datasets(datasetIdx1), datasets(datasetIdx2),...
                    samplesIdx1, samplesIdx2, geneList1, rawReads, caxisLimits{ll});
                title(['probDe' num2str(ii) 'fpkm' num2str(jj) 'fc' num2str(kk) 'nGenes=' num2str(numel(geneList2))]);
            end
            
            
            geneList.ensemblId_separate = getfield(lineageGenes_ebseq(idx), 'ensemblId_separate');
            for mm = 1:numel(geneList.ensemblId_separate)
                [~, idx1, ~] = intersect(table2cell(geneCellTableNew(:,1)), geneList.ensemblId_separate{mm});
                 geneList1 = table2cell(geneCellTableNew(idx1,2));
                 [~, idx1,~] = intersect(geneList1, geneList2); % only keep genes expressed in both human & monkey
                 geneList_De_hM{counter}{mm} = sort(geneList1(idx1));
            end
            counter = counter+1;
        end
    end
end
%%
%%
%% save genelists as a text file
% 1) individual gene lists for each lineage
lineages = {'ctb', 'epi', 'intermediate'};
for ii = 1:3
    writetable(cell2table(geneList_De_hM{1}{ii}), [lineages{ii} '.txt'], 'WriteVariableNames', 0);
end
%%
% 2) all lineage combined
for ii = 1:5
    fc = stats.foldChange_unique(ii);
    geneList1 = unique(cat(1, geneList_De_hM{ii}{:}));
    writetable(cell2table(geneList1), ['fc' num2str(fc) '.txt'], 'WriteVariableNames', 0);
end


%% ---------------------------------------------------------------------------------------------------------
%% ------------------------- -------------------------------------------------------------------------------
%% ---------------------------------------------------------------------------------------------------------


%% ---------------------------------------------------------------------------------------------------------
%% ------------------------- -------------------------------------------------------------------------------
%% ---------------------------------------------------------------------------------------------------------
%%
% how similar are single cells to mean expression?
stats.probDe = cell2mat({lineageGenes_ebseq.probDe});
stats.fpkm = cell2mat({lineageGenes_ebseq.fpkm});
stats.foldChange = cell2mat({lineageGenes_ebseq.foldChange});

idx = find(stats.probDe == 0.95 & stats.fpkm == 50 & stats.foldChange == 5);
geneList.ensemblId = getfield(lineageGenes_ebseq(idx), 'ensemblId');
[~, idx1, ~] = intersect(table2cell(geneCellTableNew(:,1)), geneList.ensemblId);
geneList1 = table2cell(geneCellTableNew(idx1,2));

%%
% single cells on a specific day and lineage with (a)others in the same group
% (datasetIndices2{1} samplesIndices2{1}); (b) mean expression in different
% lineages, days (datasetIndices2{2} samplesIndices2{2}).

geneList1 = geneLists.humanEmbryoLineage_revised;
datasetIdx1 = 6;
datasetIndices2 = {6,7};
samplesIndices2{2} = [1:35];
caxisLimit = [0 1];
rawReads = 0;

datasetIdx1_lineageIdx = [6];

for ii = datasetIdx1_lineageIdx
    dayIdx = find(~cellfun(@isempty, ldm.cellIdx(ii,:)));
    
    for jj = dayIdx
       
        cellIdx = ldm.cellIdx{ii,jj};
        [~, ~, idx2] = intersect(cellIdx, datasets(datasetIdx1).samples, 'stable');
        samplesIndices2{1} = idx2;
        samplesIdx1 = idx2;
        
        for kk = 1:2
            datasetIdx2 = datasetIndices2{kk};
            samplesIdx2 = samplesIndices2{kk};
            plotCorrelationCoefficients(datasets(datasetIdx1), datasets(datasetIdx2),...
                samplesIdx1, samplesIdx2, geneList1, rawReads, caxisLimit);
            title([ldm.lineages{ii} ldm.days{jj} 'nGenes' num2str(numel(geneList1))]);
        end
    end
end
%%
% single cells on a specific day and lineage with mean expression in different
% lineages on the same days
geneList1 = variableGenesPCA(1).allLineages.genes
datasetIdx1 = 7; % single cells
datasetIdx2 = 10; % mean expression
caxisLimit = [0 1]; rawReads = 1;

datasetIdx1_lineageIdx = 1;

for ii = datasetIdx1_lineageIdx
    dayIdx = find(~cellfun(@isempty, ldm.cellIdx(ii,:)));
    
    for jj = dayIdx
        cellIdx = ldm.cellIdx{ii,jj};
        [~, idx1,~] = intersect(datasets(7).samples, cellIdx);
        samplesIdx1 = idx1;
        samplesIdx2 = find(contains(datasets(datasetIdx2).samples, ldm.days{jj})...
            & ~contains(datasets(datasetIdx2).samples, ldm.lineages{ii}));
        plotCorrelationCoefficients(datasets(datasetIdx1), datasets(datasetIdx2),...
            samplesIdx1, samplesIdx2, geneList1, rawReads, caxisLimit);
        title([ldm.lineages{ii} ldm.days{jj}]);
    end
end

%%
saveInPath = ['/Users/sapnachhabra/Desktop/CellTrackercd/Experiments/rnaSeqExperiments/rnaSeq_datasets_fpkmTable/'...
    'humanEmbryoD1619/correlationPlots/variableGenes/'];
saveAllOpenFigures(saveInPath)






