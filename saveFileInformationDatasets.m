
%% 
% save zscores of group gene counts(mean count for a cell type) and single
% cell gene counts in a separate structure for performing correlation
% analyses across different samples
%% ----------------------------------------------------------------------------
% 1) Human Embryo D14 Xiang et al
%% ----------------------------------------------------------------------------%%
masterFolder = '/Users/sapnachhabra/Desktop/CellTrackercd/Experiments/rnaSeqExperiments/rnaSeq_datasets_fpkmTable';
ldmFile =[masterFolder filesep  'humanEmbryoAmnionDataD14/forReport' filesep 'forReport_singleCells.mat']; % cells with ldm information
genesFile = [masterFolder filesep 'humanEmbryoAmnionDataD14/forReport' filesep 'forReport_genes.mat'];

load(ldmFile, 'geneCellTable', 'ldm');
load(genesFile);

%% ----------------------------------------------------------------------------
% edit samples 6,7 in the structure datasets
%% ----------------------------------------------------------------------------

%% a) dataset6: readcount for single cells

toKeep = find(any(table2array(geneCellTable(:,3:end))>0, 2)); % exclude genes not expressed in any cell
geneCellTable_new = geneCellTable(toKeep,:);
%%

genes6 = table2cell(geneCellTable_new(:,2));
rawReads6 = table2array(geneCellTable_new(:,3:end));
normReads6 = zscore(log2(rawReads6+1), 0, 'all');

%%
dId = 6;
datasets(dId).id = 'humanEmbryoD14_singleCell_allLineages';
datasets(dId).genes = genes6;
datasets(dId).samples = geneCellTable_new.Properties.VariableNames(3:end);
datasets(dId).rawReads = rawReads6;
datasets(dId).normReads = normReads6;


%% b) dataset7: readcount for cell group

% compute mean readcount for each group (specificDay_specificCellType)
counter = 1;
for ii = 1:numel(ldm.lineages)
    
    dayIdx = find(~cellfun(@isempty,ldm.cellIdx(ii,:)));
    
    for jj = dayIdx
        
        cellIdx = ldm.cellIdx{ii,jj};
        [~, idx1, ~] = intersect(geneCellTable_new.Properties.VariableNames, cellIdx);
        rawReads1 = mean(table2array(geneCellTable_new(:,idx1)),2);
        rawReads7(:,counter) = rawReads1;
        sampleNames{counter} = strcat(ldm.days{jj},ldm.lineages{ii});
        sampleCounts(:,counter) = numel(cellIdx);
        counter = counter+1;
    end
end

normReads7 = zscore(log2(rawReads7+1), 0, 'all');
%%
dId = 7;
datasets(dId).id = 'humanEmbryoD14_singleCells_allLineages';
datasets(dId).genes = genes6;
datasets(dId).samples = sampleNames;
datasets(dId).rawReads = rawReads7;
datasets(dId).normReads = normReads7;

%%
dataset1 = datasets(7);
dataset2 = dataset1;
sampleIds1 = 1:35;
sampleIds2 = 1:35;
rawReads = 0;

geneList1 = geneLists.humanEmbryoLineage_revised;
geneList2 = plotCorrelationCoefficients(dataset1, dataset2, sampleIds1, sampleIds2,...
    geneList1, rawReads);






















