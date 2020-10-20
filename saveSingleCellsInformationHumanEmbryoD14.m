
%% make all file names, plotPSAcellCountPie as function arguments
%
plotPSAcellCountPie = 1;
plotAllCellCountHeatMap = 1;
%%
fpkmFile = readtable('GSE136447_555-samples-fpkm.txt');
allCellIdsFile = readtable('41586_2019_1875_MOESM10_ESM.xlsx', 'Sheet',1); % all samples
psaCellIdsFile = readtable('41586_2019_1875_MOESM10_ESM.xlsx', 'Sheet',3); % epi, amnion, ps

load('/Users/sapnachhabra/Desktop/CellTrackercd/Experiments/rnaSeqExperiments/rnaSeq_datasets_fpkmTable/humanEmbryoAmnionDataD14/forReport/cellColors.mat');
%%
%% ------------------------------------------------------------------------------------------
%% -------------------------------relabel cells (stratify PSA)-------------------------------
%% ------------------------------------------------------------------------------------------

% save indices for psa-epi
psaIdx = ismember(table2cell(allCellIdsFile(:,4)), {'PSA-EPI'});
%%
% replace d14 epi & psa-epi with d14 epi/intermediate/amnion
cellIdx1 = table2cell(allCellIdsFile(:,1));
cellIdx2 = table2cell(psaCellIdsFile(:,1));
%
[~, idx1, idx2] = intersect(cellIdx1, cellIdx2);% find intersecting cell ids in file 1 and 2
allCellIdsFile(idx1,4) = psaCellIdsFile(idx2, 3);
%%
% add ps-epi
idx2 = psaIdx & ismember(table2cell(allCellIdsFile(:,4)), {'EPI'});
allCellIdsFile(idx2,4) = {'psEPI'};

%% -----------------------------------------------------------------------------------------
if (plotPSAcellCountPie)
    % (optional)
    % a) how many cells correspond to day 14 AME, Intermediate and EPI?
    % order: [AME, EPI, Intermediate]
    idx3 = ismember(table2cell(psaCellIdsFile(:,2)), {'D14'}) & ...
        ismember(table2cell(psaCellIdsFile(:,3)), {'AME', 'EPI',  'Intermediate'});
    ameInterEpiD14 = table2cell(psaCellIdsFile(idx3,3));
    ameInterEpiD14_count=cellfun(@(x) sum(ismember(ameInterEpiD14,x)),unique(ameInterEpiD14));
    %%
    % b) what is psa composed of?
    % order: [AME, EPI, Intermediate]
    psa_breakDown = table2cell(allCellIdsFile(psaIdx,4));
    psa_breakDown = strrep(psa_breakDown, 'psEPI', 'EPI');
    psa_breakDown_count=cellfun(@(x) sum(ismember(psa_breakDown,x)),unique(psa_breakDown));
    
    %%
    figure; p = pie(categorical(psa_breakDown));
    hold on;
    counter = 1;
    for ii = [2 4 6]
        p(ii).String = [p(ii-1).DisplayName ' (' ...
            num2str(psa_breakDown_count(counter)) '/' num2str(ameInterEpiD14_count(counter)) ')'];
        p(ii).FontWeight = 'bold';
        p(ii).FontSize = 16;
        counter = counter+1;
    end
end
%%
%% ------------------------------------------------------------------------------------------
%% -------------------------------make lineageDayMatrix (ldm)--------------------------------
%% ------------------------------------------------------------------------------------------
%% structure ldm: lineageDayMatrix (lineage * day), fields: lineages, days, cellCounts,
% cellIds.
% structure colorCoden, fields lineages, days

%%
% initialize
ldm.lineages = table2cell(unique(allCellIdsFile(:,4)));
ldm.days =strcat('D', strsplit(int2str([6:10 12 14]), ' '));
ldm.cellCounts = zeros(numel(ldm.lineages), numel(ldm.days));
ldm.cellIdx = cell(numel(ldm.lineages), numel(ldm.days));


for ii = 1:numel(ldm.lineages)
    for jj = 1:numel(ldm.days)
        idx1 = find(ismember(table2cell(allCellIdsFile(:,4)), ldm.lineages{ii})& ...
            ismember(table2cell(allCellIdsFile(:,2)), ldm.days{jj}));
        ldm.cellIdx{ii,jj} = table2cell(allCellIdsFile(idx1,1));
        ldm.cellCounts(ii,jj) = numel(idx1);
    end
end

%% ------------------------------------------------------------------------------------------
if (plotAllCellCountHeatMap)
    figure; h = heatmap(ldm.cellCounts);
    h.XData = ldm.days;
    h.YData = ldm.lineages;
    ax = gca;
    ax.FontSize = 20;
end
%%
%% ------------------------------------------------------------------------------------------
%% -------------------------------make geneCellTable-----------------------------------------
%% ------------------------------------------------------------------------------------------
% retain unique genes in geneCellTable. For genes with multiple
% transcripts, retain the cumulative read count of all transcripts
% and geneId of most expressed transcript.


[~, ~, cellIdx] = intersect(cat(1,ldm.cellIdx{:}), fpkmFile.Properties.VariableNames);
geneCellTable1 = fpkmFile(:,[2 3 cellIdx']); % retain ensemblGeneIds, geneNames, geneCounts

%%
% time alert: takes ~4 minutes to run
geneList1 = unique(table2cell(geneCellTable1(:,2))); % unique genes
geneCounts = zeros(numel(geneList1), numel(cellIdx));
geneIds = cell(numel(geneList1), 2); % [ensemblGeneId, geneName]

tic;
for ii = 1:numel(geneList1)
    geneIdx1 = find(ismember(table2cell(geneCellTable1(:,2)), geneList1{ii}));
    geneCounts(ii,:) = sum(table2array(geneCellTable1(geneIdx1, 3:end)),1); % sum of all transcript counts
    
    if numel(geneIdx1) > 1
        [~, idx1] = max(sum(table2array(geneCellTable1(geneIdx1, 3:end)), 2));
        geneIdx1 = geneIdx1(idx1); % index of most expressed transcript
    end
    
    geneIds(ii,1:2) = table2cell(geneCellTable1(geneIdx1,1:2)); % ensembl id of most expressed transcript
end
toc;
%%
geneCellTable2 =  [cell2table(geneIds) array2table(geneCounts)];
for ii = 1:size(geneCellTable2,2)
    geneCellTable2.Properties.VariableNames{ii} = geneCellTable1.Properties.VariableNames{ii};
end
%%
geneCellTable = geneCellTable2;
%%
%% -----------------------------------update ldm---------------------------------------------

% add additional fields in ldm corresponding to lineage, day and idx of
% cells in the same order as columns of geneCellTable - useful for later
% analyses

ldm.lineageLabels = cell(1, numel(geneCellTable2.Properties.VariableNames(3:end)));
ldm.dayLabels = cell(1, numel(geneCellTable2.Properties.VariableNames(3:end)));
ldm.cellLabels = geneCellTable2.Properties.VariableNames(3:end);

for ii = 1:numel(ldm.lineages)
    lineage1_cells = cat(1,ldm.cellIdx{ii,:});
    [~, ~, idx1] = intersect(lineage1_cells, geneCellTable2.Properties.VariableNames);
    ldm.lineageLabels(idx1'-2) = ldm.lineages(ii);
end
    
   
for ii = 1:numel(ldm.days)
    day1_cells = cat(1,ldm.cellIdx{:,ii});
    [~, ~, idx1] = intersect(day1_cells', geneCellTable2.Properties.VariableNames);
    ldm.dayLabels(idx1-2) = ldm.days(ii);
end
%% ------------------------------------------------------------------------------------------
colorCode = colors_humanEmbryoCells;
%% ------------------------------------------------------------------------------------------
save('forReport_singleCells.mat', 'colorCode', 'ldm', 'geneCellTable');




