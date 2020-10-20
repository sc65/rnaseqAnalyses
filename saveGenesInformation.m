

load('/Users/sapnachhabra/Desktop/CellTrackercd/Experiments/rnaSeqExperiments/rnaSeq_datasets_fpkmTable/humanEmbryoAmnionDataD14/forReport/cellColors.mat');
%% ------------------------------------------------------------------------------------------
%% ----------------------------------add genelists-------------------------------------------
%% ------------------------------------------------------------------------------------------
geneLists.humanEmbryoLineage_published = {'NANOG', 'SOX2', 'KLF17', 'TDGF1', ...
    'PDGFRA', 'GATA6', 'GATA4', 'SOX17', 'GATA2', 'GATA3', 'TEAD3', 'KRT18'};
geneLists.humanEmbryoLineage_revised = {'NANOG', 'SOX2', 'KLF17', 'TDGF1', ...
    'PDGFRA', 'GATA6', 'GATA4', 'SOX17', 'GATA2', 'GATA3', 'TEAD3', 'KRT7'};

geneLists.pluripotencyCore = {'NANOG', 'POU5F1', 'SOX2'};
geneLists.humanNaiveEpi = {'SOX2', 'POU5F1', 'NANOG', 'KLF2', 'KLF4', 'KLF17',...
    'TFCP2L1', 'DPPA3', 'ESRRB', 'TBX3', 'SALL4', 'GBX2'};
geneLists.primitiveStreak = {'MIXL1', 'TBXT'};

geneLists.trophectoderm = {'CDX2', 'ELF5', 'GATA2', 'GATA3', 'KRT7', 'TEAD4', 'TFAP2A', 'TFAP2C', 'TP63'};
geneLists.trophectoderm_stb = {'PSG1', 'CSH1', 'HSD3B1', 'CYP19A1', 'SDC1', 'INHA'};
geneLists.trophectoderm_stb_hcg = {'CGA', 'CGB1' ,'CGB2', 'CGB3', 'CGB5', 'CGB7', 'CGB8'};
geneLists.trophectoder_evt = {'HLA-G'  'ITGA5'};

geneLists.ExM = {'GATA4', 'GATA6', 'COL6A1', 'CDH2', 'VIM', 'SNAI1'};

geneLists.amnionMonkey = {'TFAP2C', 'MSX2', 'BMP4', 'ID2', 'AXIN2', 'DNMT3A'};
%%
%% ------------------------------------------------------------------------------------------
%% ----------------------------------add ensemblGeneInfo-------------------------------------
%% ------------------------------------------------------------------------------------------
geneFile = readtable('ensemblGenesGeneTypeInfo.txt');
ensemblGenes.ID = table2cell(geneFile(:,1));
ensemblGenes.name = table2cell(geneFile(:,2));
ensemblGenes.type = table2cell(geneFile(:,3));
%%

geneTypes = unique(ensemblGenes.type);
colorCodeG.geneTypes = moreColors;

%%
save('forReport_genes.mat', 'ensemblGenes', 'geneLists', 'geneTypes', 'colorCodeG');
%%