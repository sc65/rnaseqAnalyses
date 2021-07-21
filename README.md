# rnaseqAnalyses
Xiang et al 2019, Data re-analyses 

1.	Cell labelling
Cell labels corresponding to each cell type were determined from supplementary table 8 (S8.1) in Xiang et al. Cells with the primitive streak anlage (PSA) label were further stratified into AME, Intermediate or Primitive streak epi cells based on information in S8.3. 
Information saved in the structure ldm (lineageDayMatrix)

2.	Read counts and gene labels
Raw read counts for each cell were obtained from supplementary file deposited in NCBI GEO under accession number GSE136447. All the analyses were performed at the level of genes. For genes with multiple transcripts, cumulative expression of all transcripts was considered as the gene read count and the ensembl gene id of most expressed transcript was considered as its gene id. This gave a total of 58683 genes. Gene ids and names of these genes, along with their read counts for all all 555 cells was saved in the table geneCellTable. 
ldm and geneCellTable saved in the file forReport_singleCells.mat

3.	Gene type specification
Human genes with the following attributes – gene stable ID, gene name and gene type, were downloaded from ensembl biomart. This gave a total of 67130 genes. 
geneLists saved in the file forReport_genes.mat

4.	Comparing across different datasets 
To compare gene expression across datasets (one dataset = 1 independent rnaseq study or 1 version (average expression across cells of a given cell type or single cell expression data for single cell rna seq) of an independent rnaseq study), all datasets were saved in the structure datasets with fields ‘id’, ‘genes’, ‘samples’, ‘rawReads’ and ‘normReads’. Genes not expressed (fpkm/rpkm=0) in any cell within a given dataset were excluded. 

‘id’: dataset identifier (informs about species, embryonicTime and single cell/average gene expression)
‘genes’: gene names
‘samples’: cell identifier
‘rawReads’: expression in fpkm/rpkm
‘normReads’: zscores of log2(rawReads)

Gene names correspond to human genes.  

The 11 datasets in the structure correspond to the following studies - Petropoulos et al 2016, Okae et al 2018, Haider et al 2018, Chhabra et al 2019, Nakamura et al 2016, Xiang et al 2019 (single cells), Xiang et al 2019 (mean), Ma et al 2019 (single cells), Ma et al 2019 (mean), Tyser et al 2020 biorxiv (single cells), Tyser et al 2020 biorxiv (mean). 


Codes

1.	saveSingleCellsInformationHumanEmbryoD14
-saves ldm (lineage day matrix), geneCellTable 
-runs (a) plotPSAcellCountsPie. Figure 1A
-runs (b) plotAllcellCountsHeatmap. Figure 1B
2.	saveGenesInformation
- saves geneLists: certain lineage specific genes, ensemblGenes: genelist from biomart, geneTypes: different gene categories of ensembl genes, colorCodeG: one color/geneType
3.	plotGeneCountHeatmapBoxScatter
4.	plotAndSavePCAresults
5.	saveFileInformationDatasets
- save raw and normalized read counts for comparison across datasets.
6.	plotCorrelationCoefficients
- computes and plots pearson correlation coefficients for samples saved in the structure datasets. 



New annotations:
ICM – correlation coefficients of expression of lineage specific genes in human cells.
Intermediate – PCA based on variable protein coding genes (PrE derived) -> Monkey dataset to infer VE/EXMC. Cells which do not maximally correlate with human PrE based on expression of lineage specific genes labelled as ambiguous. 
AME – PCA and correlation coefficients of CTB/STB/EVT specific genes. Extensive details in the text. 
