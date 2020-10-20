function geneList1 = plotCorrelationCoefficients(dataset1, dataset2, sampleIds1, sampleIds2, geneList, rawReads, caxisLimits)
%% ----- function to plot similarity matrix
% -- inputs: datasets, sampleIds, geneList
% -- output: similarity matrix
%%
% extract read count values
[genes1, idx1,~] = intersect(dataset1.genes, geneList);
if rawReads == 0
    reads1 = dataset1.normReads(idx1, sampleIds1);
else
    reads1 = log2(dataset1.rawReads(idx1, sampleIds1)+1);
end

[geneList1, idx2,idx1] = intersect(dataset2.genes, genes1);
if rawReads == 0
    reads2 = dataset2.normReads(idx2, sampleIds2);
else
    reads2 = log2(dataset2.rawReads(idx2, sampleIds2)+1);
end

reads1 = reads1(idx1,:);
%%
% compute correlation coefficients
ccMatrix = zeros(numel(sampleIds1), numel(sampleIds2));
for ii = 1:numel(sampleIds1)
    for jj = 1:numel(sampleIds2)
        cc = corrcoef(reads1(:,ii), reads2(:,jj));
        ccMatrix(ii,jj) = round(cc(1,2),2);
    end
end
%%
if ~exist('caxisLimits', 'var')
    caxisLimits = [min(ccMatrix(:)), max(ccMatrix(:))];
end
% plots
figure;
h = heatmap(ccMatrix); colormap('jet');
caxis([caxisLimits]);
h.YDisplayLabels = strrep(dataset1.samples(sampleIds1), '_', ':');
h.CellLabelColor = 'none';
h.XDisplayLabels = strrep(dataset2.samples(sampleIds2), '_', ':');
ax = gca; ax.FontSize = 20;

end
