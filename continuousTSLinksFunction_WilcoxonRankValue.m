% this is where I find the tissue specificity for each link with the wilcoxon
% Rank test. 

function[] = continuousTSLinksFunction_WilcoxonRankValue(tissue, s, n)

load(['~/data/general/tissueExpGenes/' tissue 'ExpGenes0.8.mat']);
load('~/data/general/linkExprInfo/wholeGeneExpr.mat')
load('~/data/general/linkExprInfo/dataSetProbeInf.mat')

finalMat = zeros(n, 53);
for i = 1:53
    i
    dsName = dataSetProbeInf(i).name;
    dsMat = wholeGeneExpr(:, dataSetProbeInf(i).sampleInd(1): ...
                          dataSetProbeInf(i).sampleInd(2));
    dsMat = dsMat(expGenesInd, :);
    %sib = corr(dsMat(1:100, :)');
    sib = corr(dsMat');
    upOnes = sib(logical(triu(ones(size(sib)), 1)));
    clear sib
    arr = upOnes(s:s+n-1);
    clear upOnes
    finalMat(:, i) = corrQuantileReturn(dataSetProbeInf(i).name, ...
                                        arr);
end
% find the tissueInd

tissueCell = {dataSetProbeInf(:).tissue};
dsIndCell = strfind(tissueCell, tissue);
allDSInd = find(~cellfun(@isempty, dsIndCell));

tic
allTissues = [1:53];
pVals = zeros(n, 1);
tissueInd = allDSInd;
otherInd = allTissues(~ismember(allTissues, tissueInd));
tc = length(tissueInd);
oc = length(otherInd);

for p = 1: n
    small = finalMat(p,tissueInd);
    big = finalMat(p, otherInd);
    pVals(p) = ranksum(small, big, 'tail', 'right');
end

fileName = sprintf('~/networks/continuousTSLinksResults/OneSided%s_%d_%d_pVals.mat', ...
                   tissue, s, n)
save(fileName, 'pVals');

toc
end
