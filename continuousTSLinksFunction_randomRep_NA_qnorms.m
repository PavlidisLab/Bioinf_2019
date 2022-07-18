% this is the main function for running the code in parallel and
% saving the data. 
% where tissue is the tissue and s is the starting point for the
% correlation matrix and n is the number of edges to be checked
% (cause I am running in parallel, so I choose how long each run takes)
%function[] = continuousTSLinksFunction(tissue, s, n)
% s is the starting point (cause I get in limited number of edges
% each time) - n is the count of edges , I have ~57e+6 edges, for
% brain, for example. 1e+6 takes about half an hour to run. I give
% it 10mil at a time, I guess. 10e+6 takes 5 hours, and probably
% 10G of RAM. 

function[] = continuousTSLinksFunction_randomRep_NA_qnorms(dataSetList, ...
                                                  s, n, r, rs)

load('~/data/general/tissueExpGenes/allDSExprInfo.mat')
load('~/data/general/linkExprInfo/wholeGeneExpr.mat')
load('~/data/general/linkExprInfo/dataSetProbeInf.mat')
tic

totalDSCount = length(dataSetProbeInf)
% getting the expressed genes
exprMat = zeros(18494, length(dataSetList));
[a, dsExpr] = ismember(dataSetList, {allDSExprInfo(:).GSEID});
for i = 1:length(dataSetList)
    exprMat(:, i) = allDSExprInfo(dsExpr(i)).exprGenes;
end

dsThr = ceil(.8 * length(dataSetList))
sib = sum(exprMat');
expGenesInd = sib >= dsThr;

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
    finalMat(:, i) = corrQNormReturn(dataSetProbeInf(i).name, ...
                                     arr);

end
% find the tissueInd

sib = finalMat < 0;
finalMat(sib) = 0;
clear sib

% tissueCell = {dataSetProbeInf(:).tissue};
% dsIndCell = strfind(tissueCell, tissue);
% allDSInd = find(~cellfun(@isempty, dsIndCell));

%getting the IDs for the "tisssue" and "others"
AllIDs = {dataSetProbeInf(:).name};
[a, allDSInd] = ismember(dataSetList, AllIDs);
totalDSCount = length(dataSetProbeInf);

allTissues = [1:53];
space = zeros(n, 1);
tissueInd = allDSInd;
otherInd = allTissues(~ismember(allTissues, tissueInd));
sum(tissueInd)
sum(otherInd)

for i = 1:length(tissueInd)
    for j = 1:length(otherInd)
        sib = finalMat(:, tissueInd(i)) - finalMat(:, otherInd(j));
        sib(sib < 0) = 0;
        space = space + sib;
    end
end
toc

fileName = sprintf('~/networks/continuousTSLinksResults/%d_random%d_%d_%d_NA_qnorms.mat', ...
                   rs, s, n, r)

save(fileName, 'space');
end