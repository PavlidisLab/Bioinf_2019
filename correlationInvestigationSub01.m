% this file is to get the absolute mean correlation among the
% genes. I am just interested to see how much of the network we are
% building with the current method involves the genes which are not
% expressed. (bottom 1/3)

clear
addpath('~/codes/MATLAB/myCodes/suplabel/')
addpath('~/codes/MATLAB/myCodes/general/')    

dataFolder = '/space/grp/marjan/data/';
dataFolder = '~/data/array/'
tissue = 'blood';
fileFolder = [dataFolder tissue '/matFiles/geneExprMatGemmaMapBER/outlierRemoved/'];
files = dir([fileFolder '*.mat'])

meanAbsCorrAll = zeros(length(files), 4);

for i = 1:length(files)
    i
    GSEID = getGSEIDfromStr(files(i).name)
    % load the file
    load([fileFolder files(i).name])
    
    dataSet.mat = 2.^dataSet.mat;
    gCount = size(dataSet.mat, 1);
    sampleCount = size(dataSet.mat, 2);
    
    % get the 1/3 genes, and the rest
    meanExpr = mean(dataSet.mat, 2);
    %  hist(meanExpr,40);
    q = quantile(meanExpr, 2);
    
    topGenesInd = find(meanExpr > q(1));
    botGenesInd = find(meanExpr <= q(1));
            
    sib = corr(dataSet.mat', 'type', 'Pearson');
      
    absCorr = abs(sib);
    absCorr(absCorr == 0) = nan;
    
    meanAbsCorrAll(i, 1) = mean(nanmean(absCorr));
    meanAbsCorrAll(i, 2) = mean(nanmean(absCorr(topGenesInd, ...
                                                topGenesInd)));
    meanAbsCorrAll(i, 3) = mean(nanmean(absCorr(botGenesInd, ...
                                                botGenesInd)));
    meanAbsCorrAll(i, 4) = mean(nanmean(absCorr(topGenesInd, botGenesInd)));
    
end

expID = cell(1,length(files));
for i = 1:length(files)
    expID(i) = {getGSEIDfromStr(files(i).name)}
end

h = figure; 
heatmap(meanAbsCorrAll, [{'all'}, {'top2/3'}, {'bot1/3'}, {'between'}], expID, meanAbsCorrAll, 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'hot')
'heatmap done'
title(['the mean correlation among different groups of genes, ' ...
       'tissue: ' tissue])
figFolder = ['~/data/array/' tissue '/figures/']
fileName = ['SPcorrelationInv_HM_' tissue ];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf', [figFolder fileName '.pdf']);
