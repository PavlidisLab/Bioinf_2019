% the code is to investigate the correlation among genes and it's
% correlation to expression level
% Paul's idea:
% I'm interested in the relationship between mean expression level and correlation. If you take the correlation matrix for one data set, and make distributions of correlations for the following cases:
% - between genes that are in the bottom 1/3 of expression
% - between genes that are in the top 2/3 of expression levels
% - between those two groups
% My hypothesis is that the correlations among the bottom 1/3 should be closer to zero, and the correlations between the high and low groups should be intermediate (on average).
% Also, when you get your thresholded network from a single data
% set, what iXs the distribution of links among those three
% categories (compared to expectation).
% 5. Find the correlation significance for each of the 1000
% correlation ranks, using the shuffled data. 

% 1. general variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear
addpath('~/codes/MATLAB/myCodes/suplabel/')
addpath('~/codes/MATLAB/myCodes/general/')    

%dataFolder = '/space/grp/marjan/data/';
dataFolder = '~/data/array/'
tissue = 'liver';
fileFolder = [dataFolder tissue '/matFiles/geneExprMatGemmaMapBER/'];
files = dir([fileFolder '*.mat'])

% 2. results for each dataset - 
% NOTICE! this code has a shuffling part for
% the expression level, it is commented. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for each dataset
meanExprAll = zeros(length(files), 2); % for the two group 
meanCorrAll = zeros(length(files), 4); % for the two groups and
meanPosCorrAll = zeros(length(files), 4); % taking only positive
                                             % correlations 

% this code gets the mean correlation values and does the
% correlation distribution plots
for i = 1:length(files)
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
    
    %shuffling the bottom genes expression values
    %%%%%%%%%%%%%
    % for j = 1:length(botGenesInd)
    %     currentGeneInd = botGenesInd(j);
    %     shuffledInd = datasample([1:sampleCount], sampleCount, 'Replace', false);
    %     dataSet.mat(currentGeneInd,:) = dataSet.mat(currentGeneInd, shuffledInd);
    % end
    % first plot
    % get the correlation hist for all of them and the data
        
    %sib = corr(dataSet.mat', 'type', 'Spearman');
    sib = corr(dataSet.mat', 'type', 'Pearson');
    
    %all values with positive correlation
    %%%%%%%%%%%%%
    posCorr = sib;
    posCorr(posCorr<=0) = nan;
    %getting the positive correlation values
    meanPosCorrAll(i, 1) = mean(nanmean(posCorr));
    meanPosCorrAll(i, 2) = mean(nanmean(posCorr(topGenesInd, ...
                                                  topGenesInd)));
    meanPosCorrAll(i, 3) = mean(nanmean(posCorr(botGenesInd, ...
                                                  botGenesInd)));
    meanPosCorrAll(i, 4) = mean(nanmean(posCorr(topGenesInd, ...
                                                  botGenesInd)));
        
    % ignoring the diagonal of correlation matrix
    for k = 1:gCount
        sib(k,k) = nan;
    end
       
    h = figure();
    subplot(2, 2, 1);
    hold on
    hist(sib(:), 40);
    meanCorrAll(i, 1) = nanmean(sib(:));
    title('all genes')
    s = skewness(sib(:));
    mu = nanmean(nanmean(sib));
    xlabel(sprintf('skewness=%0.3f, mean=%0.3f', s, mu));
     plot([0,0], ylim, 'w-', 'LineWidth', 1)
    plot([mu,mu],ylim,'r--','LineWidth',1)
    hold off
    
    subplot(2, 2, 2);
    holu = sib(topGenesInd, topGenesInd);
    hold on
    hist(holu(:), 40)
    meanCorrAll(i, 2) = nanmean(holu(:));
    title('top 2/3 genes')
    s = skewness(holu(:));
    mu = nanmean(nanmean(holu));
    xlabel(sprintf('skewness = %0.3f, mean=%0.3f', s,mu));
    plot([0,0], ylim, 'w-', 'LineWidth', 1)
    plot([mu,mu],ylim,'r--','LineWidth',1)
    hold off
    
    subplot(2, 2, 3);
    holu = sib(botGenesInd, botGenesInd);
    hold on
    hist(holu(:), 40)
    meanCorrAll(i, 3) = nanmean(holu(:));
    title('bottom 1/3 genes')
    s = skewness(holu(:));
    mu = nanmean(nanmean(holu));
    xlabel(sprintf('skewness = %0.3f, mean=%0.3f', s,mu));
    plot([0,0], ylim, 'w-', 'LineWidth', 1)
    plot([mu,mu],ylim,'r--','LineWidth',1)
    hold off
    
    subplot(2, 2, 4);
    holu = sib(botGenesInd, topGenesInd);
    hold on
    hist(holu(:), 40)
    meanCorrAll(i, 4) = nanmean(holu(:));
    title('bottom and top genes')
    s = skewness(holu(:));
    mu = nanmean(nanmean(holu));
    xlabel(sprintf('skewness = %0.3f, mean=%0.3f', s,mu));
    plot([mu,mu],ylim,'r--','LineWidth',1)
    plot([0,0], ylim, 'w-', 'LineWidth', 1)
    hold off
    
    
    suplabel(sprintf(['correlation distribotion in groups ' ...
                      'of genes - dataSet: %s, tissue: %s'], GSEID, ...
                     tissue), 't') 
    
    suplabel(sprintf(['red line: mean, white line: zero, ' ...
             'maximum bottom expression : %0.3f'], log2(q(1))),'x')
    
    
    fileName = sprintf('SPcorrelationInv_%s_%s', tissue, GSEID)
    figFolder = [dataFolder 'general/']
    print(h, '-depsc', [figFolder fileName '.eps']);
    print(h, '-dpdf', [figFolder fileName '.pdf']);
  
end    

% getting a heatmap of mean correlations for different groups
expID = cell(1,length(files));
for i = 1:length(files)
    expID(i) = {getGSEIDfromStr(files(i).name)}
end

h = figure; 
heatmap(meanCorrAll, [{'all'}, {'top2/3'}, {'bot1/3'}, {'between'}], expID, meanCorrAll, 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'hot')
'heatmap done'
title(['the mean correlation among different groups of genes, ' ...
       'tissue: ' tissue])
figFolder = ['~/data/array/' tissue '/figures/']
fileName = ['SPcorrelationInv_HM_' tissue ];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf', [figFolder fileName '.pdf']);


%for positive correlation values
h = figure; 
heatmap(meanPosCorrAll, [{'all'}, {'top2/3'}, {'bot1/3'}, {'between'}], expID, meanPosCorrAll, 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'hot')
'heatmap done'
title(['the mean positive correlation among different groups of genes, ' ...
       'tissue: ' tissue])
figFolder = '/space/grp/marjan/data/general/'
fileName = ['SPcorrelationInvPos_HM_' tissue ];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf', [figFolder fileName '.pdf']);

%saving the results for all these correlation values
resultFolder = '~/MATLABResults/';
file = sprintf('%sSPcorrelationInvestigation_meanCorrAll_%s.mat', ...
               resultFolder, tissue);
save(file, 'meanCorrAll')

file = sprintf('%sSPcorrelationInvestigation_meanPosCorrAll_%s.mat', ...
               resultFolder, tissue);
save(file, 'meanPosCorrAll')

%putting together all the results from tissues and do the heatmap
%for them. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
dataFolder = '~/MATLABResults/correlationInv/normalSpearman/'
files = dir(dataFolder);
load([dataFolder files(3).name])
correlationResults = meanCorrAll
for i = 4:length(files)
    load([dataFolder files(i).name])
    correlationResults = [correlationResults; meanCorrAll];
end

h = figure; 
heatmap(correlationResults, [{'all'}, {'top2/3'}, {'bot1/3'}, {'between'}], [], [], 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'hot')
'heatmap done'
title('mean correlation among different sets of genes for each dataset')
figFolder = '/space/grp/marjan/data/general/'
fileName = ['SPallCorrelationMeans'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf', [figFolder fileName '.pdf']);

% new version of shuffled one
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
addpath('~/codes/MATLAB/myCodes/suplabel/')
addpath('~/codes/MATLAB/myCodes/general/')    

dataFolder = '~/data/array/'
tissue = 'skeletalMuscle';
fileFolder = [dataFolder tissue '/matFiles/geneExprMatGemmaMapBER/'];
files = dir([fileFolder '*.mat'])

% for each dataset
meanExprAll = zeros(length(files), 2); % for the two group 
meanCorrAll = zeros(length(files), 4); % for the two groups and
meanPosCorrAll = zeros(length(files), 4); % taking only positive
                                             % correlation values
meanAbsCorrAll = zeros(length(files), 4);

%meanPosCorrAll = zeros(length(files) -2, 4); % taking only positive

% this code gets the mean correlation values and does the
% correlation distribution plots
for i = 1:length(files)

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
    
    %shuffling the genes
    %%%%%%%%%%%%%
    shuffledDataSet = zeros(gCount, sampleCount);
    for j = 1:gCount
        currentGeneInd = j;
        shuffledInd = datasample([1:sampleCount], sampleCount, 'Replace', false);
        shuffledDataSet(currentGeneInd,:) = dataSet.mat(currentGeneInd, shuffledInd);
    end
    % first plot
    % get the correlation hist for all of them and the data
        
    sib = corr(dataSet.mat', 'type', 'Pearson');
    sufSib = corr(shuffledDataSet');
    
    %all values with positive correlation
    
    posCorr = sib;
    posCorr(posCorr<=0) = nan;
    %getting the positive correlation values
    meanPosCorrAll(i, 1) = mean(nanmean(posCorr));
    meanPosCorrAll(i, 2) = mean(nanmean(posCorr(topGenesInd, ...
                                                  topGenesInd)));
    meanPosCorrAll(i, 3) = mean(nanmean(posCorr(botGenesInd, ...
                                                  botGenesInd)));
    meanPosCorrAll(i, 4) = mean(nanmean(posCorr(topGenesInd, ...
                                                  botGenesInd)));
    
    meanAbsCorr = abs(sib);
        
    % ignoring the diagonal of correlation matrix
    for k = 1:gCount
        sib(k,k) = nan;
        sufSib(k,k) = nan;
    end
           
    % h = figure();
    % hold on
    % binxos = -1:0.01:1;
    % hSib = hist(sib(:), bins);
    % hSufSib = hist(sufSib(:), bins);
    % plot(bins, hSib, 'linewidth', 2, 'color', [180, 10, 20]/255);
    % plot(bins, hSufSib, '-.', 'linewidth', 2, 'color', [180, 10, 20]/255);
    % mu = nanmean(sib(:));
    % meanCorrAll(i, 1) = mu;
    % plot([mu,mu],ylim,'LineWidth',2, 'color', [180, 10, 20]/ ...
    %      255)

    % holu = sib(topGenesInd, topGenesInd);
    % sufHolu = sufSib(topGenesInd, topGenesInd);
    % hHolu = hist(holu(:), bins);
    % hSufHolu = hist(sufHolu(:), bins);
    % plot(bins, hHolu, 'linewidth', 2, 'color', [20, 10, 180]/255);
    % plot(bins, hSufHolu, '-.', 'linewidth', 2, 'color', [20, 10, 180]/255);
    % mu = nanmean(holu(:));
    % meanCorrAll(i, 2) = mu;
    % plot([mu,mu],ylim,'LineWidth',2, 'color', [20, 10, 180]/ ...
    %      255)
        
    % holu = sib(botGenesInd, botGenesInd);
    % sufHolu = sufSib(botGenesInd, botGenesInd);
    % hHolu = hist(holu(:), bins);
    % hSufHolu = hist(sufHolu(:), bins);
    % plot(bins, hHolu, '-.','linewidth', 2, 'color', [20, 180, 10]/255);
    % plot(bins, hSufHolu, 'linewidth', 2, 'color', [20, 180, 10]/255);
    % mu = nanmean(holu(:));
    % meanCorrAll(i, 3) = mu;
    % plot([mu,mu],ylim,'LineWidth',2, 'color', [20, 180, 10]/ ...
    %      255)
    
    % holu = sib(botGenesInd, topGenesInd);
    % sufHolu = sufSib(botGenesInd, topGenesInd);
    % hHolu = hist(holu(:), bins);
    % hSufHolu = hist(sufHolu(:), bins);
    % plot(bins, hHolu, 'linewidth', 2, 'color', 'k');
    % plot(bins, hSufHolu, '-.','linewidth', 2, 'color', 'k');
    % mu = nanmean(holu(:));
    % meanCorrAll(i, 4) = mu;
    % plot([0,0], ylim, '-.','LineWidth', 2, 'color', 'k')
    % plot([mu,mu],ylim,'LineWidth',2, 'color', 'k')
    
    % hold off
  
       
    s1 = gCount * gCount;
    h = figure();
    subplot(2, 2, 1);
    hold on
    bins = -1:0.01:1;
    hSib = hist(sib(:), bins);
    hSufSib = hist(sufSib(:), bins);
    plot(bins, hSib, 'linewidth', 2, 'color', [180, 10, 20]/255);
    plot(bins, hSufSib, 'linewidth', 2, 'color', 'k');
    mu = nanmean(sib(:));
    meanCorrAll(i, 1) = mu;
    title('all genes')
    s = skewness(sib(:));
    xlabel(sprintf('skewness=%0.3f, mean=%0.3f', s, mu));
    plot([0,0], ylim, 'k-', 'LineWidth', 2)
    plot([mu,mu],ylim,'--','LineWidth',2, 'color', [180, 10, 20]/ ...
         255)
    
    s2 = length(topGenesInd)*length(topGenesInd)
    subplot(2, 2, 2);
    holu = sib(topGenesInd, topGenesInd);
    sufHolu = sufSib(topGenesInd, topGenesInd);
    hold on
    hHolu = hist(holu(:), bins);
    hSufHolu = hist(sufHolu(:), bins);
    plot(bins, hHolu, 'linewidth', 2, 'color', [180,10,20]/255);
    plot(bins, hSufHolu, 'linewidth', 2, 'color', 'k');
    mu = nanmean(holu(:));
    meanCorrAll(i, 2) = mu;
    title('top 2/3 genes')
    s = skewness(holu(:));
    xlabel(sprintf('skewness = %0.3f, mean=%0.3f', s,mu));
    plot([0,0], ylim, 'k-', 'LineWidth', 2)
    plot([mu,mu],ylim,'--','LineWidth',2,'color', [180,10,20]/255);
    hold off
    
    s3 = length(botGenesInd)*length(botGenesInd)
    subplot(2, 2, 3);
    holu = sib(botGenesInd, botGenesInd);
    sufHolu = sufSib(botGenesInd, botGenesInd);
    hold on
    hHolu = hist(holu(:), bins);
    hSufHolu = hist(sufHolu(:), bins);
    plot(bins, hHolu, 'linewidth', 2, 'color', [180,10,20]/255);
    plot(bins, hSufHolu, 'linewidth', 2, 'color', 'k');
    mu = nanmean(holu(:));
    meanCorrAll(i, 3) = mu;
    title('bottom 1/3 genes')
    s = skewness(holu(:));
    xlabel(sprintf('skewness = %0.3f, mean=%0.3f', s,mu));
    plot([0,0], ylim, 'k-', 'LineWidth', 2)
    plot([mu,mu],ylim,'--','LineWidth',2, 'color', [180,10,20]/255)
    hold off
    
    subplot(2, 2, 4);
    holu = sib(botGenesInd, topGenesInd);
    sufHolu = sufSib(botGenesInd, topGenesInd);
    hold on
    hHolu = hist(holu(:), bins);
    hSufHolu = hist(sufHolu(:), bins);
    plot(bins, hHolu, 'linewidth', 2, 'color', [180,10,20]/255);
    plot(bins, hSufHolu, 'linewidth', 2, 'color', 'k');
    mu = nanmean(holu(:));
    meanCorrAll(i, 4) = mu;
    title('bottom and top genes')
    s = skewness(holu(:));
    xlabel(sprintf('skewness = %0.3f, mean=%0.3f', s,mu));
    plot([0,0], ylim, 'k-', 'LineWidth', 2)
    plot([mu,mu],ylim,'--','LineWidth',2, 'color', [180,10,20]/255)
    hold off
    
    
    
    suplabel(sprintf(['correlation distribotion in groups ' ...
                      'of genes - dataSet: %s, tissue: %s'], GSEID, ...
                     tissue), 't') 
    
    suplabel(sprintf(['black : shuffled, red: not shuffled, ' ...
             'maximum bottom expression : %0.3f'], log2(q(1))),'x')
    
        
        
    fileName = sprintf('correlationInv_%s_%s', tissue, GSEID)
    figFolder = [dataFolder tissue '/figures/']
    print(h, '-depsc', [figFolder fileName '.eps']);
    print(h, '-dpdf', [figFolder fileName '.pdf']);

  
end    


expID = cell(1,length(files));
for i = 1:length(files)
    expID(i) = {getGSEIDfromStr(files(i).name)}
end

h = figure; 
heatmap(meanCorrAll, [{'all'}, {'top2/3'}, {'bot1/3'}, {'between'}], expID, meanCorrAll, 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'hot')
'heatmap done'
title(['the mean correlation among different groups of genes, ' ...
       'tissue: ' tissue])
figFolder = ['~/data/array/' tissue '/figures/']
fileName = ['SPcorrelationInv_HM_' tissue ];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf', [figFolder fileName '.pdf']);


% the network section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getting the true networks
% for each network, we want to keep : 
% 1. topMeanND,  
% 2. botMeanND,
% 3. allMeanND, 
% 4. conEdges, 
% 5. expectedConEdges : 2/3*botMeanND*1/3*N
% where N is the node count

clear
addpath('~/codes/MATLAB/myCodes/suplabel/')
addpath('~/codes/MATLAB/myCodes/general/')    

dataFolder = '/space/grp/marjan/data/';
dataFolder = '~/data/array/'
tissue = 'blood';
fileFolder = [dataFolder tissue '/matFiles/geneExprMatGemmaMapBER/outlierRemoved/'];
files = dir([fileFolder '*.mat'])

results = zeros(length(files), 5);

for i = 1:length(files)
    
    GSEID = getGSEIDfromStr(files(i).name)
    % load the file
    load([fileFolder files(i).name])
    
    dataSet.mat = 2.^ dataSet.mat;
    gCount = size(dataSet.mat,1);
    % get the 1/3 genes, and the rest
    meanExpr = mean(dataSet.mat, 2);
    %  hist(meanExpr,40);
    q = quantile(meanExpr, 2);
    
    topGenesInd = find(meanExpr > q(1));
    botGenesInd = find(meanExpr <= q(1));
    
    % first plot
    % get the correlation hist for all of them and the data
    sib = corr(dataSet.mat');
    sib = sib - eye(gCount);

    upperSingle = triu(sib, 1);
    QSingle = quantile(upperSingle(:), (1 - 0.005))
    singleNet = sib > QSingle;
    ND = sum(singleNet);
    
    topNet = singleNet(topGenesInd, topGenesInd);
    botNet = singleNet(botGenesInd, botGenesInd);
    conNet = singleNet(topGenesInd, botGenesInd);
    
    topNetND = sum(topNet);
    botNetND = sum(botNet);
    conCount = sum(sum(conNet)); % for connection I only care about
                                 % number of edges
    
    topND = ND(topGenesInd);
    botND = ND(botGenesInd);
    
    topMeanND = mean(topND);
    botMeanND = mean(botND);
    allMeanND = mean(ND);
    
    results(i, 1) = topMeanND;
    results(i, 2) = botMeanND;
    results(i, 3) = allMeanND;
    results(i, 4) = conCount;
    results(i, 5) = botMeanND * (2/3) * (length(botND));
        
    q = quantile(ND, 2);
    lowNDInd = find(ND <= q(1));
    p = intersect(lowNDInd, botGenesInd);
    p = length(p)/length(botGenesInd);
    e = sum(botND)/sum(ND);
    %the three groups ND distribotion in the whole network
    h = figure;
    subplot(1, 3,1)
    hist(ND, 40);
    mu = allMeanND;
    xlabel(sprintf('mean=%0.3f', mu));
    title('all genes')
    
    subplot(1,3,2)
    hist(topND, 40);
    mu = topMeanND;
    xlabel(sprintf('mean=%0.3f', mu));
    title('top 2/3 genes')
    
    subplot(1,3,3)
    hist(botND, 40);
    mu = botMeanND;
    xlabel(sprintf('mean=%0.3f', mu));
    title('bottom 1/3 genes')
    
    suplabel(sprintf(['node degree distribotion in groups ' ...
                      'of genes - dataSet: %s, tissue: %s, min correlation: %0.3f'], GSEID, ...
                     tissue, QSingle), 't') 
    
    suplabel(sprintf(['%0.3f of the 1/3 genes with lowest node ' ...
                      'degree were also in the bottom 1/3 \n '...
                      '%0.3f of the edges were in the bottom 1/3 ' ...
                      'network'], p, e) ,'x')
    
    fileName = sprintf('correlationInv_NDdist_%s_%s', tissue, GSEID)
    figFolder = [dataFolder 'general/']
    print(h, '-depsc', [figFolder fileName '.eps']);
    print(h, '-dpdf', [figFolder fileName '.pdf']);
    
end

%putting up the two heatmpas for the result

h = figure;

expID = cell(1,length(files));
for i = 1:length(files)
    expID(i) = {getGSEIDfromStr(files(i).name)}
end

subplot(1, 2, 1)

heatmap(results(:, 1:3), [{'top mean ND'}, {'bottom mean ND'}, {'all mean ND'}], expID, [], 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'summer')
title('the mean node degree')

subplot(1,2,2)

heatmap(results(:, 4:5), [{'true connection edges'}, {'expected connection edges'}], expID, [], 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'summer')
title('connection edges count')

suplabel(['tissue: ' tissue], 't')

figFolder = '/space/grp/marjan/data/general/'
fileName = ['NLcorrelationInvNDdist_HM_' tissue ];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf', [figFolder fileName '.pdf']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% keeping the index of bottom 1/3 and see how they overlap in one
% tissue


clear
addpath('~/codes/MATLAB/myCodes/suplabel/')
addpath('~/codes/MATLAB/myCodes/general/')    

dataFolder = '/space/grp/marjan/data/';
dataFolder = '~/data/array/'
tissue = 'liver';
fileFolder = [dataFolder tissue '/matFiles/geneExprMatGemmaMapBER/outlierRemoved/'];
files = dir([fileFolder '*.mat'])

% for each dataset
botCount = 6126;
topCount = 18377 - botCount;
botInd = zeros(botCount, length(files));
topInd = zeros(topCount, length(files));

%meanPosCorrAll = zeros(length(files) -2, 4); % taking only positive

% this code gets the mean correlation values and does the
% correlation distribution plots
for i = 1:length(files)

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
    topInd(:, i) = topGenesInd;
    botInd(:, i) = botGenesInd;

end

blood.botInd = botInd;
blood.topInd = topInd;

lung.botInd = botInd;
lung.topInd = topInd;

liver.botInd = botInd;
liver.topInd = topInd;

sm.botInd = botInd;
sm.topInd = topInd;

totalBotInd = [blood.botInd lung.botInd liver.botInd sm.botInd];
interBotInd = zeros(44, 44);
for i = 1:size(totalBotInd, 2)
    for j = i+1:size(totalBotInd, 2)
        interBotInd(i,j) = length(intersect(totalBotInd(:,i), ...
                                     totalBotInd(:,j)));
        interBotInd(j,i) = interBotInd(i,j);
    end
    interBotInd(i, i) = size(totalBotInd, 1);
end

interBotInd = interBotInd/max(max(interBotInd));

h = figure;
heatmap(interBotInd, [], [], [], 'Colorbar', true, 'Colormap', 'hot')

print(h, '-dpdf', 'bottomeGenes.pdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clustering datasets with their bottom genes in heatmap
% DID NOT WORK. DONE

% I got the function myDist which returns what the pdist function
% wants, and this is how linkage work with it

book = linkage(totalBotInd', 'average',  '@myDist');
dendrogram(book, 50)

CGobj = clustergram(totalBotInd', 'Cluster', 2, 'Standardize', 'row','RowPDist', '@myDist')

test = clustergram(dataSet.mat(1:100, :), 'Standardize', 'Row')

load filteredyeastdata
clustergram(yeastvalues(1:10, :))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem : filtering the genes or , which genes should we filter. 
% proposed answer : the bot 1/3, since their expression value is
% supposed to be noise, so their correlation with genes would be
% more closer to zero. 

% that is not the case. it seems that the lower expressed genes
% have their own correlation distribution, as in most datasets
% their mean correlation is bigger than the mean correlation
% between the top genes(test?) However, the absolute mean correlation is
% ofcourse lower. So they are in the network, but we don't know
% how they affect the network. 

% should we filter these genes
% are there any differences between these genes for different
% tissues? 

% for each tissue,
% for each dataset : make and save individual networks. 
% keep the bot and top indexes, and how many edges are from those
% networks. 
% later, how many of those edges find their way to the final
% network. 

h = figure
hold all

for i = 1:100
    bins = -1:0.01:1;
    %    htot = hist(sib(botGenesInd(i),:), bins);
    htop = hist(sib(botGenesInd(i), topGenesInd), bins);
    hbot = hist(sib(botGenesInd(i), botGenesInd), bins);
    %    plot(bins, htot, 'linewidth', 0.5, 'color', 'k')
    plot(bins, htop, '.', 'linewidth', 1, 'color', 'r')
    plot(bins, hbot, '.', 'linewidth', 1, 'color', 'b')
    
    htop = hist(sib(topGenesInd(i), topGenesInd), bins);
    hbot = hist(sib(topGenesInd(i), botGenesInd), bins);
    %    plot(bins, htot, 'linewidth', 0.5, 'color', 'k')
    plot(bins, htop, '.', 'linewidth', 1, 'color', 'g')
    plot(bins, hbot, '.', 'linewidth', 1, 'color', 'k')
    
end

% section : 4? sample correlation distribution of shuffled and real data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
addpath('~/codes/MATLAB/myCodes/suplabel/')
addpath('~/codes/MATLAB/myCodes/general/')    

dataFolder = '~/data/array/'
tissue = 'skeletalMuscle';
fileFolder = [dataFolder tissue '/matFiles/geneExprMatGemmaMapBER/'];
files = dir([fileFolder '*.mat'])
resFolder = ['~/data/array/' tissue '/matFiles/thresholdResults/']

res = zeros(length(files), 3, 15);

% this code gets the mean correlation values and does the
% correlation distribution plots
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
    
    %shuffling the genes
    %%%%%%%%%%%%%
    tic
    shuffledDataSet = zeros(gCount, sampleCount);
    for j = 1:gCount
        currentGeneInd = j;
        shuffledInd = datasample([1:sampleCount], sampleCount, 'Replace', false);
        shuffledDataSet(currentGeneInd,:) = dataSet.mat(currentGeneInd, shuffledInd);
    end
    toc
    % first plot
    % get the correlation hist for all of them and the data
    
    tic
    sib = corr(dataSet.mat', 'type', 'Pearson');
    sufSib = corr(shuffledDataSet');
    toc
    %all values with positive correlation
    
    % ignoring the diagonal of correlation matrix
    for k = 1:gCount
        sib(k,k) = nan;
        sufSib(k,k) = nan;
    end

    for k = 1: 15
        k
        thr = 0.004 + (k * 0.001)
        upperSingle = triu(sib, 1);
        QSingle = quantile(upperSingle(:), (1 - thr));
        sufUpperSingle = triu(sufSib, 1);
        res(i, 1, k) = thr;
        res(i, 2, k) = QSingle;
        res(i, 3, k) = sum(sufUpperSingle(:) >= QSingle);
    end    
    
    % h = figure();
    % hold all
    % bins = -1:0.01:1;
    % hSib = hist(sib(:), bins);
    % totalC = sum(hSib);
    % hSib = hSib/totalC;
    % hSufSib = hist(sufSib(:), bins);
    % hSufSib = hSufSib/totalC;
    % plot(bins, hSib, 'linewidth', 2, 'color', 'k');
    % plot(bins, hSufSib, '-','linewidth', 2, 'color', 'r');
    
    % for each network, how does the shuffle and the not shuffle
    % values change
    
    % holu = sib(topGenesInd, topGenesInd);
    % sufHolu = sufSib(topGenesInd, topGenesInd);
    % hSib = hist(holu(:), bins);
    % hSib = hSib/totalC;
    % hSufSib = hist(sufHolu(:), bins);
    % hSufSib = hSufSib/totalC;
    % plot(bins, hSib, 'linewidth', 2, 'color', 'b');
    % plot(bins, hSufSib, '-','linewidth', 2, 'color', 'b');
    

    % holu = sib(botGenesInd, botGenesInd);
    % sufHolu = sufSib(botGenesInd, botGenesInd);
    % hSib = hist(holu(:), bins);
    % hSib = hSib/totalC;
    % hSufSib = hist(sufHolu(:), bins);
    % hSufSib = hSufSib/totalC;
    % plot(bins, hSib, 'linewidth', 2, 'color', 'g');
    % plot(bins, hSufSib, '-','linewidth', 2, 'color', 'g');


    % holu = sib(botGenesInd, topGenesInd);
    % sufHolu = sufSib(botGenesInd, topGenesInd);
    % hSib = hist(holu(:), bins);
    % hSib = hSib/totalC;
    % hSufSib = hist(sufHolu(:), bins);
    % hSufSib = hSufSib/totalC;
    % plot(bins, hSib, 'linewidth', 2, 'color', 'r');
    % plot(bins, hSufSib, '-','linewidth', 2, 'color', 'r');

    %  fileName = sprintf('correlationInv_%s_%s', tissue, GSEID)
    % figFolder = [dataFolder tissue '/figures/']
    % print(h, '-depsc', [figFolder fileName '.eps']);
    % print(h, '-dpdf', [figFolder fileName '.pdf']);

    
end

fileName = [resFolder  '_thresholdInvRes.mat'];
save(fileName, 'res')

%%%%%%
load([])

% 5. Find the correlation significance 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 

addpath('~/codes/MATLAB/myCodes/suplabel/')
addpath('~/codes/MATLAB/myCodes/general/')    

% path the files 

tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}

% for each dataset (in each tissue)
d = 1
for t = 1:5
    d
    % read the dataset files
    dataFolder = '~/data/affyArray/tissues/'
    tissue = tissues{t};
    fileFolder = [dataFolder tissue '/matFiles/geneExprMatGemmaMapBER/'];
    files = dir([fileFolder '*.mat'])

    % for each dataset
    for i = 1:length(files)
        i

        GSEID = getGSEIDfromStr(files(i).name)
        % load the file
        load([fileFolder files(i).name])
        
        gCount = size(dataSet.mat, 1)
        sampleCount = size(dataSet.mat, 2)
        
        % get the 1/3 genes, and the rest
        meanExpr = mean(dataSet.mat, 2);
        %  hist(meanExpr,40);
        q = quantile(meanExpr, 2);
        
        topGenesInd = find(meanExpr > q(1));
        botGenesInd = find(meanExpr <= q(1));
        
        %shuffling the genes
        %%%%%%%%%%%%%
        tic
        shuffledDataSet = zeros(gCount, sampleCount);
        for j = 1:gCount
            currentGeneInd = j;
            shuffledInd = datasample([1:sampleCount], sampleCount, 'Replace', false);
            shuffledDataSet(currentGeneInd,:) = dataSet.mat(currentGeneInd, shuffledInd);
        end
        toc

        % get the correlation hist for all of them and the data
        tic
        sib = corr(dataSet.mat');
        sufSib = corr(shuffledDataSet');
        toc

        %getting the file's network
        upperTrio = sib(logical(triu(ones(size(sib)), 1)));
        upperShTrio = sufSib(logical(triu(ones(size(sufSib)), 1)));
        Qbase = quantile(upperTrio, [1:10000]/10000);
        Qshuff = quantile(upperShTrio, [1:10000]/10000);
        
        tic
        shuffLcount = zeros(1, 1000);
        for j = 1:1000
            shuffCount(1000 - j +1) = sum(upperShTrio >= Qbase(10000 - j + 1));
        end
        toc

        % save teh Qs for later! 
        DScorrSig(d).name = GSEID;
        DScorrSig(d).tissue = tissue;
        DScorrSig(d).sampleCount = sampleCount;
        DScorrSig(d).Qbase = Qbase;
        DScorrSig(d).Qshuff = Qshuff;
        DScorrSig(d).shuffCount = shuffCount;
        DScorrSig(d).description = ['shuffCount is count of shuffled data' ...
                            'correlations for the top 1000 correlation ' ...
                            'bins thresholds']
        d = d  + 1;
        % for saving: dsname, tissue, sample count, the
        % significance thr, the correlaion value at each
        % significance thr. 
        
    end
end

save('~/data/general/shuffledSignificanceCorrResults.mat', 'DScorrSig')



