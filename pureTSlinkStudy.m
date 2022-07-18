% This is to look at the expression level of the coexpressed links,
% make some plots and check the reproducibility of the links. 
% 
% 1. get the count of pure links as those expressed in the
% datasets, compare their TS distribution with others and make dist
% plots
% 
% 10. testing the expression level for the TS links and finding
% examples

% 11. comparing the pValues vs the distances

% 11.5 getting the TS links - obsolete

% 11.75. getting the TS links - new

% 12. new pure link study: getting the everywhere expressed genes
% using ANOVA

% 12.5 doing the regression for different tissue expression

% 13. Getting the TS genes using ANOVA

% 15. get the final table
%15.1 get the genral info for ts
% 15.2 get the count of pure and other ts links

% 16: get the expression level for the regression identified pure
% links. get the ATN for those genes. 

% 17: get the ATN and TSN for the top 1/3 links.

% Draft 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

addpath('~/codes/MATLAB/myCodes/general/')
load([ '~/data/general/GPL570GemmaMapNEW.mat'])

% BIG FAT NOTE : I am losing 20-30% of the links based on the
% expression filter, do I care? I can use the .6 expression
% filter, I think it is worth it. 

% getting the TS networks
defThr = .6
thr = defThr;

load(['~/networks/continuousTSLinksResults/' ...
      'brainSpaceMat_negAdjusted.mat'])
load(['~/networks/continuousTSLinksResults/random21_smallMat_negAdjusted.mat'])
brainTS = spaceMat >thr;
sum(sum(brainTS))
brainTS = brainTS .* spaceMat;

load('~/networks/continuousTSLinksResults/bloodSpaceMat.mat')
bloodTS = spaceMat >thr;
sum(sum(bloodTS))
bloodTS = bloodTS .* spaceMat;

load('~/networks/continuousTSLinksResults/liverSpaceMat.mat')
liverTS = spaceMat >thr;
sum(sum(liverTS))
liverTS = liverTS .* spaceMat;
load('~/networks/continuousTSLinksResults/lungSpaceMat.mat')
lungTS = spaceMat >thr;
sum(sum(lungTS))
lungTS = lungTS .* spaceMat;

load('~/networks/continuousTSLinksResults/skeletalMuscleSpaceMat.mat')
muscleTS = spaceMat >thr;
sum(sum(muscleTS))
muscleTS = muscleTS .* spaceMat;

%%%%%%%%%%%%%%%%%%%%%
clear
addpath('~/codes/MATLAB/myCodes/general/')
load([ '~/data/general/GPL570GemmaMapNEW.mat'])

tissues = {'brain', 'blood', 'liver', 'lung', 'skeletalMuscle'}

% get their expression levels, 
%load('~/data/general/linkExprInfo/wholeGeneExpr.mat')
%load('~/data/general/linkExprInfo/dataSetProbeInf.mat')
load('~/data/general/linkExprInfo/dataSetExpInf.mat', 'dataSetExpInf')
%load('~/data/general/linkExprInfo/totalMeanExp.mat', ...
%    'totalMeanExp')
load('~/data/general/linkExprInfo/totalMeanExpCentered.mat', ...
     'totalMeanExpCentered')

%%% NOTE: REDO THIS PART WITH ANOVA FOR ALL GENES AT THE SAME TIME: PKADO
% this gets the genes for each tissue 
defThr = .6
for t = 1:5 
    tissue = tissues{t};
    % load(['~/networks/tissues/' tissue '/' ...
    %       'binaryNet_0.5DSCountCorr_0.8Expr.mat'])
    %tempNet/tissue Networks
 
    load(['~/networks/tissues/' tissue '/' ...
         'binaryNet_4QDSCounter_0.8Expr_Ind0.05.mat'])
   
    
    load(['~/networks/continuousTSLinksResults/' tissue ...
          'SpaceMat_negAdjusted.mat']) % spaceMat

    totalTissueNet = spaceMat .* binNet;
    sum(sum(totalTissueNet > defThr))

    sib = pValSelected + negAdjSelected;
    
    % get the unique genes for the TS links
    sum(sum(totalTissueNet > defThr))
    
    TSNet = totalTissueNet >defThr;
    sum(sum(TSNet))
    TSNet = TSNet .* spaceMat;
    
    book = sum(TSNet) + sum(TSNet');
    uniqueGenes = book > 0;
    uniqueInd = find(book > 0);
    uniqueGenesIDs = gpl570.uniqueSymbols(uniqueGenes);
    sum(uniqueGenes)

    % get the pure links... 
    
    % Just a visualization for the links.
    % clustering the genes for pure links 
    uniqueGenesExpr = totalMeanExpCentered(uniqueGenes, :);
    Y = pdist(uniqueGenesExpr);
    Z = linkage(Y, 'complete');
    clusterCount = 10
    T = cluster(Z, 'maxclust', clusterCount);
    
    % sorting based on the cluster
    [a, b] = sort(T);
    sortedExpr = uniqueGenesExpr(b, :);
    sortedIDs = uniqueGenesIDs(b);
    sortedInd = uniqueInd(b);
    
    % h = figure; 
    % heatmap(sortedExpr, [], [],  [], 'TickAngle', 45, ...
    %         'Colorbar', true, 'colormap', bone)
    
    % title(sprintf('expression level of %d genes from %s TS links', ...
    %               length(sortedExpr), tissue))
    
    % fileName = [tissue 'ExprLevel_TSlinks_HM04_BF']
    % figFolder = ['~/data/affyArray/figures/TSplots/']
    % print(h, '-depsc', [figFolder fileName '.eps']);
    % print(h, '-dpdf', [figFolder fileName '.pdf'])

    % getting the genes with similar centered expression level
    % amongst all the datasets. 
    p = zeros(1, size(sortedExpr, 1));
    for j = 1:size(sortedExpr,1)
        p(j) = anova1(sortedExpr(j,:), {dataSetExpInf(:).tissue}, 'off' );
    end
    
    %just a test:
    % pkado = zoers(1, 18494);
    % for j = 1:18494
    %    pkado(j) = anova1(totalMeanExpCentered, {dataSetExpInf(:).tissue}, 'off');
    % end
    
    eqGenes = p > .00005/(length(sortedExpr));
    sortedEQExpr = sortedExpr(eqGenes, :);
    sortedEQInd = sortedInd(eqGenes);
    sortedEQIDs = sortedIDs(eqGenes);
    %    anova1(sortedEQExpr(2,:), {dataSetExpInf(:).tissue});
    length(sortedEQInd)
    
    linkGenesInfo(t).tissue = tissue;
    linkGenesInfo(t).pValues = p;
    linkGenesInfo(t).geneInds = sortedInd;
    linkGenesInfo(t).geneIDs = sortedIDs;
    linkGenesInfo(t).exp = sortedExpr;

    % h = figure; 
    % heatmap(sortedEQExpr(:, :), [], [],  [], 'TickAngle', 45, ...
    %         'Colorbar', true, 'colormap', bone)
    %     title(sprintf('expression level of %d genes from %s TS links', ...
    %               length(sortedEQExpr), tissue))
    % eqGeneTemplate = zeros(18494, 18494);
    % eqGeneTemplate(sortedEQInd, sortedEQInd) = 1;
    % eqGeneTemplate = sparse(eqGeneTemplate);
    % myTSEQlinks = TSNet .* eqGeneTemplate;
    % tissue
    % count = full(sum(sum(myTSEQlinks>0)))
    % str = sprintf('%d links between these genes', count)
    % consistentLinks.tissue = tissue;
    % consistentLinks.links = myTSEQlinks;
    % text(5,length(sortedEQExpr)/10, str, 'color', 'r' , 'FontSize', 14)    
    % fileName = [tissue 'ExprLevel_TSlinks_consistent_HM04_BF']
    % figFolder = ['~/data/affyArray/figures/TSplots/']
    % print(h, '-depsc', [figFolder fileName '.eps']);
    % print(h, '-dpdf', [figFolder fileName '.pdf'])
end

% from this, we have genes which are expressed at a steady level
% between the tissues. 

% eqGenes = p > .005;
% sortedEQExpr = sortedExpr(eqGenes, :);
% sortedEQInd = sortedInd(eqGenes);
% sortedEQIDs = sortedIDs(eqGenes);
% anova1(sortedEQGenes(2,:), {dataSetExpInf(:).tissue})
% length(sortedEQInd)

% h = figure; 
% heatmap(sortedEQExpr(:, :), [], [],  [], 'TickAngle', 45, ...
%         'Colorbar', true, 'colormap', bone)

% % save the plots for that 

% % get the links between them 
% eqGeneTemplate = zeros(18494, 18494);
% eqGeneTemplate(sortedEQInd, sortedEQInd) = 1;
% eqGeneTemplate = sparse(eqGeneTemplate);
% myTSEQlinks = TSNet .* eqGeneTemplate;
% myTSEQlinks

%% >>>>>> Results
% for genes which are consistantly expressed amongst all the tissues
whos linkGenesInfo

tissues = {linkGenesInfo(:).tissue};
for t = 1:length(tissues)
    tissue = tissues{t};
    tissue
    tissueCell = {dataSetExpInf(:).tissue};
    dsIndCell = strfind(tissueCell, tissue);
    allDSInd = find(cellfun(@isempty, dsIndCell));
    labels = cell(length(tissueCell),1);
    labels(:) = {'that'};
    labels(allDSInd) = {'this'};
    indexing = ones(size(dsIndCell));
    allDSInd
    indexing(allDSInd) = 0;
    indexing = logical(indexing);
    indexing
    p = zeros(size(linkGenesInfo(t).pValues));
    fChange = zeros(size(linkGenesInfo(t).pValues));
    linkGenesInfo(t).exp = linkGenesInfo(t).exp + abs(min(min(linkGenesInfo(t).exp)))
    for j = 1: length(linkGenesInfo(t).pValues)
        exp = linkGenesInfo(t).exp(j, :);
        p(j) = anova1(exp, labels , 'off' );
        fChange(j) =  (mean(exp(indexing)) + 1e-3)/ (mean(exp(~indexing)) + ...
                                            1e-3);
        lfc(j) = log((mean(exp(indexing))+1e-3)/(mean(exp(~indexing)) + ...
                                            1e-3));
    end
    linkGenesInfo(t).allVSTissuePval = p;
    linkGenesInfo(t).foldChange = fChange;
    linkGenesInfo(t).lfc = lfc;
end

%save(sprintf('~/networks/tissues/linkGenesInfo_TSThr%.2f.mat', defThr), 'linkGenesInfo')
save(sprintf('~/networks/tissues/linkGenesInfo_Thr%.2f_4QDS_0.8Expr_Ind0.05.mat', ...
             defThr), 'linkGenesInfo')

clear
load('~/networks/tissues/linkGenesInfo.mat')
defThr = .35
%load(sprintf('~/networks/tissues/linkGenesInfo_TSThr%.2f.mat', defThr))
load(sprintf('~/networks/tissues/linkGenesInfo_Thr%.2f_4QDS_0.8Expr_Ind0.05.mat', ...
             defThr))
% so I got two p-values, now the numbers from them. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for each tissue, i can get the network 

% p-value for the Bonferrani correction
pThr01 = 5e-4
pThr02 = 5e-4
tsThr = defThr
for t = 1:length(linkGenesInfo)
    t
    this = linkGenesInfo(t);
    tissue = this.tissue;
    
    % load(['~/networks/tissues/' tissue '/' ...
    %       'binaryNet_0.5DSCountCorr_0.8Expr.mat'])
    % %tempNet/tissue binary Networks
    
        load(['~/networks/tissues/' tissue '/' ...
         'binaryNet_4QDSCounter_0.8Expr_Ind0.05.mat'])
        %binNet/ tissueBinaryNetwork

    load(['~/networks/continuousTSLinksResults/' tissue ...
          'SpaceMat.mat']) % spaceMat

    % I am just going to keep the network for the later, cause I
    % already have the gene list
    totalTissueNet = spaceMat .* binNet;
    tsNet = totalTissueNet > tsThr;
    tsNet = tsNet .*spaceMat;
    sum(sum(tsNet))

    % getting the links between the equally expressed genes
    p01 = (this.pValues > 5e-4/(length(this.pValues)));
    fc01 = (this.foldChange <=1.3);
    eqGenes = (fc01 + p01)  == 2;
    % sortedEQExpr = sortedExpr(eqGenes, :);
    eqGenesInd = this.geneInds(eqGenes);
    % sortedEQIDs = sortedIDs(eqGenes);
    % length(sortedEQInd)
    'eqGenes count'
    length(eqGenesInd)
    
    eqGeneTemplate = zeros(18494, 18494);
    eqGeneTemplate(eqGenesInd, eqGenesInd) = 1;
    eqGeneTemplate = sparse(eqGeneTemplate);
    myTSEQlinks = tsNet .* eqGeneTemplate;
    eqGeneLinkCount = full(sum(sum(myTSEQlinks>0)))
    
    % getting the links between between tissue specific links
    fcInd = this.foldChange >= 1.3;
    pInd = this.allVSTissuePval < (pThr02/ ...
                                        length(this.foldChange));
    
    tsGenes = (pInd + fcInd) == 2;
    %tsGenesInd = this.geneInds(tsGenes);
    %    tsGenes = pInd + fcInd;
    %tsGenes = tsGenes == 2;
    tsGenesInd = this.geneInds(tsGenes);
    
    'tsGenes Count'
    length(tsGenesInd)
    
    tsGeneTemplate = zeros(18494, 18494);
    tsGeneTemplate(tsGenesInd, tsGenesInd) = 1;
    tsGeneTemplate = sparse(tsGeneTemplate);
    myTSTSlinks = tsNet .* tsGeneTemplate;
    tsGeneLinkCount = full(sum(sum(myTSTSlinks>0)))
    
    othGenes = ~((tsGenes + eqGenes) == 1);
    othGenesInd = this.geneInds(othGenes);
    
    othGeneTemplate = zeros(18494, 18494);
    othGeneTemplate(othGenesInd, othGenesInd) = 1;
    othGeneTemplate = sparse(othGeneTemplate);
    myOTHTSlinks = tsNet .* othGeneTemplate;
    othGeneLinkCount = full(sum(sum(myOTHTSlinks > 0)))
    
    % getting the links between TS and common genes
    tcTemplate = zeros(18494, 18494);
    tcTemplate(tsGenesInd, eqGenesInd) = 1;
    tcTemplate(eqGenesInd, tsGenesInd) = 1;
    tcTemplate = sparse(tcTemplate);
    mytcTSlinks = tsNet .* tcTemplate;
    tcGeneLinkCount = full(sum(sum(mytcTSlinks > 0)))

    % getting the links betweeen TS and other genes
    toTemplate = zeros(18494, 18494);
    toTemplate(othGenesInd, tsGenesInd) = 1;
    toTemplate(tsGenesInd, othGenesInd) = 1;
    toTemplate = sparse(toTemplate);
    mytoTSlinks = tsNet .*toTemplate;
    toGeneLinkCount = full(sum(sum(mytoTSlinks > 0)))
    
    % getting the links betwen other and common geens
    ocTemplate = zeros(18494, 18494);
    ocTemplate(othGenesInd, eqGenesInd) = 1;
    ocTemplate(eqGenesInd, othGenesInd) = 1;
    ocTemplate = sparse(ocTemplate);
    myocTSlinks = tsNet .*ocTemplate;
    ocGeneLinkCount = full(sum(sum(myocTSlinks > 0)))
    
    ocGeneLinkCount + toGeneLinkCount + tcGeneLinkCount + othGeneLinkCount ...
        + eqGeneLinkCount + tsGeneLinkCount
    
    testNet = tsNet + myTSTSlinks + myTSEQlinks + myOTHTSlinks + ...
            myocTSlinks + mytoTSlinks + mytcTSlinks;
    
    sum(tsGenes) + sum(othGenes) + sum(eqGenes)
    
    finalTable(t).wholeNet = tsNet;
    finalTable(t).tissue = tissue;
    finalTable(t).cGenePval = pThr01;
    finalTable(t).tsGenePval = pThr02;
    finalTable(t).eqGenesInd = eqGenesInd;
    finalTable(t).tsGenesInd = tsGenesInd;
    finalTable(t).othGenesInd = othGenesInd;
    finalTable(t).othGenesCount = length(othGenesInd);
    finalTable(t).cGenesCount = length(eqGenesInd);
    finalTable(t).tsGenesCount = length(tsGenesInd);
    finalTable(t).tsEdges = myTSTSlinks;
    finalTable(t).pureEdges = myTSEQlinks;
    finalTable(t).othEdges = myOTHTSlinks;
    finalTable(t).tcEdges = mytcTSlinks;
    finalTable(t).toEdges = mytoTSlinks;
    finalTable(t).coEdges = myocTSlinks;
    finalTable(t).tsEdgeCount = tsGeneLinkCount;
    finalTable(t).pureEdgeCount = eqGeneLinkCount;
    finalTable(t).othEdgeCount = othGeneLinkCount;
    finalTable(t).tcEdgeCount = tcGeneLinkCount;
    finalTable(t).toEdgeCount = toGeneLinkCount;
    finalTable(t).coEdgeCount = ocGeneLinkCount;
end

save(sprintf('~/data/general/TSlinksTable_4QDSCounter_0.8Expr_Ind0.05%.2f_TSValueSaved.2f.mat', defThr), ...
     'finalTable')


save(sprintf('~/data/general/TSlinksTable_4QDSCounter_0.8Expr_Ind0.05%.2f.mat', defThr), ...
     'finalTable')

clear
load([ '~/data/general/GPL570GemmaMapNEW.mat'])
defThr = .58
load(sprintf('~/data/general/TSlinksTable_4QDSCounter_0.8Expr_Ind0.05%.2f_TSValueSaved.2f.mat', defThr))

load(sprintf('~/data/general/TSlinksTable%.2f.mat', defThr))
load('~/data/general/TSlinksTable.mat') % NOTE: default is .63

myexp = this.exp(b, :);
h = figure; 
heatmap(sortedEQExpr(:, :), [], [],  [], 'TickAngle', 45, ...
        'Colorbar', true, 'colormap', bone)
title(sprintf('expression level of %d genes from %s TS links', ...
              length(sortedEQExpr), tissue))

brainFC = linkGenesInfo(1).foldChange;
brainP = linkGenesInfo(1).allVSTissuePval;
brainExp = linkGenesInfo(1).exp;

smallPInd = find(brainP < (0.05/length(brainP)));
length(smallPInd)

bigFCInd = find(brainFC > 2);

[a, b] = ismember(bigFCInd, smallPInd);
finalInd = bigFCInd(a);
sum(a)

h = figure; 
heatmap(brainExp(finalInd, :), [], [],  [], 'TickAngle', 45, ...
        'Colorbar', true, 'colormap', bone)
h = figure
plot(sign(brainFC(kado)).* log(abs(brainFC(kado))), brainP(kado), 'o')
heatmap(log(abs(brainFC(kado))), [], [], [], 'TickAngle', 45, ...
        'Colorbar', true, 'colormap', bone)

% TODO: a table generated from the finalTable in pu
% 2. check them in GTEx and Illumina 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gene to edges count for each group

% TS for illumina links ...

% Illumina links are from the file illArrayMain01.m
load('~/data/general/brainTSlinksRankedInIllData.mat')
% getting the median of each group: 
bloodMed = quantile(rankedCorrMatStr.mat(:, 2:5)', .5);
brainMed = quantile(rankedCorrMatStr.mat(:, 13:15)', .5);
brainDiff = brainMed - bloodMed;
geneIDs = rankedCorrMatStr.geneIDs;

sib = corr(rankedCorrMatStr.mat);

load('~/data/general/commonLinksRankedInIllData.mat')
bloodMed = quantile(rankedCorrMatStr.mat(:, 2:5)', .5);
brainMed = quantile(rankedCorrMatStr.mat(:, 13:15)', .5);
commonDiff = brainMed - bloodMed;

% finding the p-values 
pVals = zeros(size(brainDiff));
for i = 1:length(brainDiff)
    pVals(i) = sum(commonDiff > brainDiff(i)) / length(commonDiff); % 500 is where the p-value of null
                                                        % gets smaller than .05
end

hist(pVals, 40)

% getting the FDR

[a, b] = sort(pVals);
sortedPVals = pVals(b);
sortedBrainDiff = brainDiff(b);
sortedGenes = geneIDs(b, :);
counter = 0;
alpha = .1;
for i = 1:length(pVals)
    if (sortedPVals(i) <= i * alpha / length(pVals))
        counter = counter + 1;
    end
end
counter

% is it related to TSVALUE?
links = sortedGenes(1:counter,:);
linksTemplate = zeros(18494, 18494);
for i = 1:1755
    linksTemplate(links(i, 1), links(i, 2)) = 1;
    linksTemplate(links(i, 2), links(i, 1)) = 1;
end
linksTemplate = sparse(linksTemplate);

sib = finalTable(1).pureEdges .* linksTemplate;
sum(sum(sib))

genes = unique(sortedGenes(1:1755,:));

h = figure
heatmap([sortedBrainMed', sortedBloodMed'], [],[],[], 'Colorbar', true, ...
        'colormap', bone)

h = figure;
subplot(1, 2, 1)
hist(abs(kado(kado <0)), 30)
total = sum(kado < 0)
title(sprintf('%d', total))
subplot(1, 2, 2)
hist((kado(kado >= 0)), 30)
total = sum(kado >= 0)
title(sprintf('%d', total))

sum(sortedBrainMed > 600)
brainCorr = (sortedBrainMed > 600)


h = figure;
heatmap([sortedBrainMed', brainSortedBloodMed'], [],[],[], ...
        'Colorbar', true, 'colormap', bone)

% think of some plot later, for now, just look for the confirmed
% Illumina links

% reproducibiliby for each edge group

% confirmed GTEx links... 

% 3. plot them in a plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% draft testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sum(eqGenes)
h = figure
testExp = this.exp(eqGenes, :);
heatmap(testExp(b(1:10), :), [], [], [], 'Colorbar', true, 'Colormap', ...
        bone)

plot(this.lfc(eqGenes), this.pValues(eqGenes), 'o')
lfcEQ = this.lfc(eqGenes);
indEQ = this.geneInds(eqGenes);
expEQ = this.exp(eqGenes, :);
toBChecked = find(lfcEQ > 2);

testExp = expEQ(toBChecked, :);
heatmap(testExp, [], [], [], 'Colorbar', true, 'Colormap', ...
        bone)


[a, b]= sort(this.pValues(eqGenes));

h = figure
testNotExp = this.exp(~eqGenes, :);
heatmap(testNotExp(1:50, :), [], [], [], 'Colorbar', true, 'Colormap', ...
        bone)

heatmap(this.exp(tsGenes, :), [], [], [], 'Colorbar', true, 'Colormap', ...
        bone)

fcInd = find(this.foldChange >= .8);
pInd = find(this.allVSTissuePval < (pThr02/length(fcInd)));
tsGenes = pInd + fcInd;
tsGenes = tsGenes == 2;
tsGenesInd = this.geneInds(tsGenes);

%genes with fcInd OK
heatmap(this.exp(holu, :), [], [], [], 'Colorbar', true, 'Colormap', ...
        bone)

h = figure
heatmap(this.exp(pInd,:),  [], [], [], 'Colorbar', true, 'Colormap', ...
        bone)


% genes which are NOT in pInd and ARE in FC
[a, b]= ismember(fcInd, pInd);
kado = fcInd(~a);
h = figure
heatmap(this.exp(kado, :),  [], [], [], 'Colorbar', true, 'Colormap', ...
        bone)
% I think they should not be kept. so, we go for pInd, but how
% about those which are in PInd and not in FC? almost half...
[a, b]= ismember(pInd, fcInd);

[a, b] = ismember(tsGenesInd, eqGenesInd);
sum(a)
kado = pInd(~a);
h = figure
heatmap(this.exp(tsGenesInd, :),  [], [], [], 'Colorbar', true, 'Colormap', ...
        bone)

this.foldChange(kado(50:60))

mean(this.exp(kado(50), 1:41))
mean(this.exp(kado(50), 42:53))

h = figure
heatmap(this.exp(tsGenes, :), [],[],[], 'Colorbar', true, ...
        'Colormap', bone)

h = figure
heatmap(this.exp(eqGenes, :), [],[],[], 'Colorbar', true, ...
        'Colormap', bone)

kado = eqGenes + tsGenes;
otGenes = ~kado;

otGenes01 = find(othGenes);
h = figure
heatmap(this.exp(otGenes01(42:79), :), [],[],[], 'Colorbar', true, ...
        'Colormap', bone)

h = figure
heatmap(this.exp(p01(~a),:), [],[],[], 'Colorbar', true, ...
        'Colormap', bone)
this.geneIDs(p01(~a))

% what to do with the links 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defThr = .58
load(sprintf(['~/data/general/' ...
              'TSlinksTable_4QDSCounter_0.8Expr_Ind0.05%.2f.mat'], defThr))

geneSums = zeros(1, 5)
linkSums = [7828 1466 625 294 671]
myBar = zeros(3, 5)
for i = 1:5
    finalTable(i).tissue
    a = finalTable(i).othGenesCount;
    b = finalTable(i).cGenesCount
    c = finalTable(i).tsGenesCount
    genesSums(i) = a + b + c;
    myBar(3, i) = finalTable(i).pureEdgeCount /linkSums(i);
    myBar(2, i) = finalTable(i).tsEdgeCount/linkSums(i);
    myBar(1, i) = finalTable(i).tcEdgeCount/linkSums(i);
end

colormap = 'bone'
h = figure

bar(myBar', 'stack')
P=findobj(gca,'type','patch');
colors = [230, 85, 13;
         99, 99, 99;
         49, 130, 189]/256

for n=1:length(P) 
set(P(n),'facecolor',colors(n, :));
end
set(gca, 'XTick', [1 2 3 4 5])
set(gca, 'XTickLabel', {'brain', 'blood', 'liver', 'lung', 'muscle'})
grid on

line([3.7 4.2], [.65 .65], 'color', colors(1, :), 'LineWidth', 6)
line([3.7 4.2], [.6 .6], 'color', colors(2, :), 'LineWidth', 6)
nline([3.7 4.2], [.55 .55], 'color', colors(3, :), 'LineWidth', 6)
title('Percent of the TS links within each group of genes')

file = '~/data/general/tempFigures/linksShared.pdf'
print(h, '-dpdf', file)

% 10. testing the expression level for the TS links and finding examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
load([ '~/data/general/GPL570GemmaMapNEW.mat'])
defThr = .58
load(sprintf('~/data/general/TSlinksTable_4QDSCounter_0.8Expr_Ind0.05%.2f_TSValueSaved.2f.mat', defThr))
load('~/data/general/tissueExpGenes/allDSExprInfo.mat')

tissue = 'brain'
load(['~/networks/tissues/' tissue '/' ...
      'binaryNet_4QDSCounter_0.8Expr_Ind0.05.mat'])

tissue = 'lung'
load(['~/networks/continuousTSLinksResults/' tissue ...
      'SpaceMat_negAdjusted.mat']) % spaceMat

newNet = spaceMat .* binNet;
holu1 = newNet > .76;
whos holu1
sum(sum(holu1))

book = holu1 + holu2;
sum(sum(book == 2))

i = 1
finalTable(i).tissue
sib = finalTable(i).pureEdges + finalTable(i).tsEdges + ...
      finalTable(i).othEdges + finalTable(i).tcEdges + ...
      finalTable(i).toEdges + finalTable(i).coEdges;
sum(sum(sib > 0))
holu2 = sib > 0; 
whos holu2
sum(sum(holu2))


[a, b, c] = find(finalTable(i).tcEdges);

sym = 'MOG'
[a, ID1] = ismember(sym, gpl570.uniqueSymbols(:))

% finding the half and half links: 
pureSum = zeros(18494, 18494);
tsSum = zeros(18494, 18494);
tcSum = zeros(18494, 18494);
for i = 1:5
    i
    tempNet = (finalTable(i).pureEdges >0) .* i;
    pureSum = pureSum + tempNet;
    
    tempNet = (finalTable(i).tsEdges>0) .* i;
    tsSum = tsSum + tempNet;
    
    tempNet = (finalTable(i).tcEdges>0) .* i;
    tcSum = tcSum + tempNet;
end

sum(sum(tcSum > 0))
sum(sum(tsSum > 0))
sum(sum(pureSum > 0))

% see if a gene has two partners
multiple = zeros(1, 18494);
net = tcSum;
for i = 1:18494
    v = (net(i, :)) +( net(:, i)');
    [a, b] = ismember([1:5], v);
    if(sum(a) >= 2)
        multiple(i) = 1;
    end
end

tcMultiple = multiple;
pureMultiple = multiple;


[a, b] = find(tcMultiple);
ntcSum = sparse(tcSum);

% getting the expression levels for these genes and saving them
myGenes = [4129, 3049, 4970, 711]

%brains %
myGenes = [965 1000]
myGenes = [9601 12053]
myGenes = [5906 4273]

%bloods%
myGenes = [263 3343]

myGenes = [965 1000 9601 12053 5906 4273 263 3343]


sib = allDSExprInfo(1).exprLevel';
holu = mean(sib)
sib = sib - holu;
for i = 2:21
    kado = allDSExprInfo(i).exprLevel;
    holu = mean(kado)
    kado = kado - holu;
    sib = [sib kado'];
end

brainMat = sib(myGenes, 1:12);
bloodMat = sib(myGenes, 13:21);

h = figure
heatmap(brainMat, [], [], [], 'Colorbar', true, 'Colormap', ...
        bone)

print(h, '-dpdf', '~/data/general/tempFigures/randomThing.pdf')

h = figure
heatmap(brainMat, [], [], [], 'Colorbar', true, 'Colormap', ...
        bone)

print(h, '-dpdf', '~/data/general/tempFigures/brainMatPure.pdf')
h = figure
heatmap(bloodMat, [], [], [], 'Colorbar', true, 'Colormap', ...
        pink)
title([gpl570.uniqueSymbols(myGenes)])
print(h, '-dpdf', '~/data/general/tempFigures/bloodMatPureGeneNames.pdf')


gpl570.uniqueSymbols(478)
finalTable(i).tsEdges(478, :)
sib(478, :)
gpl570.uniqueSymbols(8798)


% 11. comparing the p-values vs the distance. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this whole part is moved to theContinuousTSlinkSelection.m part
% 11. 

% 11.5. getting the TS links - old version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}

R = 0;

for t = 1:5
t
    tissue = tissues{t};
    % load(['~/networks/tissues/' tissue '/' ...
    %       'binaryNet_0.5DSCountCorr_0.8Expr.mat'])
    %tempNet/tissue Networks

    % load(['~/networks/tissues/' tissue '/' ...
    %       'binaryNet_4QDSCounter_0.8Expr_Ind0.05.mat'])


    % this is the final binary network with FDR
    expThr = '0.8'
    load( ['~/networks/tissues/' tissue '/' ...
           'binaryNet_FDR5e-5_' expThr 'Expr_Ind0.10.mat'])
    sum(sum(binNet))
    
    % this is the semi-binary network 
    %    load(['~/networks/tissues/' tissue '/topSumNet_5_0.8Expr_Ind0.10.mat'])
    % max(max(binNet))
    % modBinNet = binNet;
    % binNet = binNet == 5;
    % sum(sum(binNet))

    load(['~/networks/continuousTSLinksResults/' tissue ...
          'SpaceMat_negAdjusted.mat']) % spaceMat
    spaceMat = spaceMat + 0.001;
    negAdjNet =  spaceMat .* binNet;
    [naa, nab, na] = find(negAdjNet);
    naSelected = na > .63;
    sum(naSelected)
    clear spaceMat;

    load(['~/networks/continuousTSLinksResults/' tissue ...
          'SpaceMat_pVals.mat']) % spaceMat
    spaceMat = spaceMat + 1e-10;
    pVals =  spaceMat .* binNet; 
    [pva, pvb, pv] = find(pVals);
    pv = pv - 10e-12;

    % load(['~/networks/continuousTSLinksResults/' ...
    %       'random10_smallMat_negAdjusted.mat'])
    % smallMat = smallMat / (10*43*500);
    % max(max(smallMat))
    % [sma, smb, smna] = find(smallMat);
    % clear sma, smb;
    % sum(smna > .6)/length(smna)

    load('~/networks/continuousTSLinksResults/random2_smallMat_pVals.mat')
    [rpva, rpvb, rpv] = find(smallMat);
    clear rpva, rpvb;
    sum(rpv < 5e-4)/length(rpv)

    pvSelected = pv <= 5e-4;
    sum(pvSelected)
    pvSelected = pvSelected .*2;

    sib = pvSelected + naSelected;

    R = R + sum(sib == 3)
    
    % testing 
    % testa = naa(sib == 1);
    % testb = nab(sib == 1);
    % testna = na(sib == 1);
    % k = 6
    % [testa(k), testb(k), testna(k)]
    
    % [sorteda, sortedb] = sort(testna, 'descend');
    % stesta = testa(sortedb);
    % stestb = testb(sortedb);
    
    % k = 100
    % [stesta(k), stestb(k)]
    % sorteda(k)
    % testing
    
    % another test
    % kado = find(sib == 3);
    % testa = naa(kado);
    % testb = nab(kado);
    % testna = na(kado);
    
    % [sortedna sinds] = sort(testna, 'descend');
    % stestb = testb(sinds);
    % stesta = testa(sinds);
    
    % k = length(testa)
    % [stestb(k) stesta(k)]
    % sortedna(k)
    % another test
    
    kado = find(sib == 3);
    TSlinks{t}.a = naa(kado);
    TSlinks{t}.b = nab(kado);
    TSlinks{t}.pv = pv(kado);
    TSlinks{t}.na = na(kado);
    TSlinks{t}.tissue = tissue;

    w = zeros(1, length(kado));
    for i = 1:length(TSlinks{t}.a)
        w(i)= modBinNet(TSlinks{t}.a(i), TSlinks{t}.b(i));
    end
    TSlinks{t}.linkWeight = w;
end

% 11.75. getting the TS links - new
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}

R = 0;
fdrthr = .001

for t = 1:5
t
    tissue = tissues{t};
 
    % this is the final binary network with FDR
    expThr = '0.8'
    load( ['~/networks/tissues/' tissue '/' ...
           'binaryNet_FDR5e-5_' expThr 'Expr_Ind0.10.mat'])
    sum(sum(binNet))
    
    % this is the semi-binary network 
    %    load(['~/networks/tissues/' tissue '/topSumNet_5_0.8Expr_Ind0.10.mat'])
    % max(max(binNet))
    % modBinNet = binNet;
    % binNet = binNet == 5;
    % sum(sum(binNet))

    load(['~/networks/continuousTSLinksResults/qnorms/' tissue  ...
          '/SpaceMat_NA_qnorms.mat']) % spaceMat
    spaceMat = spaceMat + 0.0001;
    negAdjNet =  spaceMat .* binNet;
    [naa, nab, na] = find(negAdjNet);
        % load the random tables for NA
    % from .40 to 1.00
        
%             [distr(t).fi distr(t).xi] = ksdensity(na,'Support', [0, ...
%                                 1]);
% end
% save('~/resultsAndFigures/TSlinks/WholeTissueNetNADists.mat', 'distr')

    load('~/resultsAndFigures/TSlinks/NegAdjusted_ExpAdjusted_FDR.mat')
    myNAFDR = tFDR(t, :);
    ind = min(find(myNAFDR <= fdrthr))
    myNA = .4 + ind/100;
    naSelected = na >= myNA;
        sum(naSelected)

    % top = quantile(na(naSelected), .66)
    naSelected = naSelected .* 2;
    % naSelected(na >= top) = 3;
    % sum(naSelected > 0)

    load(['~/networks/continuousTSLinksResults/negAdjusted/' tissue '/' tissue ...
          'SpaceMat_pVals_oneSided.mat']) % spaceMat
    pVals =  spaceMat .* binNet; 
    [pva, pvb, pv] = find(pVals);
    % load the random tables for pv
    thePs = [1e-6, 5e-6, 1e-5, 5e-5, 1e-4, ...
                        5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, ...
                        5e-1];
    load('~/resultsAndFigures/TSlinks/pVals12Counts_FDR.mat')
    myFDR = tFDR(t, :);
    ind = max(find(myFDR <= fdrthr))
    myP = thePs(ind)
    pvSelected = pv <= myP;
    % top = quantile(pv(pvSelected), .33)
    % pvSelected = pvSelected + 0;
    % pvSelected(pv <= top) = 3;
    sum(pvSelected >0)
    
    sib = pvSelected + naSelected;

    R = sum(sib > 0)
    
    % testing 
    % testa = naa(sib == 1);
    % testb = nab(sib == 1);
    % testna = na(sib == 1);
    % k = 6
    % [testa(k), testb(k), testna(k)]
    
    % [sorteda, sortedb] = sort(testna, 'descend');
    % stesta = testa(sortedb);
    % stestb = testb(sortedb);
    
    % k = 100
    % [stesta(k), stestb(k)]
    % sorteda(k)
    % testing
    
    % another test
    % kado = find(sib == 3);
    % testa = naa(kado);
    % testb = nab(kado);
    % testna = na(kado);
    
    % [sortedna sinds] = sort(testna, 'descend');
    % stestb = testb(sinds);
    % stesta = testa(sinds);
    
    % k = length(testa)
    % [stestb(k) stesta(k)]
    % sortedna(k)
    % another test
    
    kado = find(sib > 0);
    sib = sib(kado);
    TSlinks{t}.a = naa(kado);
    TSlinks{t}.b = nab(kado);
    TSlinks{t}.pv = pv(kado);
    TSlinks{t}.na = na(kado);
    TSlinks{t}.sib = sib;
    TSlinks{t}.pFDR = myP;
    TSlinks{t}.naFDR = myNA;
    TSlinks{t}.tissue = tissue;

    % w = zeros(1, length(kado));
    % for i = 1:length(TSlinks{t}.a)
    %     w(i)= modBinNet(TSlinks{t}.a(i), TSlinks{t}.b(i));
    % end
    % TSlinks{t}.linkWeight = w;
end


save('~/resultsAndFigures/TSlinks/TSlinks_FDR001_NA_EA.mat', ...
     'TSlinks')

save('~/resultsAndFigures/TSlinks/TSlinks_FDR001_NA_EA.mat', ...
     'TSlinks')

kado = 1;
i = 2;
ts = filTS{i};
tissues{i}
book = find(ts.pv > .001);
book = find(ts.na < .15);
j = 1
kado = book(j)
[ts.a(kado) ts.b(kado)]
[ts.pv(kado) ts.na(kado)]
j = j + 1

[ts.a(kado) ts.b(kado) ts.pv(kado) ...
 ts.na(kado)]

[TSlinks{i}.a(kado) TSlinks{i}.b(kado) TSlinks{i}.pv(kado) ...
 TSlinks{i}.na(kado)]
kado = kado + 1

save('~/resultsAndFigures/TSlinks/TSlinks_FDR0012.mat', 'TSlinks')
save('~/resultsAndFigures/TSlinks/TSlinks_FDR001.mat', 'TSlinks')
save('~/resultsAndFigures/TSlinks/TSlinks_AllNetworkTSValues.mat', 'TSlinks')
save('~/resultsAndFigures/TSlinks/TSlinks_FDR0012Union33Inter66.mat', ...
     'TSlinks')
save('~/resultsAndFigures/TSlinks/TSlinks_FDR0012Inter.mat', 'TSlinks')
save('~/resultsAndFigures/TSlinks/TSlinksExtended.mat', 'TSlinks')
save('~/resultsAndFigures/TSlinks/TSlinks.mat', 'TSlinks')
save('~/resultsAndFigures/TSlinks/tissues.mat', 'tissues')
load('~/resultsAndFigures/TSlinks/TSlinksExtended.mat')
load('~/resultsAndFigures/TSlinks/TSlinks_FDR0012.mat')
filTS = TSlinks
load('~/resultsAndFigures/TSlinks/TSlinks_AllNetworkTSValues.mat')
allTS = TSlinks;
load(['~/resultsAndFigures/TSlinks/' ...
      'TSlinks_FDR0012Union33Inter66.mat'])

% getting sum of all the TS links:
s = 0;
for i = 1:5
    TSlinks{i}.
end

% plotting the values
clear
folder = '~/resultsAndFigures/TSlinks/figures/'
load('~/resultsAndFigures/TSlinks/TSlinks_AllNetworkTSValues.mat')
for t = 1:5
    tissue = TSlinks{t}.tissue
    h = figure
    plot(log10(TSlinks{t}.pv),TSlinks{t}.na, '.', 'color', 'k')
    xlabel('log10 p-value')
    ylabel('distance measure')
    hold on
    line([log10(min(TSlinks{t}.pv)) 0], [TSlinks{t}.naFDR, TSlinks{t}.naFDR], ...
         'LineWidth', 2,'Color', 'r')
    line([log10(TSlinks{t}.pFDR), log10(TSlinks{t}.pFDR)], [0 , 1], ...
         'LineWidth', 2,'Color', 'r')
    title(sprintf('Correlation of the two methods for tissue %s', ...
                  TSlinks{t}.tissue))
    r = corr(TSlinks{t}.pv, TSlinks{t}.na, 'type','Spearman');
    text(-2.5, .96, sprintf('corr: %f', r))
    text(-2.5, .92, sprintf('total edges: %d', length(TSlinks{t}.na)))
    first = TSlinks{t}.na > TSlinks{t}.naFDR;
    text(-2.5, .88, sprintf('accepted by distance: %d', sum(first)))
    second = TSlinks{t}.pv < TSlinks{t}.pFDR;
    text(-2.5, .84, sprintf('accepted by p-value: %d', sum(second)))
    text(-2.5, .8, sprintf('intersect: %d', sum((first + second) == ...
                                                2)))
    
    file = sprintf('twoMethodsCorr_%s', tissue)
    print(h, '-dpdf', [folder file '.pdf'])
        print(h, '-depsc', [folder file '.eps'])
end

book1 = find(selectTSlinks{1}.pv>1e-4);
book2 = find(selectTSlinks{1}.pv<.4);
[a, b] = ismember(book2, book1);
a = ~a;
book = book2(a);
book = book1;
i = 1
selectTSlinks{1}.a(book(i))
selectTSlinks{1}.b(book(i))
selectTSlinks{1}.na(i)
selectTSlinks{1}.pv(i)
i = i + 1

load('~/resultsAndFigures/TSlinks/TSlinks_AllNetworkTSValues.mat')
allTSlinks = TSlinks;
h = figure
t = 1
plot(log10(1-allTSlinks{t}.pv), allTSlinks{t}.na, '.')

book = find(allTSlinks{t}.na<.15)
allTSlinks{t}.a(book)
allTSlinks{t}.b(book)
allTSlinks{t}.na(book)
allTSlinks{t}.pv(book)

% get naFDR which is not specific for any other tissue.
clear
load('~/resultsAndFigures/TSlinks/TSlinks_FDR0012.mat')

totalNet = zeros(18494, 18494);
for t = 1:5
    myTS = TSlinks{t}
    kado = myTS.sib > 1;
    myNet = sparse(myTS.a(kado), myTS.b(kado), ones(sum(kado), 1), ...
                   18494, 18494);
    totalNet = totalNet + myNet;
end

% 80 links ovrelap in other tissues. Find the maximum of the link
% and which tissue does it belong to, for each tissue. This is from
% theContinuous....m , lines ~1330 

clear
load('~/resultsAndFigures/TSlinks/TSlinks_FDR0012.mat')

load('~/resultsAndFigures/TSlinks/TSlinks_FDR001_NA_EA.mat')
TSinOthers = zeros(5, 5);
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'};
for i = 1:5
    i
    tissue = TSlinks{i}.tissue
    tslCount = sum((TSlinks{i}.sib >  0));
    
    tempNet = sparse(TSlinks{i}.a, TSlinks{i}.b, 1, 18494, 18494);

    otherInd = [1:5];
    otherInd(i) = [];
    maxTS = zeros(tslCount, 1);
    nextTissue = zeros(tslCount, 1);
    
    minTSp = ones(tslCount, 1);
    nextPTissue = zeros(tslCount, 1);
    for j = 1:4
        ctissue = tissues{otherInd(j)}
        
        
        load(['~/networks/continuousTSLinksResults/negAdjusted_ExpAdjusted/' ctissue  ...
          '/SpaceMat_NA_EA.mat']) % spaceMat
        % load(['~/networks/continuousTSLinksResults/' ctissue '/' ctissue ...
        %       'SpaceMat_negAdjusted.mat']) % spaceMat
        
        spaceMat = spaceMat + .001;
        book = (spaceMat) .* tempNet;

        [a, b, c] = find(book);
        maxTS = max(c, maxTS);
        book = maxTS - c;
        nextTissue(book == 0) = otherInd(j);

        load(['~/networks/continuousTSLinksResults/negAdjusted/' ctissue '/' ctissue ...
              'SpaceMat_pVals_oneSided.mat']) % spaceMat
        book1 = spaceMat .* tempNet;
        book2 = tempNet;
        book2(book1 > 0) = book1(book1 >0);
        
        [a, b, c] = find(book2);
        minTSp = min(c, minTSp);
        book = minTSp - c;
        nextPTissue(book == 0) = otherInd(j);
    end
    TSlinks{i}.minTSp = minTSp;
    TSlinks{i}.nextPTissue = nextPTissue;
    TSlinks{i}.maxTSna = maxTS;
    TSlinks{i}.nextNATissue = nextTissue;
end

save(['~/resultsAndFigures/TSlinks/' ...
      'TSlinks_FDR001_NA_EA_withNextTissueValues.mat'], 'TSlinks')

save(['~/resultsAndFigures/TSlinks/' ...
      'TSlinks_FDR0012_withNextTissueValues.mat'], 'TSlinks')

% looking at the next high value links . 
sib = corr(log10(TSlinks{1}.minTSp), TSlinks{1}.maxTSna);
k = TSlinks{1}.na > TSlinks{1}.naFDR;
plot((TSlinks{1}.na(k)), TSlinks{1}.maxTSna(k), '.') 
myTSna = TSlinks{1}.maxTSna(k); 
myna = TSlinks{1}.na(k); 
hist(myTSna,30)
sum(myTSna > .4)
h = figure
hist(myna, 30)

book = find(TSlinks{1}.maxTSna > .5);
i = book(5)

[TSlinks{1}.a(i), TSlinks{1}.b(i)]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getting the links and networks and etc...
clear 
load('~/resultsAndFigures/TSlinks/TSlinks.mat')
links = TSlinks{2}

% load the brain network and study the network. 
tissue = 'brain'
expThr = '0.8'
load( ['~/networks/tissues/' tissue '/' ...
       'binaryNet_FDR5e-5_' expThr 'Expr_Ind0.10.mat'])

load('~/data/general/GPL570GemmaMapNEW.mat')
sumTot = sum(binNet) + sum(binNet');
inBrain = sumTot > 0;
hist(sumTot(inBrain), 30)

sib = unique([links.a, links.b]);
links.a = [links.a', 18494];
links.b = [links.b', 18494];
links.na = [links.na', 0];
holu = sparse(links.a, links.b, ones(1, length(links.b)));
holu(18494, 18494) = 0;

% node distribution of the TS network
s1 = sum(holu);
s2 = sum(holu');
sumTS = s1 + s2;

sum(sumTS == 4)

inBrainTS = (sumTS > 0);
hold on
hist(sumTS(inBrainTS), 30, 'Color', 'r')

[ma, mb] = max(sumTS)

maxNeighbours = find(holu(mb ,:) + holu(:, mb)');

% 12. New pure link study
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

% load the list of tissues
load('~/resultsAndFigures/TSlinks/tissues.mat')

% 1. Identify the equally expressed genes. 

% 1.1 reading the genes
myGenes = zeros(1, 18494);
for t = 1:5
    tissue = tissues{t};
    load(['~/data/general/tissueExpGenes/' tissue ...
          'ExpGenes0.8.mat'])
    myGenes = myGenes + expGenesInd;
end
expGenes = myGenes > 0;
sum(expGenes)

% 1.2 get their exp levels
%clear
load('~/data/general/linkExprInfo/totalMeanExpCentered.mat', ...
     'totalMeanExpCentered')
load('~/data/general/linkExprInfo/dataSetExpInf.mat')

tExpMat = vec2mat([dataSetExpInf(:).meanExp], 18494)';
expMat = zscore(tExpMat);

tic
p = zeros(1, 18494);
for j = 1:length(p)
    p(j) = anova1(expMat(j,:), {dataSetExpInf(:).tissue}, 'off');
end
toc

save(['~/resultsAndFigures/expAndCorrRegression/' ...
      'TSExpressionAnova.mat'], 'p')
load('~/resultsAndFigures/expAndCorrRegression/TSExpressionAnova.mat')

% 12.5: trying to get TS links using regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ps = p > (.1 / 1849);
sum(ps)
commonGenes = ps; 
book = (ps + expGenes) == 2;
sum(book)
sib = find(book);
sib(1:10)
p(sib(1:10))

myExp = expMat(book, :);
Y = pdist(myExp);
Z = linkage(Y, 'complete');
clusterCount = 6
T = cluster(Z, 'maxclust', clusterCount);

% sorting based on the cluster
[a, b] = sort(T);
sortedExpr = myExp(b, :);
sortedIDs = uniqueGenesIDs(b);
sortedInd = uniqueInd(b);

% figure
heatmap(sortedExpr, [],  [],  [], 'Colorbar', true, 'Colormap', ...
         'summer')

% >>>>>>>>>>>>>>>>>>>>>>>>> Just trying the reg
featureMat = zeros(53, 5);
featureMat(1:9, 1) = 1;
featureMat(10:24, 2) = 1;
featureMat(25:31, 3) = 1;
featureMat(32:41, 4) = 1;
featureMat(42:53, 5) = 1;

coefs = cell(1, 18494);
RSred = zeros(1, 18494);
for i = 1:18494
    shuffInd = datasample([1:53], 53, 'Replace', false);
    YS = expMat(i, shuffInd);
    lm = LinearModel.fit(featureMat, YS, 'Intercept', false);
    %    lm = LinearModel.fit(featureMat, expMat(i, :), 'Intercept', false);
    %    anova1(expMat(i, :), {dataSetExpInf(:).tissue})
    coefs{i} = lm.Coefficients;
    RSred(i) = lm.Rsquared.Ordinary;
    % anova1(expMat(2,:), {dataSetExpInf(:).tissue})
    % RSEd(i) = lm2.Rsquared.Ordinary;
    % SSE(i) = lm2.SSE;
    % SSR(i) = lm2.SSR;
    % SST(i) = lm2.SST;
    % x1int1 = lm2.Coefficients(2,1);
    % x1int2 = dataset2struct(x1int1);
    % coef1(i) = x1int2.Estimate;
    % x2int1 = lm2.Coefficients(3,1);
    % x2int2 = dataset2struct(x2int1);
    % coef2(i) = x2int2.Estimate;
end

res.coefs = coefs;
res.RSred = RSred;
res.shCoefs = coefs;
res.shRSred = RSred;

save(['~/resultsAndFigures/expAndCorrRegression/' ...
      'TSExpressionRegression.mat'], 'res')
load(['~/resultsAndFigures/TSGenes/' ...
      'TSExpressionRegression.mat'])

% meanExp = mean(expMat');
% noDiffInd = res.RSred < .25;
% expInd = meanExp > -.5;
% selected = (noDiffInd + expInd) == 2;
% inds = find(selected);
% length(inds)

% getting each of the pvalues
cfMat = zeros(18494, 5);
pValMat = zeros(18494, 5);
for i = 1:18494
    x = res.coefs{i};
    x1 = dataset2struct(x);
    cfMat(i, :) = [x1.Estimate];
    pValMat(i, :) = [x1.pValue];
end

% the average of each tissue
bloodMean = mean(expMat(:, 1:9)');
lungMean = mean(expMat(:, 10:24)');
muscleMean = mean(expMat(:, 25:31)');
liverMean = mean(expMat(:, 32:41)');
brainMean = mean(expMat(:, 42:53)');

tissueMeanExp = [bloodMean', lungMean', muscleMean', liverMean', ...
                 brainMean'];

[maxTa, maxTb] = max(tissueMeanExp');

[pa, pb] = min(pValMat');
% because of the unbalanced design, I get more of the lung tissues
% here, so I am not sure how to trust this. 
% tsInd = a < 1e-7;
% sum(tsInd)

[cfa, cfb] = max(cfMat'); 

cfMat = tissueMeanExp;
ccfmat= tissueMeanExp;
for i = 1:18494
    ccfmat(i, cfb(i)) = -3;
end
[cfa2, cfb2] = max(ccfmat');

fc = cfa - cfa2;
kado = find(fc >= 1);
length(kado)

fcPval = zeros(1, length(kado));
for i = 1:length(kado)
    fcPval(i) = pValMat(kado(i), cfb(kado(i)));
end

selectedTS = kado(fcPval < .1/18494);


% now getting the tissue for the TS genes: 
% looking for relatively high mean and low RSQRD
% sort the mean coeffs with the low RSQRD. 

anova1(expMat(5,:), {dataSetExpInf(:).tissue})

% <<<<<<<<<<<<<<<<<<<<<<<<< 

pVals = p;
[a, b] = sort(pVals);
sortedPVals = pVals(b);
counter = 0;
alpha = .001;
pos = zeros(1, length(pVals));
for i = 1:length(pVals)
    if (sortedPVals(i) <= i * alpha / length(pVals))
        counter = counter + 1;
    end
end
pos(b(1:counter)) = 1;

pos = ~logical(pos);

% adjP = p * length(p);
% holu = adjP > 5e-2;

final = pos + expGenes;
eqGenes = final == 2;
sum(eqGenes) % I have about 7k genes in total which are considered as
           % commonly expressed between the tissues. 

% cluster and sort - just checking 
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
uniqueGenesExpr = totalMeanExpCentered(final, :);
Y = pdist(uniqueGenesExpr);
Z = linkage(Y, 'complete');
clusterCount = 10
T = cluster(Z, 'maxclust', clusterCount);

% sorting based on the cluster
[a, b] = sort(T);
sortedExpr = uniqueGenesExpr(b, :);
h = figure; 
heatmap(sortedExpr, [], [], [], 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'bone')
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 

% >>>>>>>>>>>>>>>>>>>>>>>>> This is only for exploratory plotting. 
selectedTSTissue = cfb(selectedTS);
length(selectedTSTissue)
[sta, stb] = sort(selectedTSTissue);

semExp = expMat(selectedTS(stb), :);

Y = pdist(semExp);
Z = linkage(Y, 'complete');
clusterCount = 6
T = cluster(Z, 'maxclust', clusterCount);

[tsa, tsb] = sort(T);
sortedExpr = semExp(tsb, :);

heatmap(semExp, [],  [],  [], 'Colorbar', true, 'Colormap', ...
         'summer')

% <<<<<<<<<<<<<<<<<<<<<<<<<<< end of exploratory plotting. 
% This part is obsolete, I have a new one.
% 13_obsolete. Getting the TS genes using ANOVA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< not regression no
% more. Back to anova for tissues: 
clear
load('~/data/general/linkExprInfo/dataSetExpInf.mat')
%load('~/resultsAndFigures/TSlinks/TSlinks.mat') % 'TSlinks'
load('~/resultsAndFigures/TSlinks/TSlinks_FDR001.mat') % 'TSlinks'

% get the tissues
tissues = unique({dataSetExpInf(:).tissue})

% get the expression mat
expMat = vec2mat([dataSetExpInf(:).meanExp], 18494)';
fcThr = 1;
pValThr = .05/(18494);

% get the tissue Index for each tissue. 
for t = 1:length(tissues)
    tissue = tissues{t}
    tissueCell = {dataSetExpInf(:).tissue};
    dsIndCell = strfind(tissueCell, tissue);
    allDSInd = find(cellfun(@isempty, dsIndCell));
    labels = cell(length(tissueCell),1);
    labels(:) = {'that'};
    labels(allDSInd) = {'this'};
    indexing = ones(size(dsIndCell));
    allDSInd
    indexing(allDSInd) = 0;
    indexing = logical(indexing);
    indexing
    tp = zeros(1, 18494);
    fChange = zeros(1, 18494);
    
    for j = 1: 18494
        exp = expMat(j, :);
        tp(j) = anova1(exp, labels , 'off' );
        fChange(j) = mean(exp(indexing)) - mean(exp(~indexing));
        fcMed(j) = median(exp(indexing)) - median(exp(~indexing));
    end
    
    tsGenes(t).anp = tp;
    tsGenes(t).fc = fChange;
    tsGenes(t).fcMed = fcMed;
    tsGenes(t).tissue = tissue;
end

t = 5
k = 4
tsGenes(t).tissue
[tsGenes(t).anp(k), tsGenes(t).fc(k)]

save('~/resultsAndFigures/TSGenes/tsGenes_TSlinksFDR001.mat', 'tsGenes')
load('~/resultsAndFigures/TSGenes/tsGenes.mat')

% 13_new. I use the anova for all 5 of the tissues. 
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< not regression no
% more. Back to anova for tissues: 
clear
load('~/data/general/linkExprInfo/dataSetExpInf.mat')
%load('~/resultsAndFigures/TSlinks/TSlinks.mat') % 'TSlinks'
load('~/resultsAndFigures/TSlinks/TSlinks_FDR0012.mat') % 'TSlinks'

% get the tissues
tissues = unique({dataSetExpInf(:).tissue})

% get the expression mat
expMat = vec2mat([dataSetExpInf(:).meanExp], 18494)';

%expMat = zscore(expMat);

bloodMean = mean(expMat(:, 1:9)');
lungMean = mean(expMat(:, 10:24)');
muscleMean = mean(expMat(:, 25:31)');
liverMean = mean(expMat(:, 32:41)');
brainMean = mean(expMat(:, 42:53)');

tissueMeanExp = [bloodMean', lungMean', muscleMean', liverMean', ...
                 brainMean'];

[maxTa, maxTb] = max(tissueMeanExp');

% fcThr = 1;
% pValThr = .05/(18494);

maxSecTb = zeros(1, 18494);
p = zeros(1, 18494);
topTwoFC = zeros(1, 18494);
% for each gene:
for i = 1:18494
    % 
    labels = {dataSetExpInf(:).tissue};
    
    % get the anova p value between the five tissues. 
    p(i) = anova1(expMat(i,:), labels, 'off');
    
    % get the FC between the maximum top two. 
    tempExp = tissueMeanExp(i, :);
    tempExp(maxTb(i)) = 0;
    [seca, secb] = max(tempExp);
    maxSecTb(i) = secb;
    topTwoFC(i) = maxTa(i) - seca;
    % record the p-value, record the tissue of the top exp. 
end

tsGenes.tissues = {'blood', 'lung', 'skeletalMuscle', 'liver', 'brain'};
tsGenes.p = p;
tsGenes.maxTissue = maxTb;
tsGenes.maxSecTissue = maxSecTb;
tsGenes.topTwoFC = topTwoFC;

k = 13174
tsGenes.p(k)
tsGenes.maxTissue(k)
tsGenes.maxSecTissue(k)
tsGenes.topTwoFC(k)

fcThr = 3;
pValThr = .1/(18494);

sib1 = (topTwoFC > fcThr);
sib2 = (p < pValThr);

TSGenes = find((sib1 + sib2) == 2);

tsGenes3fc = TSGenes;
tsGenes2fc = TSGenes;

sum(maxTb(TSGenes) == 5)

save('~/resultsAndFigures/TSGenes/tsGenes_raw.mat', 'tsGenes')
load('~/resultsAndFigures/TSGenes/tsGenes_raw.mat')

% I am checking this with the regression results >>> 
% fcMat = vec2mat([tsGenes(:).fc], 18494);
% fcMatBin = fcMat > 1;
% g1s = sum(fcMatBin);
% g1s = g1s > 0;
% anpMat = vec2mat([tsGenes(:).anp], 18494);
% anpMatBin = anpMat <= pValThr;
% g2s = sum(anpMatBin);
% g2s = g2s > 0;

% sib = (g1s + g2s);
% sum(sib == 2)
% sib = sib == 2;
% sum(sib)

% binInd = zeros(1, 18494);
% binInd(selectedTS) = 1;
% sum(binInd)
% sum(sib)
% final = binInd + sib;

% sum(final == 2)

% finalTSGenesInd = find(final == 2);

% [a, b] = max(fcMat);
% tissueInd = b(finalTSGenesInd);
% hist(tissueInd)

clear
load('~/resultsAndFigures/TSGenes/tsGenes_raw.mat')

fcThr = 2;
pValThr = .1/(18494);

sib1 = (tsGenes.topTwoFC > fcThr);
sib2 = (tsGenes.p < pValThr);

tsg = find((sib1 + sib2) == 2);

TSGenes.tissueInd = tsGenes.maxTissue(tsg);
TSGenes.ind = tsg;
TSGenes.tissues = tsGenes.tissues;;

%save('~/resultsAndFigures/TSGenes/finalTSGenes_FDR.1_test.mat', 'TSGenes');
save('~/resultsAndFigures/TSGenes/finalTSGenes_FDR.1_FC2log.mat', ...
     'TSGenes');
save('~/resultsAndFigures/TSGenes/finalTSGenes_FDR.1_FC1.5.mat', 'TSGenes');
save('~/resultsAndFigures/TSGenes/finalTSGenes_FDR.1_FC1.mat', ...
     'TSGenes');
load('~/resultsAndFigures/TSGenes/finalTSGenes_FDR.1_FC2.mat');
load('~/resultsAndFigures/TSGenes/finalTSGenes_FDR.1_FC1.5.mat');
%<<<<<<< 

figure
i = 2
hist(tsGenes(i).fcMed, 40)
title(tissues{i})
% save the results
save('~/resultsAndFigures/expAndCorrRegression/tsGenes_fcMed_pValue_fc.mat', 'tsGenes')

% >>>>>>>>>>>>>>>>>>>>>>>>> This is only for exploratory plotting. 
selectedTSTissue = cfb(selectedTS);
length(selectedTSTissue)
[sta, stb] = sort(selectedTSTissue);

semExp = expMat(selectedTS(stb), :);

Y = pdist(semExp);
Z = linkage(Y, 'complete');
clusterCount = 6
T = cluster(Z, 'maxclust', clusterCount);

[tsa, tsb] = sort(T);
sortedExpr = semExp(tsb, :);

heatmap(semExp, [],  [],  [], 'Colorbar', true, 'Colormap', ...
         'summer')

% 15. Getting the final table
% <<<<<<<<<<<<<<<<<<<<<<<<<<< end of exploratory plotting. 
% Get the TS links, half and half and everything
clear

tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}
% load  the links
load('~/resultsAndFigures/TSlinks/tissues.mat')
%load('~/resultsAndFigures/TSlinks/TSlinksExtended.mat')

FDR = '0012'
FC = '3'

load(sprintf(['~/resultsAndFigures/TSlinks/' ...
              'TSlinks_FDR%s_withNextTissueValues.mat'], FDR))


FDR = '001'
FC = '3'

load(sprintf(['~/resultsAndFigures/TSlinks/' ...
              'TSlinks_FDR%s_NA_EA_withNextTissueValues.mat'], FDR))


% load the TS and common genes. 
load(sprintf('~/resultsAndFigures/TSGenes/finalTSGenes_FDR.1_FC%slog.mat', ...
             FC));
% load(['~/resultsAndFigures/TSGenes/' ...
%       'TSExpressionAnova.mat'])
% cGenes = p > .1/1849;

% another way of getting the everywhere expressed genes, using the
% genes expressed everywhere and leave the specificity alone.
load('~/data/general/linkExprInfo/dataSetExpInf.mat')
mat12 = zeros(18494, 53);
mat23 = zeros(18494, 53);
mat13 = zeros(18494, 53);
for i = 1:53
    q = quantile(dataSetExpInf(i).meanExp, 5);
    mat12(:, i) = dataSetExpInf(i).meanExp > q(3);
    mat23(:, i) = dataSetExpInf(i).meanExp > q(4);
    mat13(:, i) = dataSetExpInf(i).meanExp > q(2);
end

% sum over the tissues: 
t12 = zeros(18494, 5);
t12(:,1) = sum(mat12(:,1:9)');
t12(:,2) = sum(mat12(:,10:24)');
t12(:,3) = sum(mat12(:,25:31)');
t12(:,4) = sum(mat12(:,32:41)');
t12(:,5) = sum(mat12(:,42:53)');
sib = max(t12);
thr = sib * .8

t23 = zeros(18494, 5);
t23(:,1) = sum(mat23(:,1:9)');
t23(:,2) = sum(mat23(:,10:24)');
t23(:,3) = sum(mat23(:,25:31)');
t23(:,4) = sum(mat23(:,32:41)');
t23(:,5) = sum(mat23(:,42:53)');

t13 = zeros(18494, 5);
t13(:,1) = sum(mat13(:,1:9)');
t13(:,2) = sum(mat13(:,10:24)');
t13(:,3) = sum(mat13(:,25:31)');
t13(:,4) = sum(mat13(:,32:41)');
t13(:,5) = sum(mat13(:,42:53)');

expMat12 = zeros(18494, 5);
expMat23 = zeros(18494, 5);
expMat13 = zeros(18494, 5);
for i = 1:5
    expMat12(:, i) = t12(:, i) >= thr(i);
    expMat23(:, i) = t23(:, i) >= thr(i);
    expMat13(:, i) = t13(:, i) >= thr(i);
end
cGenes = sum(expMat12') == 5;
sum(cGenes)
tempInds = find(cGenes);
cGenes = sum(expMat23') == 5;
sum(cGenes)
% tempInds = find(cGenes);
cGenes = sum(expMat13') == 5;
sum(cGenes)
tempInds = find(cGenes);
save('~/resultsAndFigures/tissueNetworkStudies/the13AllTissueGenes.mat', ...
     'tempInds')

%%%% IMPORTANT: modify the file name based on which cgenes you used 
% for each tissue: 
% 
rsqs = zeros(5, 1);
for t = 1:length(tissues)
    tissue = tissues{t}
    
    % get the tissue
    [a, tempTissueID] =  ismember(tissue, TSGenes.tissues(:));
    
    % get the ts genes
   ttsg = TSGenes.tissueInd == tempTissueID;
   ttsgInd = TSGenes.ind(ttsg);
   clear tempTissueID
   % get the common genes
   cgInd = find(cGenes);
   
   % get the links:
   na1 = TSlinks{t}.na >= TSlinks{t}.naFDR;
   na2 = (TSlinks{t}.maxTSna <= .4);
   na = (na1 + na2) == 2;

   myTSlinks.a = TSlinks{t}.a(na);
   myTSlinks.b = TSlinks{t}.b(na);
   myTSlinks.tsValue = TSlinks{t}.na(na);
   totalGeneCount = length(unique([myTSlinks.a, myTSlinks.b]));
   links = zeros(length(myTSlinks.b), 2);

   % mark the TS genes
   [ats, bts] = ismember(myTSlinks.a, ttsgInd);
   links(ats, 1) = 2;
   sum(ats)
   list1 = (unique(ttsgInd(bts(ats))));
   length(list1)
   
   [ats, bts] = ismember(myTSlinks.b, ttsgInd);
   links(ats, 2) = 2;
   sum(ats)
   length(unique(bts))
   list2 = unique(ttsgInd(bts(ats)));

   myTSGenes = unique([list1, list2]);
   
   % mark the common genes
   [ats, bts] = ismember(myTSlinks.a, cgInd);
   links(ats, 1) = 3;
   sum(ats)
   list1 = (unique(cgInd(bts(ats))));
   length(list1)
   
   [ats, bts] = ismember(myTSlinks.b, cgInd);
   links(ats, 2) = 3;
   sum(ats)
   list1 = (unique(cgInd(bts(ats))));
   length(list1)

   myCGenes =  unique([list1, list2]);

   llabs = sum(links');
   
   % % Just some observations
   % sum(llabs == 4)
   % ml = find(llabs == 0);
   % i = ml(2)
   % TSlinks{tissueID}.b(i)
   % TSlinks{tissueID}.a(i)
   
   % getting the networks...
   wholeNet = sparse(myTSlinks.a, myTSlinks.b, ...
                   1, 18494, 18494);
   
   ml = find(llabs == 4);
   tsNet = sparse(myTSlinks.a(ml), myTSlinks.b(ml), ...
                   1, 18494, 18494);
   length(ml)

   ml = find(llabs == 5);   
   hafnhafNet = sparse(myTSlinks.a(ml), myTSlinks.b(ml), ...
                   1, 18494, 18494);
   length(ml)
   
   ml = find(llabs == 6);   
   cgNet = sparse(myTSlinks.a(ml), myTSlinks.b(ml), ...
                   1, 18494, 18494);

   length(ml)
   % >>> just trying : getting the links with minimum next value.
   % nextNA = TSlinks{t}.maxTSna(na);
   % cgNextNa = nextNA(ml);
   % book = [myTSlinks.a(ml), myTSlinks.b(ml), nextNA(ml)];
   % [a, b] = sort(book(:, 3));
   % newBook = [book(b, 1), book(b, 2), book(b, 3)];
   % f = 1
   % newBook(f, 1:2)
   % f = f + 1
   % <<< just trying
   
   % Identify the pure links using FDR from the RSRed 
   load(['~/resultsAndFigures/expAndCorrRegression/' tissue ...
         'TSReg_UnionFDR0012.mat'])
   tReg = TSReg;
   clear TSReg
   % load the shuffled for the tissue to find FDR
   load(['~/resultsAndFigures/expAndCorrRegression/' tissue ...
        'TSReg_UnionFDR_shuffled0012.mat'])
   rReg = TSReg;
   clear TSReg
   
   ht = hist(tReg.RSEd, [0:.01:.99])/length(tReg.RSEd);
   hr = hist(rReg.RSEd, [0:.01:.99])/length(rReg.RSEd);
   
   for k = 1:100
       fdr(k) = sum(hr(k:end))/sum(ht(k:end));
   end
   
   rsqThr = (min(find(fdr<=.01)))/100

   rsqs(t) = rsqThr;
   
   % [tfi, txi] = ksdensity(tReg.RSEd, 'Support', [0,1]);
   % [rfi, rxi] = ksdensity(rReg.RSEd, 'Support', [0,1]);
   
   % plot(txi, tfi)
   % hold on
   % plot(rxi, rfi)
   
   rsed = tReg.RSEd(na);
   myTSlinks.rsed = rsed;
   %coef1 = tReg.coef1(na);
   
   noTSeff = (rsed < rsqThr);
   %   noTSGTExNet(i).net(g1, g2);eff(end)= [];
   sum(noTSeff)
   
   % the no predictable net
   noTSNet = sparse(myTSlinks.a(noTSeff), myTSlinks.b(noTSeff), ...
                   1, 18494, 18494);
   
   % coefTSNet = sparse(myTSlinks.a(noTSeff), myTSlinks.b(noTSeff), ...
   %                 coef1, 18494, 18494);
   
   % build the networks for the links and the same dataset. 
   
   finalTable(t).wholeNet = wholeNet;
   finalTable(t).TSVNet = sparse(myTSlinks.a, myTSlinks.b, myTSlinks.tsValue, ...
                                 18494, 18494)
   finalTable(t).RSEDNet = sparse(myTSlinks.a, myTSlinks.b, myTSlinks.rsed, ...
                                 18494, 18494)
   finalTable(t).totGCount = totalGeneCount;
   finalTable(t).tissue = tissue;
   finalTable(t).cGenesInd = myCGenes;
   finalTable(t).tsGenesInd = myTSGenes;
   finalTable(t).totalTSGenes = ttsgInd;
    finalTable(t).tsNet = tsNet;
    finalTable(t).hafNhafNet = hafnhafNet;
    finalTable(t).noTSNet = noTSNet;
    finalTable(t).cgNet = cgNet;

    % get the labels net, with zeros as I donno. 
    llabs(llabs == 0) = 1;
    llabsNet = sparse(myTSlinks.a, ...
                         myTSlinks.b, [llabs], 18494, 18494);
    finalTable(t).llabs = llabsNet;
end

save(sprintf('~/resultsAndFigures/TSlinks/finalTable_CG13_NA_EA_FC%slog_FDR%s.mat', FC, FDR), ...
     'finalTable')


save(sprintf('~/resultsAndFigures/TSlinks/finalTable_CG13_FC%slog_FDR%s.mat', FC, FDR), ...
     'finalTable')

clear
FC = '3'
FDR = '0012'
load(sprintf(['~/resultsAndFigures/TSlinks/finalTables/' ...
              'finalTable_CG13_FC%slog_FDR%s.mat'], FC, FDR))

book = finalTable(1).cgNet + finalTable(1).noTSNet;
% 15. 1 saving the general info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
FC = '3'
FDR = '0012'
load(sprintf(['~/resultsAndFigures/TSlinks/finalTables/' ...
              'finalTable_CG13_FC%slog_FDR%s.mat'], FC, FDR))
myft = finalTable;
% getting the network info

totGC = zeros(1, 5);
totEC = zeros(1, 5);
netD =  zeros(1, 5);
inGenes = zeros(1, 18494);
plinks =  zeros(1, 5);
hafnhaf =  zeros(1, 5);
explinks =  zeros(1, 5);
for t = 1:5
    totGC(t) = myft(t).totGCount;
    totEC(t) = sum(sum(myft(t).wholeNet));
    netD(t) = totEC(t) / (totGC(t) * (totGC(t) - 1))/2;
    nd = sum(myft(t).wholeNet') + sum(myft(t).wholeNet);
    inGenes = inGenes + (nd > 0);
    plinks(t) = sum(sum(myft(t).cgNet));
    hafnhaf(t) = sum(sum(myft(t).hafNhafNet));
    explinks(t) = sum(sum(myft(t).tsNet));
end

file = ['~/resultsAndFigures/tissueNetworkStudies/' ...
        'generalNetInfo_TSNet_FDR0012.txt']
fid = fopen(file, 'w')

fprintf(fid, ['attribute\t blood\t brain\t liver\t lung\t skeletalMuscle ' ...
              '\n'])
fprintf(fid, ['%s\t' repmat('%d\t', 1, 5) '\n'],...
            'Network Genes', totGC)
fprintf(fid, ['%s\t' repmat('%d\t', 1, 5) '\n'],...
            'Edge Count', totEC)
fprintf(fid, ['%s\t' repmat('%.5f\t', 1, 5) '\n'],...
            'Network Density', netD)
fclose(fid)

% clear finalTable
% load(sprintf(['~/resultsAndFigures/TSlinks/finalTables/' ...
%               'obsolete_error/finalTable_CG23_FC2_FDR0012.mat']))
p
% 15.2 getting the count of different links and genes. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

clear
FC = '3'
FDR = '0012'
load(sprintf(['~/resultsAndFigures/TSlinks/finalTables/' ...
              'finalTable_CG12_FC%slog_FDR%s.mat'], FC, FDR))
myft = finalTable;
book = myft(2).noTSNet + myft(2).hafNhafNet;
[a, b] = find(book == 2);
length(a)
k = 4
[a(k), b(k)]

save('~/resultsAndFigures/TSlinks/finalTable23_FDR0012.mat', 'finalTable')
save('~/resultsAndFigures/TSlinks/finalTable13.mat', 'finalTable')
save('~/resultsAndFigures/TSlinks/finalTable.mat', 'finalTable')

load('~/resultsAndFigures/TSlinks/finalTable12.mat')
load('~/resultsAndFigures/TSGenes/tsGenes.mat')
load('~/resultsAndFigures/TSlinks/finalTable12_FDR0012.mat')
load('~/resultsAndFigures/TSlinks/finalTables/finalTable23_FDR0012.mat')

save('~/resultsAndFigures/TSlinks/rsqThrs.mat', 'rsqs')

tspMat = vec2mat([tsGenes(:).anp], 18494);

% just adding the node degree for the TS net
%%%%

t = 5 % for t = 1:5
tnet = finalTable(t).wholeNet;
tnd = sum(tnet) + sum(tnet');
finalTable(t).nd = tnd;

%%%%%%%

wholeNet = finalTable(2).wholeNet;

sib = sum(wholeNet);
sib = sib + sum(wholeNet');
[a, b] = sort(sib, 'descend');

bins = zeros(1500, 7);
for i = 1:1500
    kado = wholeNet(b(i), :) + wholeNet(:, b(i))';
    myps = tspMat(2, kado > 0);
    bins(i, :) = hist(log10(myps), [-6:1:0])/length(myps);
end
figure
mypsSorted = sort(myps);
heatmap(bins, [],  [],  [], 'Colorbar', true, 'Colormap', ...
         'summer')

find(bins(:, 1) > .7)

whos sib 
b(1493)

% let's identify the brain specificity of these 1500 genes. 

tind = find(kado);
load('~/data/general/GPL570GemmaMapNEW.mat')

[as, bs] = ismember({'LTBR', 'MAVS', 'TBK1', 'CD40'},gpl570.uniqueSymbols)

kado(b)

% load the blood network to check these geens. 
tissue = 'brain'
expThr = '0.8'
load( ['~/networks/tissues/' tissue '/' ...
       'binaryNet_FDR5e-5_' expThr 'Expr_Ind0.10.mat'])

nn = binNet(b, :) + binNet(:, b)';
sum(nn)

sum((nn + kado)==2)

% this gene is a purely TS gene and therfore so many TS
% links. Also, as my lung data is mixed supposed with inflamatory
% lung, this gene has many links in the lung too. amazingly, no
% overlap with the blood. Something like this gene should be
% checked in GTEx, cause in the gene card seems to have expressione
% everywhere. 

% getting list of random genes
sib = finalTable(5).llabs == 6;
finalTable(5).tissue
[a, b, c] = find(sib);
whos a

r = randi([1 length(a)], 1, 10)

[a(r), b(r)]

% for a given pair of links, give me the RSQred and the anova and
% look at random. 

load('~/resultsAndFigures/TSlinks/TSlinksExtended.mat')

book = TSlinks{2};

pvNet = sparse(book.a, book.b, book.pv);

naNet = sparse(book.a, book.b, book.na);

sib =[     6953       15440;
        1215        2211;
        8861        9995;
          18        2987;
        1229        4861;
        2336       16574;
        2765       15229;
       10343       10599;
       10991       17708;
         785        2336]

i = 7
p = pvNet(sib(i, 1), sib(i, 2))
load(['~/networks/continuousTSLinksResults/' ...
      'random5_smallMat_pVals.mat'])

myps = smallMat(logical(triu(ones(size(smallMat)), 1)));

sum(sum(myps < p))


% another looking into it. 

myNet = finalTable(2).noTSNet;

mySNet = sparse(TSlinks{2}.a, TSlinks{2}.b, TSlinks{2}.na);
mySNet(18494, 18494) = 1;

holu = myNet .* mySNet;

[a, b, c] = find(holu);
[as, bs] = sort(c, 'descend');

b = b(bs);
c = c(bs);
a = a(bs);

i = 5
[a(i), b(i)]


% 16. the part about regression identified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 17. the part about CG genes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}
% load  the links
load('~/resultsAndFigures/TSlinks/tissues.mat')
%load('~/resultsAndFigures/TSlinks/TSlinksExtended.mat')

FDR = '0012'
FC = '3'

load(sprintf(['~/resultsAndFigures/TSlinks/' ...
              'TSlinks_FDR%s_withNextTissueValues.mat'], FDR))

% load the TS and common genes. 
load(sprintf('~/resultsAndFigures/TSGenes/finalTSGenes_FDR.1_FC%slog.mat', ...
             FC));
% load(['~/resultsAndFigures/TSGenes/' ...
%       'TSExpressionAnova.mat'])
% cGenes = p > .1/1849;

% another way of getting the everywhere expressed genes, using the
% genes expressed everywhere and leave the specificity alone.
load('~/data/general/linkExprInfo/dataSetExpInf.mat')
mat12 = zeros(18494, 53);
mat23 = zeros(18494, 53);
mat13 = zeros(18494, 53);
for i = 1:53
    q = quantile(dataSetExpInf(i).meanExp, 5);
    mat12(:, i) = dataSetExpInf(i).meanExp > q(3);
    mat23(:, i) = dataSetExpInf(i).meanExp > q(4);
    mat13(:, i) = dataSetExpInf(i).meanExp > q(2);
end

% sum over the tissues: 
t12 = zeros(18494, 5);
t12(:,1) = sum(mat12(:,1:9)');
t12(:,2) = sum(mat12(:,10:24)');
t12(:,3) = sum(mat12(:,25:31)');
t12(:,4) = sum(mat12(:,32:41)');
t12(:,5) = sum(mat12(:,42:53)');
sib = max(t12);
thr = sib * .8

t23 = zeros(18494, 5);
t23(:,1) = sum(mat23(:,1:9)');
t23(:,2) = sum(mat23(:,10:24)');
t23(:,3) = sum(mat23(:,25:31)');
t23(:,4) = sum(mat23(:,32:41)');
t23(:,5) = sum(mat23(:,42:53)');


t13 = zeros(18494, 5);
t13(:,1) = sum(mat13(:,1:9)');
t13(:,2) = sum(mat13(:,10:24)');
t13(:,3) = sum(mat13(:,25:31)');
t13(:,4) = sum(mat13(:,32:41)');
t13(:,5) = sum(mat13(:,42:53)');


expMat12 = zeros(18494, 5);
expMat23 = zeros(18494, 5);
expMat13 = zeros(18494, 5);
for i = 1:5
    expMat12(:, i) = t12(:, i) >= thr(i);
    expMat23(:, i) = t23(:, i) >= thr(i);
    expMat13(:, i) = t13(:, i) >= thr(i);
end
cGenes = sum(expMat12') == 5;
sum(cGenes)
tempInds = find(cGenes);
cGenes = sum(expMat23') == 5;
sum(cGenes)
% tempInds = find(cGenes);
cGenes = sum(expMat13') == 5;
sum(cGenes)
tempInds = find(cGenes);

% get the atn ND for these genes 
FC = '3'
FDR = '0012'
load(sprintf(['~/resultsAndFigures/TSlinks/finalTables/' ...
              'finalTable_CG13_FC%slog_FDR%s.mat'], FC, FDR))

tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'};
atnds = zeros(5, 18494);
tsnds = zeros(5, 18494);
for t = 1 : 5
    tissue = tissues{t};
    load( ['~/networks/tissues/' tissue '/' ...
           'binaryNet_FDR5e-5_0.8Expr_Ind0.10.mat'])
    atn{t} = binNet;
    atnds(t, :) = sum(binNet) + sum(binNet');
    tsnet = finalTable(t).wholeNet;
    tsnds(t, :) = sum(tsnet) + sum(tsnet');
end

% getting the distri
for t = 1:5
    mytsnds = tsnds(t, :);
    inGenes = mytsnds > 0;
    kado = inGenes + cGenes;
    cgts = hist(mytsnds(kado == 2), [0:20:2000]);
    ndhists(t).cgts = cgts / sum(cgts);
    [fi, xi] = ksdensity(mytsnds(kado == 2), 'Support', [0, 2500]);    
    nddists(t).tsCG = [fi; xi];
    ogts = hist(mytsnds(kado == 1), [0:20:2000]);
    ndhists(t).ogts = ogts / sum(ogts);
    [fi, xi] = ksdensity(mytsnds(kado == 1), 'Support', [-.001, 2500]);    
    nddists(t).tsOG = [fi; xi];
    
    
    myatnds = atnds(t, :);
    inGenes = myatnds > 0;
    kado = inGenes + cGenes;
    cgat = hist(myatnds(kado == 2), [0:20:2000]);
    ndhists(t).cgat = cgat / sum(cgat);
    [fi, xi] = ksdensity(myatnds(kado == 2), 'Support', [0, 2500]);    
    nddists(t).atCG = [fi; xi];
    ogat = hist(myatnds(kado == 1), [0:20:2000]);
    ndhists(t).ogat = ogat / sum(ogat);
    [fi, xi] = ksdensity(myatnds(kado == 1), 'Support', [-0.001, 2500]);    
    nddists(t).atOG = [fi; xi];
end

h = figure
plot([0:20:2000], ndhists(1).cgts, 'k')
hold on
plot([0:20:2000], ndhists(1).cgat)

h = figure
plot([0:20:2000], ndhists(1).ogts, 'k')
hold on
plot([0:20:2000], ndhists(1).ogat)


plot(nddists(1).tsCG(2, :), nddists(1).tsCG(1,:), 'k')
hold on
plot(nddists(1).atCG(2, :), nddists(1).tsCG(1,:))

h = figure
plot(nddists(1).tsOG(2, :), nddists(1).tsOG(1,:), 'k')
hold on
plot(nddists(1).atOG(2, :), nddists(1).tsOG(1,:))

% first try: just by the count of the ts vs atn links for the
% gene. This kind of doesnt work, since the density of the networks
% don't match and I don't know how to normalize for that: just go
% for tissue based, reduce the count of hists to 20. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% In the following, I am looking at the ND distribution of three
%% groups of genes in the tissues: From this, I am hoping to be
%% able to capture the effect of the common genes in the TS
%% networks. 

totTSnds = sum(tsnds);
totATnds = sum(atnds);

genesInInd = (totATnds > 0) + (~cGenes);
genesInCount = length(genesInInd);

% book has the portion of ts links for the genesInInd
book = totTSnds(genesInInd) ./ (totATnds(genesInInd));

% 1. values has the count TS links for genesInInd 
tlpInGenes = zeros(5, genesInCount); % tissue link portion in genes
histInGenes = zeros(5, 20);
for t = 1:5
    tlpInGenes(t, :) = tsnds(t, genesInInd) ./ (atnds(t, genesInInd) ...
                                                  + 1);
    % getting the hist
    histInGenes(t, :) = hist(tlpInGenes(t, :), [.025:.05:1]);
end
h = figure
bar(log10(histInGenes), 'k')
hold on

% 2. same thing for the cGenes tissue link
cGenesInd = find(cGenes);
cGenesCount = length(cGenesInd);
tlpCGenes = zeros(5, length(cGenesInd)); % tissue link portion
                                         % common genes
histCGenes = zeros(5, 10);
for t = 1:5
    % getting the values
    tlpCGenes(t, :) = tsnds(t, cGenesInd) ./ (atnds(t, cGenesInd) ...
                                                  + 1);
    % getting the hist
    histCGenes(t, :) = hist(tlpCGenes(t, :), [.025:.1:1]);
end
whos tlpCGenes
h = figure
bar(log10(histCGenes), 'b')
bar((histCGenes), 'r')

% 3. Also, for all the genes except for the common genes
otherGenes = (totATnds > 0 )- cGenes;
otherGenesInd = find(otherGenes);
otherGenesCount = length(otherGenesInd)
tlpOtherGenes = zeros(5, length(otherGenesInd)); % tissue link
                                                 % portion common
                                                 % genes
histOtherGenes = zeros(5, 10);
for t = 1:5
    tlpOtherGenes(t, :) = tsnds(t, otherGenesInd) ./ (atnds(t, otherGenesInd) ...
                                                  + 1);
    % getting the hist
    histOtherGenes(t, :) = hist(tlpOtherGenes(t, :), [.025:.1:1]);
end
h = figure
bar(log10(histOtherGenes), 'k')
hold on

title('to be explained')
ylabel('log10 count of genes')
set(gca, 'XTickLabel', {'blood', 'brain', 'liver', 'lung', 'muscle'})
file = '~/resultsAndFigures/contribution.pdf'

print(h, '-dpdf', file)
% NOTE: so I have the bar plot of 2 over 3, which presenets how
% little these genes contribute to tissue specific
% rewiring (the hist plots should be normalized). However, among these genes, we looked at thos who did
% have multiple TS partners : get the genes labeled belonging to
% each group for each tissue, then study the genes interesting. 

ndinfo.atnds = atnds;
ndinfo.tsnds = tsnds;
ndinfo.cgenes = cGenes;

save('~/resultsAndFigures/tissueNetworkStudies/ndinfo.mat', 'ndinfo')

%% Get the box plots for the cgenes , for each tissue
%
% which group does it belong to
groups = zeros(size(tlpCGenes));
s = 0;
for i = 1:20
    e = s + .05;
    for t = 1:5
        temp = (tlpCGenes(t, :) > s) + (tlpCGenes(t, :) <= e);
        myIns = (temp > 1);
        groups(t, myIns) = i;
    end
    s = e;
end
groups(groups == 0)= 1;

% the box plots show  how little do they contribute to the TS
% network. 
values = tsnds(:, cGenes);
for t = 1:5
    h = figure
    boxplot(values(t, :), groups(t, :))
    title([tissues{t} ' common'])
    xlim([0, 21])
end

% count of TS linkso
tempC = values > 40;

% portion of TS links
tempP = tlpCGenes > .30;

potential = tempC + tempP;

[a, b, c] = find(potential == 2);
length(a)
% so I found 259 genes : what are they? 
myGenes = cGenesInd(unique(b)); 

save('~/resultsAndFigures/pureGeneEnrichment/hitList.mat', 'myGenes')

%% Get the box plot for the other genes. Notice that we care
%% for the trend here, not the counts : but, those with the less
%% ATN and more TSN are cool: notice that I am missing a lot of
%% them because of the expression thing: get the new networks!!! 
% which group does it belong to
groups = zeros(size(tlpOtherGenes));
s = 0;
for i = 1:20
    e = s + .05;
    for t = 1:5
        temp = (tlpOtherGenes(t, :) > s) + (tlpOtherGenes(t, :) <= e);
        myIns = (temp > 1);
        groups(t, myIns) = i;
    end
    s = e;
end
groups(groups == 0)= 1;

% the box plots show  how little do they contribute to the TS
% network. 
values = tsnds(:, logical(otherGenes));
for t = 1:5
    h = figure
    boxplot(values(t, :), groups(t, :))
    title([tissues{t} ' other'])
    xlim([0, 21])
end

%%%%% Done with the plots. 

% hist 1
h = figure;
hist(book, 100)
barInGenes = hist(book, [0:.01:1]);

% getting the box plot of the TSN links
% count of ts links
ovalues = totTSnds(genesInInd);
% which group does it belong to
groups = zeros(1, length(genesInInd));
for i = 1:100
    s = (i - 1) /100
    e = i /100
    temp = (book > s) + (book <= e);
    myIns = (temp > 1);
    groups(myIns) = i;
end
groups(groups == 0)= 1;
h = figure
boxplot(values, groups)

% all the common genes have ATN links
cgenesInd = find(cGenes);
book1 = totTSnds(cGenes) ./ totATnds(cGenes);
barcGenes = hist(book, [0:.01:1]);
hold on

% hist 2
h = figure;
hist(book1, 100)

values1 = totTSnds(cGenes);
% which group does it belong to
groups1 = zeros(1, sum(cGenes));
for i = 1:100
    s = (i - 1) /100
    e = i /100
    temp = (book1 > s) + (book1 <= e);
    myIns = (temp > 1);
    groups1(myIns) = i;
end
groups1(groups1 == 0)= 1;
boxplot(values1, groups1)

cGenesHighTS = find(groups1 > 30)
cGenesInd = find(cGenes);
myfinal = cGenesInd(cGenesHighTS);
finalMat = [tsnds(:, myfinal); atnds(:, myfinal)]
figure
heatmap(finalMat, [],  [],  [], 'Colorbar', true, 'Colormap', ...
         'bone')

heatmap(finalMat);

h = figure;
hist(totATnds(genesIn), 100)

h = figure;
hist(totTSnds(cGenes), 100)

h = figure;
hist(totATnds(cGenes), 100)

tsGenesIn = totTSnds > 0;
hist(totTSnds(tsGenesIn), 100)

% my list of genes: cGenes. 

% Draft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
load('~/data/general/GPL570GemmaMapNEW.mat')
load('~/data/general/ribosomeIDs.mat')




