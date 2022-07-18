% In this file I am getting the genes which are expressed in 80% of
% the datasets in a tissue. 

clear
load('~/data/general/linkExprInfo/dataSetProbeInf.mat')
load('~/data/general/linkExprInfo/wholeGeneExpr.mat')
load('~/data/general/corrBins_expStd.mat')

% getting the tissue indexes
bloodInd = [1:9];
lungInd = [10:24];
smInd = [25:31];
smInd = [32:41];
brainInd = [42:53];

gCount = length(wholeGeneExpr);
qExpr = 2 .* ones(gCount, 53);
for i = 1:length(dataSetProbeInf)
    i 
    s = dataSetProbeInf(i).sampleInd(1);
    e = dataSetProbeInf(i).sampleInd(2);
    expMat = wholeGeneExpr(:, s:e);

    avgExpr = mean(expMat');
    stdExpr = std(expMat');
    
    kado = [0:1/6:1]
    qCount = 6;
    Qs = quantile(avgExpr, kado)
    qSize= floor(gCount/qCount)

    addVal = [1, 2, 2, 2, 2];
    for j = 2:qCount
        tempInd = avgExpr >= Qs(j);
        qExpr(tempInd, i) = qExpr(tempInd, i) + addVal(j-1);
    end
end

tissue = 'skeletalMuscle'
indX = smInd;
percs = [.4, .6, .8, 1]
npercs = 1/3;

for i =1:length(percs)

    book = qExpr(:, indX);
    sib = book > 3;

    holu = sum(sib');

    dsCount = length(indX);
    expDSC = ceil(percs(i) * dsCount)

    nonExpDSC = floor(npercs * dsCount)

    expCounts = sum(holu >= expDSC )
    expGenesInd = holu >= expDSC;

    nonExpCounts = sum(holu <= nonExpDSC);
    nonExpGenesInd = holu <= nonExpDSC;

    save(sprintf('~/data/general/tissueExpGenes/%sExpGenes%.1f.mat', tissue, percs(i)), ...
         'expGenesInd')
    save(sprintf('~/data/general/tissueExpGenes/%sNonExpGenes%.1f.mat', ...
                 tissue, npercs), 'nonExpGenesInd')
end

load('~/data/general/hkgInd.mat')

% getting and saving genes which were expressed in all the tissues,
% 80%
tissues = {'blood', 'liver', 'lung', 'skeletalMuscle', 'brain'}
allExp = zeros(18494, 1);
for i = 1:5
    tissue = tissues{i};
    load(['~/data/general/tissueExpGenes/' tissue 'ExpGenes0.8.mat'])
    allExp = allExp + expGenesInd';
end
allExp = allExp == 5;
sum(allExp)

save('~/data/general/tissueExpGenes/allExp0.8.mat', 'allExp')

% here, check if jackknife of expression changes in different tissues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% getting random set of DS for the continuousTSLinksRandomFunction.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% I am gonna save the random sets and save the ds IDs in a file.

% random sets are generated in a file, I generate 5 random sets 

bloodInd = [1:9];
lungInd = [10:24];
liverInd = [25:34];
smInd = [35:41];
brainInd = [42:53];

dsCounts = [7, 12, 11, 9, 15];

for i = 1:10
    
    temp = remeinder(i, 5) + 1;
    dCount = dsCounts(temp);

    % selecting the datasets
    dslist = zeros(1, dsCount);

    % get the div of the dcount / 5 and rem, for each tissue get
    % the count of datasets you need from that tissue. 

    expGenes = zeros(1, size(wholeGeneExpr, 1));
    for i = 1:length(dsIDlist)
        d = dsIDlist(i);
        dsMat = wholeGeneExpr(:, dataSetProbeInf(d).sampleInd(1): ...
                              dataSetProbeInf(d).sampleInd(2));
        avgExp = mean(dsMat');
        qExp = quantile(avgExp, [0:1/6:1]);
        tempExp = avgExp > qExp(3);
        expGenes = expGenes + tempExp;
    end

    expGenesInd = expGenes > (length(dsIDlist) * .8);
    sum(expGenes)
    randDSList.dsList = list;
    randDSList.expGenes = expGenesInd;
    load(sprintf('~/data/general/linkExprInfo/randDSList_%d.mat', ...
                 i), 'randDSList');
end


% getting 2/3 of genes for all the datasets, so I can select
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
addpath('~/codes/MATLAB/myCodes/general/')
tissues = {'brain', 'blood', 'skeletalMuscle', 'liver', 'lung'}

% for each tissue 
c = 1;
for t = 1:length(tissues)

    files = dir(['~/data/affyArray/tissues/' tissues{t} '/' ...
                        'matFiles/geneExprMatGemmaMapBER/*.mat'])
    for f = 1:length(files)
        load(['~/data/affyArray/tissues/' tissues{t} ...
              '/matFiles/geneExprMatGemmaMapBER/' files(f).name])
        GSEID = getGSEIDfromStr(files(f).name)
        sib = mean(dataSet.mat');
        q = quantile(sib, 1/3);
        expr = sib > q;
        sum(expr)
        
        allDSExprInfo(c).tissue = tissues{t};
        allDSExprInfo(c).GSEID = GSEID;
        allDSExprInfo(c).exprLevel = sib;
        allDSExprInfo(c).exprGenes = expr;

        c = c + 1;
    end
end

save('~/data/general/tissueExpGenes/allDSExprInfo.mat', 'allDSExprInfo')
 



