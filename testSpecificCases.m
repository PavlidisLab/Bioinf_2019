% 1. Testing individual cases in GTEx data, in GO data and also
% closer looks into networks. 

% 1. Testing individual cases in GTEx data, in GO data and also
% closer looks into networks. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
load('~/networks/GTEx/allNets.mat')

FC = '3'
FDR = '0012'
load(sprintf(['~/resultsAndFigures/TSlinks/finalTables/' ...
              'finalTable_CG13_FC%slog_FDR%s.mat'], FC, FDR))
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'};

gCount = 18494;

load('~/data/general/GOdata_GPL570_07_nonComp.mat')
%load('~/data/general/GOdata_GPL570_07.mat')
sumF = sum(GOdata.matF);
sumP = sum(GOdata.matP);
sumC = sum(GOdata.matC);
geneCounts = sumP;% + sumC;% + sumF;

% removing the functions which has < 3 genes
inFunctions1 = geneCounts >= 10;
inFunctions2 = geneCounts <= 200;
inFunctions = (inFunctions1 + inFunctions2) == 2;
inFInds = find((inFunctions1 + inFunctions2) == 2);
fCount = sum(inFunctions)
inTerms = GOdata.GOTerms(inFInds);

load('~/resultsAndFigures/functionalIdentityForGenes/sofUtil.mat')
fMat = SOFutil.fMat;
tsfs = SOFutil.tsfs;
myGenes = SOFutil.myGenes;

for t = 1 : 5
    tissue = tissues{t};
    load( ['~/networks/tissues/' tissue '/' ...
           'binaryNet_FDR5e-5_0.8Expr_Ind0.10.mat'])
    atn{t} = binNet;
end

% looking at ATF2
geneID = 1137;
termID = 358;

% how many genes have that term in fMat
sum(fMat(:, termID))
posGenes = find(fMat(:, termID));
posGenes = [3456, 3458, 3459, 3788, 13053, 1169, 1170]'

posGenes = find(geneNeighs(2, 12965:13066)) + 12964;
posGenes = find(radfamily == 2)' + 12964;
posGenes = posGenes';

% hamny genes have that term in AFT2 neighbours

geneNeighs = zeros(5, 18494);
for t = 1:5
    net = atn{t};
    geneNeighs(t, :) = net(:, 1137) + net(1137, :)';
end

sib = geneNeighs(:, posGenes);
sum(sib')

tinyNetGenes = [1137 posGenes']
net = atn{5};
GOdata.geneSymbols(tinyNetGenes)'
tinyNet = full(net(tinyNetGenes, tinyNetGenes))


gtexGeneNeighs = zeros(17, 18494);
for t = 1:17
    net = GTExNet(t).net;
    gtexGeneNeighs(t, :) = net(:, 1137) + net(1137, :)';
end

sib = gtexGeneNeighs(:, posGenes(6));
subGtexNeighs = gtexGeneNeighs([17, 5, 14 15 16], : );
sum(sib')

% but the actual correaltion values for GTEx for these genes: each
% tissue gets a distribution bins. 
addpath('~/codes/MATLAB/myCodes/general')
finalMat = zeros(length(posGenes), 17);
for t = 1:17
    fName = GTExNet(t).tissue;
    load(['~/data/GTEx/matFormat_GPL570/' fName '.mat'])
    dsGT = dataSet.mat;
    size(dsGT)
    smallDS = dsGT([1137, posGenes'], :);
    sib = corr(smallDS');
    
    myCorrs = GTExNet(t).corrQs;
    
    corrQMat = 500 .* ones(1+length(posGenes), 1+length(posGenes));
    for i = 1:(1+length(posGenes))
        for j = (i + 1):(1+length(posGenes))
            value = sib(i, j);
            if (value < myCorrs(1000))
                vrank = min(find(myCorrs > value));
            else
                vrank = 1000;
            end
            corrQMat(i, j) = vrank;
        end
    end
    
    finalMat(:, t) = corrQMat(1, 2:end);
    
    % h = figure;
    % heatmap(corrQMat, [], [], [], 'TickAngle', 45, 'Colorbar', true, ...
    %         'Colormap', 'gray')
end

h = figure
boxplot(finalMat(:, [17, 7,14:16]))
boxplot(totalMedsCorrVals(:, [17, 7, 14:16]))
ylim([-.3, 1])
title('GTEx')

figFolder = ['~/resultsAndFigures/functionalIdentityForGenes/' ...
             'figures/'']
fileName = ['ATF2_5repairthings_GTEx_box_corr'];
fileName = ['ATF2_17rabrad_GTEx_box_corr'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);
 

rows = ([1, 2, 3, 6, 10,5,7,4,  9 12 11 8 ])
 heatmap(finalMat(rows, [6 7 8, 17,14:16]), [], [], [], 'TickAngle', 45, 'Colorbar', true, ...
             'Colormap', 'gray')

posGenes(rows)

%microarrayMedians: ()
% for the DNA damage response 
affyMeds = [737, 994, 950, 891, 988; 
           904, 992, 982, 747, 900;
           911, 979, 959, 649, 873;
           427, 962, 645, 565, 883;
           957, 908, 956, 850, 985;
           949, 989, 961, 980, 999;
           303, 879, 574, 439, 726;
           593, 265, 346, 441, 693;
           235, 661, 385, 693, 525;
           876, 992, 971, 795, 990;
           728, 696, 833, 693, 720;
           251, 732, 631, 529, 574]

h = figure;
affyMeds = affyMeds([1 2 3 6 10], :);
boxplot(affyMeds)
boxplot(totalMedsCorrVals)
ylim([-.2, 1])
figFolder = '~/resultsAndFigures/functionalIdentityForGenes/figures/'
fileName = ['ATF2_5damageResponseGenes_affy_box_corr'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);


affyMeds = [717, 979, 783, 903, 970;
           929, 981, 840, 820, 986;
           950, 978, 987, 956, 974;
           809, 987, 973, 954, 988;
           733, 995, 834, 842, 898;
           981, 983, 468, 616, 230;
           898, 987, 825, 875, 977;
           903, 987, 906, 903, 911;
           849, 991, 963, 979, 984;
           974, 974, 858, 946, 961;
           869, 957, 610, 546, 809;
           968, 962, 936, 799, 965;
           876, 992, 971, 795, 990]

h = figure;
boxplot(affyMeds)
boxplot(totalMedsCorrVals)
ylim([-.2, 1])
figFolder = '~/resultsAndFigures/functionalIdentityForGenes/figures/'
fileName = ['ATF2_17RabRad_affy_box_corr'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);

% getting the ranks back to one dataset 

h = figure
 heatmap(affyMeds(rows, [2, 1, 3 4 5]), [], [], [], 'TickAngle', 45, 'Colorbar', true, ...
             'Colormap', 'gray')
 
totalMeds = [affyMeds, finalMat(:,  [17, 7,14:16])];

load('~/data/general/corrBins_expStd.mat')
% getting ranks to values... 

transferVals = affyMeds;
transferVals = finalMat;
totalMedsCorrVals = zeros(size(transferVals));
for i = 1:size(transferVals, 1)
    for j = 1:(size(transferVals, 2))
        totalMedsCorrVals(i, j) = dsCorrQInf(14).totalQ(transferVals(i,j));
    end
end

h = figure
 heatmap(totalMedsCorrVals(rows, :), [],GOdata.geneSymbols(posGenes(rows)) , [], 'TickAngle', 45,'ShowAllTicks', true, 'Colorbar', true, ...
             'Colormap', 'gray')
figFolder = '~/resultsAndFigures/functionalIdentityForGenes/figures/'
fileName = ['ATF2_12damageResponseGenes_Affy_GTEx_corr'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);
 
h = figure
 heatmap(totalMeds(rows(1:5), :), [],GOdata.geneSymbols(posGenes(rows)) , [], 'TickAngle', 45,'ShowAllTicks', true, 'Colorbar', true, ...
             'Colormap', 'gray')
figFolder = '~/resultsAndFigures/functionalIdentityForGenes/figures/'
fileName = ['ATF2_5damageResponseGenes_Affy_GTEx_Bins'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);

% at this point, I was curious to get the partners of the pos genes
% for that thing in all the networks, to see if their relatively
% high correlation values in the other tissues is the reuslt of
% partial correlation and not actual correlation

% ATMIN and ATM are not in the list: 1169, 1170
posGenes = [1137, 3456, 3458, 3459, 3788, 13053]'    

for k = 1: length(posGenes)
    tNet = zeros(18494, 5);
    for t = 1:5
        net = atn{t};
        tNet(:, t) = net(posGenes(k), :)' + net(:, posGenes(k));
    end
    geneTissueNeighs{k} = tNet;
end


for t = 1:5
    gNet = zeros(18494, length(posGenes));
    for k = 1:length(posGenes)
        net = atn{t};
        gNet(:, k) = net(posGenes(k), :)' + net(:, posGenes(k));
    end
    tissueGeneNeighs{t} = gNet;
end

sum(sum(tissueGeneNeighs{5}))

% >> how many genes in total? 1303, It is ok. it is bearable.
% >> do they have common neighbors, and is it shared with AFT2
% let's look at this in a network: 

% relatively common partners in all the tissues for each of them, 
% and does brain have a whole different cluster? so basically, what
% are the neighbors in other tissues that are not shared with
% brain. 

allTNeighs = zeros(18494, 6);
for k = 1:6
    allTNeighs(:, k) = sum(geneTissueNeighs{k}(:, [1,3 4 5]), 2);
end

sib = allTNeighs >= 3; % in most tissues, connected. 

kado = sum(sib(:, 2:end), 2); 

hist(kado(kado > 0)) % for most of them ... 

find(kado == 4) 



