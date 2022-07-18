% In this code I do all the side tests and explorations for
% continuous ts link selection. 
%
%% Part 1: pilot phases for the measurement and measuring the TS
%% for the links
% 
% 1. First step: trying out my algorithm for some random TS
% links. Just wanted to check if it behaves as I expected and code
% the core of the algorithm. 
%
% 2. Trying the code for some real links between HKG - also a
% shoortcut for making random finalMat - it is different than
% the
% 1st part's random finalMat, that you don't get to choose values
% for each tissue. I used it for checking the "everywhere
% coexpressed" or "everywhere absent" space value for links. 
%
% 3. Trying the same code differently - I get the correlations by
% for randomly selected pairs of genes. (not the first part!)
%
% 4. assembling the results and some other checkings - this is
% where I study the everywhere potential links
%
% 5. finding final mat and space for a given pair of genes and
% tissue
%
% 6. getting links expressed in all the tissues - between them,
% getting the TS links.
%
% 6.5. getting the TS links for each tissue, the links under the
% three categoreis
% 
% 7. finding final mat for a given set of g1 and g2
%
% 77. Doing a function for a whole tissue
%
%% Part 2: Determining the THR value for identification of the
%% links and putting the cut off for the tissues. 
%
% 8. TS for a random list of networks: generating random list of datasets for testing and testing
% them. (code for calling the function is here)
%
% 9. Cut off the Threshold: For each link passed the thr, what is
% the TS value for the next best tissue? Get this for all the tissues.
% 
% 10. load the extended TS links and the previous ts links and
% check the correlation of the two values. 
% % 11. put together the results from psuedo tissues and do sth. 
%
% 12. Get the hist for the various neg adjusted
% 
% 13.% 13. Get the tissue negAdjusted and pvalues. Also, modified pvalues
%
% 14. get the thresholds for the union. 
% 15. get the average values

% 1. Trying out my algorithm for TS links - the random part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
% count of datasets in that tissue
tc = 10;

% count of other values 
oc = 40;

ec = 1000;
n = 100
pn = 100
prandMat = zeros(pn, tc + oc + ec);
space = zeros(pn, 1);
for p = 1:pn
    % a random seq
    big = randi(n, 1, oc);
    % the ts of this link, 
    small = randi(n, 1, tc);

    % space 
    % the count is up to 1000 + tc + oc
    total = ec + tc + oc;

    % for each step, I should count the space under the curve. 
    allseq = zeros(1, total);
    allseq(big) = 1;
    for i = 1:length(big)
        allseq(big(i)) = allseq(big(i)) * 2; 
    end

    for i = 1:length(small)
        if allseq(small(i)) == 0
            allseq(small(i)) = 3;
        else
            allseq(small(i)) = allseq(small(i)) * 3; 
        end
    end

    perMat(p, :) = allseq;
    emptyUnit = 1/ec;
    ts = 0;
    tsUnit = 1/tc; % 1/tc
    nontsUnit = 1/oc;
    ss = 0;
    steps = 1;
    tsEmptySurface = 0;
    while steps <= total
        current = allseq(steps);
        if(current)  == 0 
            tsEmptySurface = tsEmptySurface + emptyUnit * ts;
            steps = steps + 1;
        else
            f = factor(current);
            for j = 1:length(f)
                if f(j) == 2
                    ss = ss + (tsEmptySurface * nontsUnit);
                    steps = steps + 1;
                    steps;
                else
                    ts = ts + tsUnit;
                    steps = steps + 1;
                end
            end
        end
    end
    space(p)= ss;
end

% 2. Trying the code for some real links between HKG. Since the
% links are real, I need to build the correlation rank matrix first
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
load('~/data/general/linkExprInfo/wholeGeneExpr.mat')
load('~/data/general/linkExprInfo/dataSetProbeInf.mat')
load('~/data/general/hkgInd.mat')

% count of datasets in that tissue
tic
finalMat = zeros(1000, 53);
for i = 1:53
    i
    dsName = dataSetProbeInf(i).name;
    dsMat = wholeGeneExpr(:, dataSetProbeInf(i).sampleInd(1): ...
                          dataSetProbeInf(i).sampleInd(2));
    hkgMat = dsMat(hkgInd(1:1000), :);
    sib = corr(hkgMat');
    upOnes = triu(ones(size(sib)), 1) .*sib;
    clear sib
    [a, b, c] = find(upOnes);
    finalMat(:, i) = corrQuantileReturn(dataSetProbeInf(i).name, c);
end
toc


% >>>>>>>>> a bypass code for making the final matrix at random values
n = 1;
finalMat = randi(200, n, 53);
finalMat = finalMat + 800;


% >>>>>>>>> a bypass code for making the final matrix at random
% values
n = 2
f1 = randi(1000, 1, 10);
f2 = randi(1000, 1, 43);
finalMat = [f1 f2]
finalMat = [finalMat; finalMat];
finalMat(2, :) = finalMat(2, :) - 60;
% <<<<<<<<<

% other values for eihter random finalMat or real made one. 
%tissue = 'brain'
%tissueInd = [42:53];
space = zeros(n, 1);
tissueInd = [1:10];
otherInd = [1:53];
otherInd(tissueInd)= [];
tc = 10;
oc = 43;
ec = 947;
tic
for p = 1: n
    % % a random seq
    % big = randi(n, 1, oc);
    % % the ts of this link, 
    % small = randi(n, 1, tc);
    % % space 
    % % the count is up to 1000 + tc + oc

    small = finalMat(p, tissueInd);

    big = finalMat(p, otherInd);

    total = ec + tc + oc;

    % for each step, I should count the space under the curve. 
    allseq = zeros(1, total);
    allseq(big) = 1;
    for i = 1:length(big)
        allseq(big(i)) = allseq(big(i)) * 2; 
    end

    for i = 1:length(small)
        if allseq(small(i)) == 0
            allseq(small(i)) = 3;
        else
            allseq(small(i)) = allseq(small(i)) * 3; 
        end
    end

    emptyUnit = 1/ec;
    ts = 0;
    tsUnit = 1/tc; % 1/tc
    nontsUnit = 1/oc;
    ss = 0;
    steps = 0;
    tsEmptySurface = 0;
    while steps < total
        current = allseq(1000 - steps);
        if(current)  == 0 
            tsEmptySurface = tsEmptySurface + emptyUnit * ts;
            steps = steps + 1;
        else
            f = factor(current);
            for j = 1:length(f)
                if f(j) == 2
                    ss = ss + (tsEmptySurface * nontsUnit);
                    steps = steps + 1;
                    steps;
                else
                    ts = ts + tsUnit;
                    steps = steps + 1;
                end
            end
        end
    end
    space(p)= ss;
end
toc

space1 = space; % doing the blood
space2 = space; % doing the liver
space3 = space; % the lung
space4 = space; % the sm
space5 = space; % brain

SS = [space1, space2, space3, space4 , space5];

kado = max(SS');

% >>> checking values for space. I saved it for brain and made the
% priliminary plots to see how do the top hits look. 

[sortedSpace, sortedb] = sort(space, 'descend');
sortedFinalMat = finalMat(sortedb, :);

brainMed = median(sortedFinalMat(:, 42:end)');
otherMed = median(sortedFinalMat(:, 1:41)');

sorteda(20000)

book = hist(space, 40);

[hkgInd(a(sortedb(1:10)))' , hkgInd(b(sortedb(1:10)))']

result.geneInds = [a, b];
result.finalMat = finalMat;
result.space = space;

medians.brain = brainMed;
medians.other = otherMed;
save('~/data/general/continuousTSLinks/brainMedians.mat', ...
     'medians')

save('~/data/general/continuousTSLinks/brainResult.mat', 'result')

% 3. Trying the same code differently - I get the correlations by
% for randomly selected pairs of genes. (not the first part!)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

load('~/data/general/linkExprInfo/wholeGeneExpr.mat')
load('~/data/general/linkExprInfo/dataSetProbeInf.mat')
load('~/data/general/hkgInd.mat')

sampleIndsCell = {dataSetProbeInf(:).sampleInd}
sampleInds = cell2mat(sampleIndsCell');

tic
finalMat = zeros(1000000, 53);
tissue = 'brain'
load(['~/data/general/' tissue 'ExpGenes.mat']);
s = 1
for i = 1:53
    tic
    i
    dsName = dataSetProbeInf(i).name;
    dsMat = wholeGeneExpr(:, dataSetProbeInf(i).sampleInd(1): ...
                          dataSetProbeInf(i).sampleInd(2));
    dsMat = dsMat(brainGenesInd, :);
    sib = corr(dsMat');
    upOnes = sib(logical(triu(ones(size(sib)), 1)));
    clear sib
    arr = upOnes(s:s+999999);
    clear upOnes
    finalMat(:, i) = corrQuantileReturn(dataSetProbeInf(i).name, ...
                                       arr);
    toc
end
toc

tissue = 'brain'
tissueInd = [42:53];
space = zeros(1000, 1);
n = 1000;
tc = 12;
oc = 41;
ec = 947;
g1ID = randi(18494, n);
g2ID = randi(18494, n);
tic
for p = 1: n
    % % a random seq
    % big = randi(n, 1, oc);
    % % the ts of this link, 
    % small = randi(n, 1, tc);

    % % space 
    % % the count is up to 1000 + tc + oc
    
    % building the correlation values for two genes. 
    corrArray = zeros(1, 53);
    for i = 1:53
        g1Exp = wholeGeneExpr(g1ID(p), sampleInds(i, 1):sampleInds(i, ...
                                                          2));
        g2Exp = wholeGeneExpr(g2ID(p), sampleInds(i, 1):sampleInds(i, ...
                                                          2));
        corrArray(i) = corr(g1Exp', g2Exp');
    end
    
    small = corrArray(42:53);

    big = corrArray(1:12);

    total = ec + tc + oc;

    % for each step, I should count the space under the curve. 
    allseq = zeros(1, total);
    allseq(big) = 1;
    for i = 1:length(big)
        allseq(big(i)) = allseq(big(i)) * 2; 
    end

    for i = 1:length(small)
        if allseq(small(i)) == 0
            allseq(small(i)) = 3;
        else
            allseq(small(i)) = allseq(small(i)) * 3; 
        end
    end

    emptyUnit = 1/ec;
    ts = 0;
    tsUnit = 1/tc; % 1/tc
    nontsUnit = 1/oc;
    ss = 0;
    steps = 0;
    tsEmptySurface = 0;
    while steps < total
        current = allseq(1000 - steps);
        if(current)  == 0 
            tsEmptySurface = tsEmptySurface + emptyUnit * ts;
            steps = steps + 1;
        else
            f = factor(current);
            for j = 1:length(f)
                if f(j) == 2
                    ss = ss + (tsEmptySurface * nontsUnit);
                    steps = steps + 1;
                    steps;
                else
                    ts = ts + tsUnit;
                    steps = steps + 1;
                end
            end
        end
    end
    space(p)= ss;
end
toc

% 4. assembling the results and some other checkings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for each tissue, I get a bunch of files. 
clear

tissue = 'blood'
who
% load all the brain files one by one 
load(['~/data/general/tissueExpGenes/' tissue 'ExpGenes0.8.mat']);

tissueGeneCount = sum(expGenesInd)
smallMat = zeros(tissueGeneCount, tissueGeneCount);


% modify this reading as appropriate
files = dir(['~/networks/continuousTSLinksResults/' ...
             'qnorms/'  tissue '/' tissue ...
             '_*qnorms*.mat'])

files = dir(['~/networks/continuousTSLinksResults/' ...
             'negAdjusted_ExpAdjusted/'  tissue '/' tissue ...
             '_*NegAdjusted*.mat'])

files = dir(['~/networks/continuousTSLinksResults/' tissue '_*NegAdjusted.mat'])

files = dir(['~/networks/continuousTSLinksResults/One*' tissue '_*' ...
             'pVals.mat'])

secondTissue  = 'skeletalMuscle'
files = dir(['~/networks/continuousTSLinksResults/' tissue '*_JK_' ...
             secondTissue '.mat'])

files = dir(['~/networks/continuousTSLinksResults/' tissue '/OneSided' tissue '_*' ...
             'pVals.mat'])

% fill in the small tissue matrix 
for i = 1:length(files)
    i
    kado = strsplit(files(i).name, {'_', '.'})
    s = str2num(kado{2});
    e = str2num(kado{3});
    
    load(['~/networks/continuousTSLinksResults/' ...
          'qnorms/' tissue '/' files(i).name])

    % load(['~/networks/continuousTSLinksResults/' ...
    %          'negAdjusted_ExpAdjusted/'  tissue '/' files(i).name])

    %load(['~/networks/continuousTSLinksResults/' files(i).name]);
    pVals = space;
    sl = length(pVals);
    
    % finding the beginning in the mat
    col = 2;
    for j = 1:tissueGeneCount
        if s > j
            s = s - j;
            col = col + 1;
        else
            row = s;
        end
    end
    
    % filling up the matrix! 
    c = 0;
    while (c < sl)
        c = c + 1;
        smallMat(row, col) = pVals(c);
        if (row < (col - 1))
            row = row + 1;
        else
            col = col + 1;
            row = 1;
        end
    end
end

% for the raw neg adjusted  
%sib = smallMat ./ (43 * 10 * 500);

% for the quantile normalized 
sib = smallMat ./ (41 * 12 * .83);

smallMat = sib;
clear sib
% extend to the large tissue matrix 
gCount = 18494;
bigMat = zeros(gCount, gCount);
sib = find(expGenesInd == 1);

% US! cause it was the upper triangle, I donno why! 
for i = 1:tissueGeneCount
    for j = i+1: tissueGeneCount
        bigMat(sib(i), sib(j)) = smallMat(i, j);
    end
end

spaceMat = sparse(bigMat);

% modify the file name as you wish. 

save(['~/networks/continuousTSLinksResults/qnorms/' tissue '/SpaceMat_NA_qnorms.mat'], ...
     'spaceMat')

save(['~/networks/continuousTSLinksResults/negAdjusted_ExpAdjusted/' tissue '/SpaceMat_NA_EA.mat'], ...
     'spaceMat')

save(['~/networks/continuousTSLinksResults/' tissue 'SpaceMat_pVals_oneSided.mat'], ...
     'spaceMat')

save(['~/networks/continuousTSLinksResults/' tissue 'SpaceMat_pVals.mat'], ...
     'spaceMat')

save(['~/networks/continuousTSLinksResults/' tissue 'SpaceMat_negAdjusted.mat'], ...
     'spaceMat')

load(['~/networks/continuousTSLinksResults/brain/brainSpaceMat_negAdjusted.mat'])


save(['~/networks/continuousTSLinksResults/' tissue 'SpaceMat_JK_' secondTissue '.mat'], ...
     'spaceMat')
% 4.5 some random stuff 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bloodMat = bigMat;
bloodNet = wholeNet .* bloodMat;

kado = bloodNet > 0;

values = bloodNet(logical(kado));

[a, b, c] = find(bloodNet);

[ss, sa] = sort(c, 'descend');

g1s = a(sa);
g2s = b(sa);
finalValues = ss;
d = 10

[g1s(1:d), g2s(1:d), finalValues(1:d)]
tsgCount = sum(expGenesInd);
realValues = smallMat(logical(triu(ones(tsgCount, tsgCount))));
realValues = realValues(1:50000000);
book = hist(realValues, 100)
book = quantile(realValues, 1000);
sum(realValues > 0.65)
bigMat = sparse(bigMat);

[a, b] = find(bigMat >= 0.65);

s = bigMat(bigMat >= .65);
[ssor,sInd]= sort(s, 'descend');

[ssor(1:10) sInd(1:10)]

d = 2101
[a(d) b(d)]
i = d
bigMat(a(i), b(i))

% I will get the edges with median > 800 as TS. 

% I will have a part about links which are coexpressed everywhere -
% ribosome links and also the everywhere expressed links in my
% datasets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finding the correlation bins for a specific link - and also the
% space
% checked it for blood links, seems to be working! 
%%%%% 
clear
tissue = 'blood'

load(['~/data/general/' tissue 'ExpGenes.mat']);
gCount = sum(expGenesInd)
(gCount * (gCount - 1))/2

load(['~/networks/continuousTSLinksResults/' tissue 'SpaceMat.mat'])

load('~/data/general/linkExprInfo/wholeGeneExpr.mat')
load('~/data/general/linkExprInfo/dataSetProbeInf.mat')
load('~/data/general/GPL570GemmaMapNEW.mat')

spaceMat(1:10, 1:10)

n = c%1; % count of edges
finalMat = zeros(n, 53);
thisG1s = thisG1s;%g1 = 79
thisG2s = thisG2s;%g2 = 409
coG1s;
coG2s;

for j = 1:n
    j
    for i = 1:53
        i;
        dsName = dataSetProbeInf(i).name;
        dsMat = wholeGeneExpr(:, dataSetProbeInf(i).sampleInd(1): ...
                              dataSetProbeInf(i).sampleInd(2));
        %dsMat = dsMat(expGenesInd, :);
        %sib = corr(dsMat(1:100, :)');
        sib = corr(dsMat(thisG1s(j), :)', dsMat(thisG2s(j), :)');
        %        upOnes = sib(logical(triu(ones(size(sib)), 1)));
        %        clear sib
        finalMat(j, i) = corrQuantileReturn(dataSetProbeInf(i).name, ...
                                            sib);
    end
end

everyWhereGenes.g1s = thisG1s;
everyWhereGenes.g2s = thisG2s;
everyWhereGenes.mat = finalMat;
save('~/data/general/everywhereExpressedLinks.mat', 'everyWhereGenes')

tissueCell = {dataSetProbeInf(:).tissue};
dsIndCell = strfind(tissueCell, tissue);
allDSInd = find(~cellfun(@isempty, dsIndCell));

tic
allTissues = [1:53];
space = zeros(n, 1);
tissueInd = allDSInd;
otherInd = allTissues(~ismember(allTissues, tissueInd));
tc = length(tissueInd);
oc = length(otherInd);
ec = 947;

for p = 1: n

    small = finalMat(p, tissueInd);

    big = finalMat(p, otherInd);

    total = ec + tc + oc;

    % for each step, I should count the space under the curve. 
    allseqBig = zeros(1, total);
    allseqSM = zeros(1, total);
    for i = 1:length(big)
        allseqBig(big(i)) = allseqBig(big(i)) + 1; 
    end

    for i = 1:length(small)
        allseqSM(small(i)) = allseqSM(small(i)) + 1; 
    end

    emptyUnit = 1/ec;
    ts = 0;
    tsUnit = 1/tc; % 1/tc
    nontsUnit = 1/oc;
    ss = 0;
    steps = 0;
    tsEmptySurface = 0;
    while steps < total
        currentBig = allseqBig(1000 - steps);
        currentSM = allseqSM(1000 - steps);
        if(currentBig + currentSM)  == 0 
            tsEmptySurface = tsEmptySurface + emptyUnit * ts;
            steps = steps + 1;
        else
            for j = 1:currentBig
                ss = ss + (tsEmptySurface * nontsUnit);
                steps = steps + 1;
            end
            for j = 1:currentSM
                ts = ts + tsUnit;
                steps = steps + 1;
            end
        end
    end
    space(p)= ss;
end

space


% 6. getting links expressed in all the tissues - between them,
% getting the TS links.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
clear

genesExp = zeros(1, 18494);

tissues = {'blood', 'brain', 'skeletalMuscle', 'liver', 'lung'} 

for i = 1:length(tissues)
    tissue = tissues{i};
    % load all the brain files one by one 
    load(['~/data/general/' tissue 'ExpGenes.mat']);
    genesExp = genesExp + expGenesInd;
end

genesExp = genesExp == 5;

gEC = sum(genesExp)

gEC * (gEC-1)/2

allGenesExpNet = zeros(18494, 18494);
allGenesExpNet(genesExp, genesExp)= 1;
allGenesExpNet = triu(allGenesExpNet, 1);

[g1s, g2s, trash] = find(allGenesExpNet) ;

tissue = 'lung'
load(['~/networks/continuousTSLinksResults/' tissue ...
      'SpaceMat.mat'])

bloodSpaceNet = spaceMat;
smSpaceNet = spaceMat;
liverSpaceNet = spaceMat;
lungSpaceNet = spaceMat;
brainSpaceNet = spaceMat;

% the five space:
c = 0;

is = zeros(1, 200000);
spaces = zeros(200000, 5);

tic
for i = 1: 1000000%length(g1s)

    g1 = g1s(i);
    g2 = g2s(i);
    sib=[bloodSpaceNet(g1, g2)
         lungSpaceNet(g1, g2)
         smSpaceNet(g1, g2)
         liverSpaceNet(g1, g2)
         brainSpaceNet(g1, g2)];

    if max(sib) < .075
        c
        c = c + 1;
        is(c) = i;
        spaces(c, :) = sib;
        %    return
    end
end
toc

% 4650 links found
spaces = spaces(1:4550, :);
maxSpaces = max(spaces');
hMS = hist(maxSpaces, 40)
holu = find(maxSpaces < 0.6);
size(holu)

k = 1
is((k))
g1s(is((k)))
g2s(is((k)))
spaces((k), :)'
k = k + 1;

thisG1s = g1s(is(1:c));
thisG2s = g2s(is(1:c));

tctslinkSelectWS.spaces = spaces;
tctslinkSelectWS.is = is;
tctslinkSelectWS.g1s = g1s;
tctslinkSelectWS.g2s = g2s;

% this first one is about ts links. 
save('~/data/general/theContinuousTSlinkSelectionWS.mat', 'tctslinkSelectWS')

% coexpressed links amongst all datasets
coG1s = g1s(is(1:c));
coG2s = g2s(is(1:c));

% 6.5. getting the TS links for each tissue, the links under the
% three categoreis
clear

genesExp = zeros(1, 18494);

tissues = {'blood', 'brain', 'skeletalMuscle', 'liver', 'lung'} 

for i = 1:length(tissues)
    tissue = tissues{i};
    % load all the brain files one by one 
    load(['~/data/general/' tissue 'ExpGenes.mat']);
    genesExp = genesExp + expGenesInd;
end

genesExp = genesExp == 5;

gEC = sum(genesExp)

gEC * (gEC-1)/2

allGenesExpNet = zeros(18494, 18494);
allGenesExpNet(genesExp, genesExp)= 1;
allGenesExpNet = triu(allGenesExpNet, 1);
tissue = 'brain'
load(['~/networks/continuousTSLinksResults/' tissue ...
      'SpaceMat.mat'])
[g1s, g2s, spaces] = find(spaceMat);

tempHist = hist(spaces, 40)

sum(spaces > 0.60)

tsG1s = g1s(spaces > 0.6);
tsG2s = g2s(spaces > 0.6);

% 7. finding final mat for a given set of g1 and g2
%%%%%%%%%%%%%%%%%

load('~/data/general/linkExprInfo/wholeGeneExpr.mat')
load('~/data/general/linkExprInfo/dataSetProbeInf.mat')

spaceMat = spaceMat .* allGenesExpNet;
[g1s, g2s, space] = find(spaceMat);

[a, b] = sort(space, 'descend');

space = space(b);
n = 1000000;
thisG1s = g1s(1:n);
thisG2s = g2s(1:n);

n = length(tsG1s);%c%1; % count of edges
finalMat = zeros(n, 53);
thisG1s = tsG1s;%coG1s;%g1 = 79
thisG2s = tsG2s;%coG2s;%g2 = 409
coG1s;
coG2s;

sib = zeros(1, n);
for i = 1:53
    i
    dsName = dataSetProbeInf(i).name;
    dsMat = wholeGeneExpr(:, dataSetProbeInf(i).sampleInd(1): ...
                          dataSetProbeInf(i).sampleInd(2));

    corrMat = corr(dsMat');
    for j = 1:n
        i;
        % dsName = dataSetProbeInf(i).name;
        % dsMat = wholeGeneExpr(:, dataSetProbeInf(i).sampleInd(1): ...
        %                       dataSetProbeInf(i).sampleInd(2));
        %dsMat = dsMat(expGenesInd, :);
        %sib = corr(dsMat(1:100, :)');

        sib(j) = corrMat(thisG1s(j),thisG2s(j));
        %        upOnes = sib(logical(triu(ones(size(sib)), 1)));
        %        clear sib
    end
    tic
     finalMat(:, i) = corrQuantileReturn(dataSetProbeInf(i).name, ...
                                        sib);
toc
    %finalMat(:, i) = sib;
end

% this part is for when we have the finalMat 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% have finalMat from somewhere
% know the total count of datasets 
% know the index of tissues

% I am getting this finalMat from pureTSlinkStudy.m, lines 380+
finalMat = rankedCorrMatStr.mat;
size(finalMat)

tic
allTissues = [1:15];
space = zeros(n, 1);
tissueInd = [11:15];
otherInd = [1:10];
tc = length(tissueInd);
oc = length(otherInd);
ec = 985;

n = size(finalMat, 1);
space = zeros(n, 1);
tic
for p = 1: n

    small = finalMat(p, tissueInd);

    big = finalMat(p, otherInd);

    total = ec + tc + oc;

    % for each step, I should count the space under the curve. 
    allseqBig = zeros(1, total);
    allseqSM = zeros(1, total);
    for i = 1:length(big)
        allseqBig(big(i)) = allseqBig(big(i)) + 1; 
    end

    for i = 1:length(small)
        allseqSM(small(i)) = allseqSM(small(i)) + 1; 
    end

    emptyUnit = 1/ec;
    ts = 0;
    tsUnit = 1/tc; % 1/tc
    nontsUnit = 1/oc;
    ss = 0;
    steps = 0;
    tsEmptySurface = 0;
    while steps < total
        currentBig = allseqBig(1000 - steps);
        currentSM = allseqSM(1000 - steps);
        if(currentBig + currentSM)  == 0 
            tsEmptySurface = tsEmptySurface + emptyUnit * ts;
            steps = steps + 1;
        else
            for j = 1:currentBig
                ss = ss + (tsEmptySurface * nontsUnit);
                steps = steps + 1;
            end
            for j = 1:currentSM
                ts = ts + tsUnit;
                steps = steps + 1;
            end
        end
    end
    space(p)= ss;
end
toc

affyBrainSpace = space;
affyCommonSpace = space;

bsHist = hist(affyBrainSpace, [0:.01:1])/length(affyBrainSpace);
csHist = hist(affyCommonSpace, [0:.01:1])/length(affyCommonSpace);

plot([0:.01:1],bsHist, '*', 'color', 'r')
hold on
plot([0:.01:1],csHist, '*', 'color', 'k')

pVals = zeros(size(affyBrainSpace));
for i = 1:length(affyBrainSpace)
    pVals(i) = sum(affyCommonSpace > affyBrainSpace(i)) / length(affyCommonSpace); % 500 is where the p-value of null
                                                        % gets smaller than .05
end

hist(pVals, 40)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

cbrain = nchoosek(11:15, 3);
cblood = nchoosek(1:10, 3);

load('~/data/general/commonLinksRankedInIllData.mat')
bigAffyCommonMat = rankedCorrMatStr.mat;

load('~/data/general/brainTSlinksRankedInIllData.mat')
bigAffyBrainMat = rankedCorrMatStr.mat;

tissueInd = 4:6
otherInd = 1:3

trialCounts = length(cbrain) * length(cblood);
chooseTracker = zeros(trialCounts, 2);
cTCounter = 0;
allPVals = zeros(length(bigAffyBrainMat), trialCounts);
for ii = 1:length(cbrain)
    ii
    for jj = 1:length(cblood)
        
        % >>> getting the common variables ready 
        cTCounter = cTCounter + 1;
        chooseTracker(cTCounter, 1) = ii;
        chooseTracker(cTCounter, 2) = jj;
        jj
        % get finalMats
        cTCounter
        finalCommonMat = [bigAffyCommonMat(:, cblood(jj, :)), ...
                          bigAffyCommonMat(:, cbrain(ii, :))];
        
        finalBrainMat = [bigAffyBrainMat(:, cblood(jj, :)), ...
                         bigAffyBrainMat(:, cbrain(ii, :))];
        
        tc = length(tissueInd);
        oc = length(otherInd);
        ec = 994;
        
        % <<<<< got it. now the two lines
        
        % get space for the bs
        
        brainSpace = zeros(1, length(finalBrainMat));

        for p = 1:length(finalBrainMat)
            small = finalBrainMat(p, tissueInd);

            big = finalBrainMat(p, otherInd);

            total = ec + tc + oc;

            % for each step, I should count the space under the curve. 
            allseqBig = zeros(1, total);
            allseqSM = zeros(1, total);
            for i = 1:length(big)
                allseqBig(big(i)) = allseqBig(big(i)) + 1; 
            end

            for i = 1:length(small)
                allseqSM(small(i)) = allseqSM(small(i)) + 1; 
            end

            emptyUnit = 1/ec;
            ts = 0;
            tsUnit = 1/tc; % 1/tc
            nontsUnit = 1/oc;
            ss = 0;
            steps = 0;
            tsEmptySurface = 0;
            while steps < total
                currentBig = allseqBig(1000 - steps);
                currentSM = allseqSM(1000 - steps);
                if(currentBig + currentSM)  == 0 
                    tsEmptySurface = tsEmptySurface + emptyUnit * ts;
                    steps = steps + 1;
                else
                    for j = 1:currentBig
                        ss = ss + (tsEmptySurface * nontsUnit);
                        steps = steps + 1;
                    end
                    for j = 1:currentSM
                        ts = ts + tsUnit;
                        steps = steps + 1;
                    end
                end
            end
            brainSpace(p)= ss;
        end
        
        % get space for the cs 
        commonSpace = zeros(1, length(finalCommonMat));
        
        tc = length(tissueInd);
        oc = length(otherInd);
        ec = 994;

        for p = 1:length(finalCommonMat)
            small = finalCommonMat(p, tissueInd);

            big = finalCommonMat(p, otherInd);

            total = ec + tc + oc;

            % for each step, I should count the space under the curve. 
            allseqBig = zeros(1, total);
            allseqSM = zeros(1, total);
            for i = 1:length(big)
                allseqBig(big(i)) = allseqBig(big(i)) + 1; 
            end

            for i = 1:length(small)
                allseqSM(small(i)) = allseqSM(small(i)) + 1; 
            end

            emptyUnit = 1/ec;
            ts = 0;
            tsUnit = 1/tc; % 1/tc
            nontsUnit = 1/oc;
            ss = 0;
            steps = 0;
            tsEmptySurface = 0;
            while steps < total
                currentBig = allseqBig(1000 - steps);
                currentSM = allseqSM(1000 - steps);
                if(currentBig + currentSM)  == 0 
                    tsEmptySurface = tsEmptySurface + emptyUnit * ts;
                    steps = steps + 1;
                else
                    for j = 1:currentBig
                        ss = ss + (tsEmptySurface * nontsUnit);
                        steps = steps + 1;
                    end
                    for j = 1:currentSM
                        ts = ts + tsUnit;
                        steps = steps + 1;
                    end
                end
            end
            commonSpace(p)= ss;
        end

        % get the pvalues 
        
        pVals = zeros(size(brainSpace));
        for i = 1:length(brainSpace)
            pVals(i) = sum(commonSpace > brainSpace(i)) / length(commonSpace);
        end
        % save the pval and the 
        allPVals(:, cTCounter) = pVals;
    end
end

% 77. Doing a function for a whole tissue 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
tissue = 'blood'

load(['~/data/general/tissueExpGenes/' tissue 'ExpGenes0.8.mat']);
c = sum(expGenesInd);
finalE = c * (c-1) /2
s = 1
n = 1e7
e = 1
while e < finalE
    %contiunousTSLinksFunction_negAdjusted_ExpAdjusted(tissue, s, n)
    %    continuousTSLinksFunction_WilcoxonRankValue(tissue, s, n)    
    continuousTSLinksFunction_NA_qnorms(tissue, s, n)
    s = s + n
    e = s + n
end
n = finalE - s + 1
continuousTSLinksFunction_NA_qnorms(tissue, s, n)
%continuousTSLinksFunction_WilcoxonRankValue(tissue, s, n)    
%contiunousTSLinksFunction_negAdjusted_ExpAdjusted(tissue, s, n)

%% Part 2:

% 8. generating random list of datasets for testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first random set
% modify it as you need different numbers
clear

tissues = {'brain', 'blood', 'liver', 'lung', 'skeletalMuscle'}
load('~/data/general/linkExprInfo/dataSetProbeInf.mat')
tissueCell = {dataSetProbeInf(:).tissue};

for randSet = 1:30
    dsList = cell(1, 10);
    c = 1;
    for t = 1:length(tissues)
        tissue = tissues{t}
        dsIndCell = strfind(tissueCell, tissue);
        allDSInd = find(~cellfun(@isempty, dsIndCell));

        selected = datasample([1:length(allDSInd)], 3, 'Replace', false);
        dsList{c} = dataSetProbeInf(selected(1) + min(allDSInd) - 1).name;
        dsList{c + 1} = dataSetProbeInf(selected(2) + min(allDSInd) - ...
                                        1).name;

        % dsList{c + 2} = dataSetProbeInf(selected(3) + min(allDSInd) - ...
        %                                 1).name;
        dsListT{c} = dataSetProbeInf(selected(1) + min(allDSInd) - 1).tissue;
        dsListT{c + 1} = dataSetProbeInf(selected(2) + min(allDSInd) - ...
                                         1).tissue;
        % dsListT{c + 2} = dataSetProbeInf(selected(3) + min(allDSInd) - ...
        %                                  1).tissue;

        c = c + 2;
    end
    
    % getting 7 samples from the list of 10 randomly selected datasets
    % othInd = [datasample([1:2:10] + 1, 2, 'replace', false) [1:2: ...
    %                     10]];
    % dsList = dsList(othInd);
    % save(sprintf(['~/data/general/randomDSselections/' ...
    %               '7randomDSselection%d.mat'], randSet), 'dsList');
    
    % getting 9 samples from the list of 10 randomly selected datasets
    othInd = [datasample([1:2:10] + 1, 4, 'replace', false) [1:2: ...
                        10]];
    dsList = dsList(othInd);
    save(sprintf(['~/data/general/randomDSselections/' ...
                  '9randomDSselection%d.mat'], randSet), 'dsList');
    
    % getting 12 samples from the list of 15 randomly selected datasets
    % othInd = [datasample([1:3:15] + 2, 2, 'replace', false) [1:3: ...
    %                     15] [2:3:15]];
    % dsList = dsList(othInd);
    % save(sprintf('~/data/general/randomDSselections/12randomDSselection%d.mat', randSet), 'dsList');
end

% getting the count of expressed genes
save('~/data/general/randomDSselection2.mat', 'dsList');
save('~/data/general/randomDSselection1.mat', 'dsList');

clear

% modify in the loop in three spots:
% 1. which random selection you want to pick from,
% the 10 (default) , the 15 or the 7
% 2. Give it to the function so it is printed in the file
% 3. pick the right function

% random selection
rsString = '9'
rs = 9
rArray = [1, 3]
for i = 1:length(rArray)
    r = rArray(i)
    randSet = r;
    load(sprintf(['~/data/general/randomDSselections/%srandomDSselection%d'... 
                  '.mat'], rsString, randSet))
    load('~/data/general/tissueExpGenes/allDSExprInfo.mat')
    dataSetList = dsList;

    % getting the expressed genes
    exprMat = zeros(18494, length(dataSetList));
    [a, dsExpr] = ismember(dataSetList, {allDSExprInfo(:).GSEID});
    for j = 1:length(dataSetList)
        exprMat(:, j) = allDSExprInfo(dsExpr(j)).exprGenes;
    end

    dsThr = ceil(.8 * length(dataSetList))
    sib = sum(exprMat');
    expGenesInd = sib >= dsThr;
    sum(expGenesInd)

    count = (sum(expGenesInd) - 1) * sum(expGenesInd)/2

    s =  1
    while (count - s) >= 10e6
        s

        continuousTSLinksFunction_randomRep_NA_qnorms(dataSetList, ...
                                                      s, 10e6, r, rs)
        % continuousTSLinksFunction_randomRep_NA_ExpAdj(dataSetList, s, 10e6, ...
        %                                               r, rs)

        % continuousTSLinksFunction_randomRep(dataSetList, s, 10e6, ...
        %                                     randSet, rs)
        % continuousTSLinksFunction_randomRep_WRV(dataSetList, s, 10e6, ...
        %                                         randSet, rs)
        % continuousTSLinksFunction_randomRep_WRV_oneSided(dataSetList, ...
        %                                                  s, 10e6, ...
        %                                                  randSet, rs)
        s = s + 10e6
    end

    e = count - s;
    
    continuousTSLinksFunction_randomRep_NA_qnorms(dataSetList, ...
                                                  s, e, r, rs)

    % continuousTSLinksFunction_randomRep_NA_ExpAdj(dataSetList, s, e, ...
    %                                               r, rs)

    % continuousTSLinksFunction_randomRep(dataSetList, s, e, ...
    %                                     randSet, rs)
    % Continuoustslinksfunction_randomRep_WRV(dataSetList, s, e, ...
    %                                         randSet, rs)
    
    % continuousTSLinksFunction_randomRep_WRV_oneSided(dataSetList, ...
    %                                                  s, e, ...
    %                                                  randSet, rs)
end

%%% assembling them
clear

type = 'NegAdjusted'
subType = ''
added = 0;

type = 'pVals'
subType = '_oneSided'
subType = ''
added = 1

rCount = '15'
rCount = ''


% saving Term:
% NOTE: change the SPACE. Also, NEGADJUSTED and PVALS in the file
% names
randSets = [1:10]
rCount = 10
rString = ''
added = 0;
for r = 1:length(randSets)

    randSet = randSets(r);
    load(sprintf(['~/data/general/randomDSselections/%srandomDSselection%d'... 
                  '.mat'], rString, randSet))
    
    % load(sprintf(['~/data/general/randomDSselections/randomDSselection%d'... 
    %               '.mat'], randSet))
    
    load('~/data/general/tissueExpGenes/allDSExprInfo.mat')
    dataSetList = dsList

    % getting the expressed genes
    exprMat = zeros(18494, length(dataSetList));
    [a, dsExpr] = ismember(dataSetList, {allDSExprInfo(:).GSEID});
    for i = 1:length(dataSetList)
        exprMat(:, i) = allDSExprInfo(dsExpr(i)).exprGenes;
    end

    dsThr = ceil(.8 * length(dataSetList))
    sib = sum(exprMat');
    expGenesInd = sib >= dsThr;
    sum(expGenesInd)
    dsSetGeneCount = sum(expGenesInd)
    smallMat = zeros(dsSetGeneCount, dsSetGeneCount);

    % fill in the small tissue matrix 

    % getting the files to assemble
    % files = dir(sprintf(['~/networks/continuousTSLinksResults/negAdjusted_ExpAdjusted/random%s/%s%s_random*_%d_%s' ...
    %                     '.mat'], rCount, rCount, subType,
    %                     randSet, type))

    files = dir(sprintf(['~/networks/continuousTSLinksResults/' ...
                        'qnorms/random%d/%d_random*_%d_NA_qnorms.mat'], ...
                        rCount, rCount, randSet))

     % NOTE!!!! : change the index
    for i = 1:length(files)
        i
        kado = strsplit(files(i).name, {'m','_', '.'})
        
        s = str2num(kado{3  + added});
        e = str2num(kado{4 + added});
        
        % load(sprintf('~/networks/continuousTSLinksResults/random%s/%s' ...
        %       ,rCount, files(i).name));

        load(sprintf(['~/networks/continuousTSLinksResults/' ...
                        'qnorms/random%d/%s'], ...
                        rCount, files(i).name))
        
        % if (strcmp(type, 'pVals'))
        %     space = pVals;
        % end
         sl = length(space);
        
        % finding the beginning in the mat
        col = 2;
        for j = 1:dsSetGeneCount 
            if s > j
                s = s - j;
                col = col + 1;
            else
                row = s;
            end
        end
        
        % filling up the matrix! 
        c = 0;
        while (c < sl)
            c = c + 1;
            smallMat(row, col) = space(c);
            if (row < (col - 1))
                row = row + 1;
            else
                col = col + 1;
                row = 1;
            end
        end
    end

    % extend to the large tissue matrix 
    % gCount = 18494;
    % bigMat = zeros(gCount, gCount);
    % sib = find(expGenesInd == 1);
    % max(max(smallMat))
    
    save(sprintf(['~/networks/continuousTSLinksResults/' ...
                  'qnorms/%drandom%d.mat'], rCount, ...
                 randSet), 'smallMat')
    
    % save(sprintf('~/networks/continuousTSLinksResults/%srandom%d_smallMat_%s%s.mat', ...
    %               rCount, randSet, type, subType),'smallMat')


    % save(sprintf('~/networks/continuousTSLinksResults/%drandom%d_smallMat_%s.mat', ...
    %               rCount, randSet, type),'smallMat')
end


% UUUUUH! cause it was the upper triangle, I donno why! 
for i = 1:dsSetGeneCount
    for j = i+1: dsSetGeneCount
        bigMat(sib(i), sib(j)) = smallMat(i, j);
    end
end

spaceMat = sparse(bigMat);

book = hist(smallMat(:), [0:.02:.56]);

% phew. Making the binary network for my random sets. 
clear
randSet = 4;
load(sprintf(['~/data/general/randomDSselections/randomDSselection%d'... 
      '.mat'], randSet))
load('~/data/general/tissueExpGenes/allDSExprInfo.mat')
dataSetList = dsList;

% 9. Next TS value for the links.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for each tissue, load the space, load the binary net, get the TS
% links and get the TS value of them for other tissues. It might be
% NAN since links might not exist in other tissues. 
clear
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}

load('~/resultsAndFigures/TSlinks/TSlinks_FDR001.mat')

for k = 1:length(thrs)
    thr = thrs(k)
    TSinOthers = zeros(5, 5);
    maxHist = zeros(5, 101);
    TSHist = zeros(5, 101);
    myTSlinkCounts = zeros(1, 5);
    for i = 1:5
        i
        tissue = tissues{i}
        load(['~/networks/continuousTSLinksResults/' tissue 'SpaceMat.mat'])
        load(['~/networks/tissues/' tissue '/' ...
              'binaryNet_0.5DSCountCorr_0.8Expr.mat'])
        totalNet = spaceMat .* tempNet;
        links = totalNet > thr;
        TSValues = spaceMat(links);
        TSHist(i, :) = hist(TSValues, [0:.01:1])/length(TSValues);
        
        myTSlinkCounts(i) = sum(sum(links));

        maxTS = zeros(1, sum(sum(links>0)));

        otherInd = [1:5];
        otherInd(i) = [];
        for j = 1:4
            tissue = tissues{otherInd(j)}

            load(['~/networks/continuousTSLinksResults/' tissue 'SpaceMat.mat'])
            load(['~/networks/tissues/' tissue '/' ...
                  'binaryNet_0.5DSCountCorr_0.8Expr.mat'])
            
            tsInOthers(i, otherInd(j)) = sum(sum(tempNet .* links));
            book = (spaceMat + .01) .* links;

            [a, b, c] = find(book);

            maxTS = max(c', maxTS);
        end
        maxTSHist(i, :) = hist(maxTS, [0:.01:1])/length(maxTS);
        max(maxTS)
    end

    nextTSLinks.thr = thr;
    nextTSLinks.tissues = tissues;
    nextTSLinks.maxTSHist = maxTSHist;
    nextTSLinks.TSHist = TSHist;
    nextTSLinks.tsInOthers = tsInOthers;
    nextTSLinks.tsLinkCounts = myTSlinkCounts;
    save(sprintf('~/data/general/continuousTSLinks/nextTSLinks_Thr%0.2f.mat', ...
                 thr), 'nextTSLinks')
end

%. 10 The random networks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}

maxVector = zeros(1, 30);

temp = zeros(36000, 3);
for r = 1:30
    r
    load(sprintf('~/networks/continuousTSLinksResults/random%d_smallMat.mat', ...
                 r));
    maxVector(r) = max(max(smallMat))
    kado = smallMat(smallMat > 0);
    sib = sort(kado, 'descend');
    temp(:, r) = sib(1:36000);
    % [a, b, c]= find(smallMat);
    % kado = sort(c, 'descend');
    % subSet = kado(1:32000);
    % meHist = hist(subSet, [0:.01:1]);
    % bar(meHist)
end

save('~/networks/theTop36kTSValues_30Random.mat', 'temp');

clear
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}

maxVector = zeros(1, 30);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thr = .3
for r = 1:30
    r
    load(sprintf('~/networks/continuousTSLinksResults/random%d_smallMat.mat', ...
                 r));
    maxVector(r) = max(max(smallMat))
    kado = smallMat(smallMat > 0);
    sib = sort(kado, 'descend');
    temp(:, r) = sib(1:36000);
    % [a, b, c]= find(smallMat);
    % kado = sort(c, 'descend');
    % subSet = kado(1:32000);
    % meHist = hist(subSet, [0:.01:1]);
    % bar(meHist)
end

% copied from pureTSlilnkStudy.m

defThr = .58
for t = 1:5 
    tissue = tissues{t};
    % load(['~/networks/tissues/' tissue '/' ...
    %       'binaryNet_0.5DSCountCorr_0.8Expr.mat'])
    
    load(['~/networks/tissues/' tissue '/' ...
         'binaryNet_4QDSCounter_0.8Expr_Ind0.05.mat'])
    %tempNet/tissue Networks
    
    load(['~/networks/continuousTSLinksResults/' tissue ...
          'SpaceMat.mat']) % spaceMat

    totalTissueNet = spaceMat .* binNet;
    TSvalues = (totalTissueNet(totalTissueNet > defThr));
    allTSValues(t).tissue = tissue;
    allTSValues(t).values = TSvalues;
end

save('~/networks/TSValuesWithThr0.58_4QDS_0.8Expr_Ind0.05.mat', 'allTSValues');
save('~/networks/TSValuesWithThr0.57.mat', 'allTSValues');
load('~/networks/TSValuesWithThr0.57.mat');

h = figure

for i = 1:2
    subplot(2, 1, i)
    l = length(allTSValues(i).values)
    tissue = allTSValues(i).tissue;
    book = [temp(1:l, :) allTSValues(i).values];
    boxplot(book)
    ylabel('TS value')
    xlabel('Selection of datasets')
    title(sprintf(['Tissue Specificity of top %d links from batches ' ...
                   'of \n randomly selected datasets versus %s ' ...
                   'datasets'], l, tissue ))
    set(gca, 'XTick', [1 5 10 15 20 25 30])
    set(gca, 'XTickLabel', [1 5 10 15 20 25 30])
    grid on
end


print(h, '-dpdf', '~/data/general/tempFigures/randomPlot.pdf')
thr = .6
load(sprintf('~/data/general/continuousTSLinks/nextTSLinks_Thr%0.2f.mat', ...
             thr))

% 11. the jackknife results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% how does it even look for a tissue?

clear
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}
thrs = [.57:.01:.63]

for h = 1:length(thrs)
    thr = thrs(h);
    jkEffects = zeros(5, 5);
    load(sprintf('~/data/general/continuousTSLinks/nextTSLinks_Thr%0.2f.mat', ...
                 thr))
    diffPerMat = zeros(5, 5);
    
    for t = 1:5
        tissue = tissues{t}
        load(['~/networks/continuousTSLinksResults/' ...
              tissue 'SpaceMat']);
        load(['~/networks/tissues/' tissue '/' ...
              'binaryNet_0.5DSCountCorr_0.8Expr.mat'])
        mat = spaceMat .* tempNet;
        
        sib = sum(sum(mat > thr));
        otherInd = [1:5];
        otherInd(t) = [];
        
        for jkT = 1:4
            jkTissue = tissues{otherInd(jkT)}
            load(['~/networks/continuousTSLinksResults/' ...
                  tissue 'SpaceMat_JK_' jkTissue '.mat']);
            jkMat = spaceMat;
            jkValues = jkMat(tempNet);
            kado = sum(sum(jkValues > thr));
            jkEffects(t, otherInd(jkT)) = (kado - sib) / sib
            
        end
    end
    result(h).thr = thr;
    result(h).jkEffect = jkEffects;
end

save(sprintf('~/data/general/continuousTSLinks/JKEffects.mat', ...
             thr), 'result')


% which tissue affect the tissue specificity of other tissues the
% most? 

%%%% Running it by the tissue. n
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%e%%%%%%%%%%%%%%%%%%%

clear
tissue = 'liver'
load(['~/data/general/tissueExpGenes/' tissue 'ExpGenes0.8.mat'])
sum(expGenesInd)

count = (sum(expGenesInd) - 1) * sum(expGenesInd)/2
s =  1
while (count - s) >= 10e6
    s
    continuousTSLinksFunction_negativeAdjusted(tissue, s, 10e6)
    %continuousTSLinksFunction_WilcoxonRankValue(tissue, s, 10e6)
    s = s + 10e6;
end

e = count - s;
continuousTSLinksFunction_WilcoxonRankValue(tissue, s, e)
continuousTSLinksFunction_negativeAdjusted(tissue, s, e)


% 10. load the extended TS links and the previous ts links and
% check the correlation of the two values. 
% % 11. put together the results from psuedo tissues and do sth. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}

t = 1
tissue = tissues{t}

% load the tissue binary network
load(['~/networks/tissues/' tissue ['/binaryNet_FDR5e-' ...
                    '5_0.8Expr_Ind0.10.mat']])
bin1 = binNet;
sum(sum(bin1))
clear binNet

load(['~/networks/tissues/' tissue '/topSumNet_5_0.8Expr_Ind0.10.mat'])
max(max(binNet))
unique(binNet(:))
sum(sum(binNet > 2))
binNet = binNet > 2;
sum(sum(binNet))

load(['~/networks/continuousTSLinksResults/' tissue ...
      'SpaceMat_negAdjusted.mat']) % spaceMat
spaceMat = spaceMat + 0.001;
negAdjNet =  spaceMat .* bin1;
[naa, nab, na] = find(negAdjNet);
closeness = na;
sum(naSelected)
clear spaceMat;

% making matrix and sorted version
sparseNa = sparse(naa, nab, na, 18494, 18494);
[nas, nasInds] = sort(na, 'descend');
naas = naa(nasInds);
nabs = nab(nasInds);
sortedCloseness = 1 - nas;

load(['~/networks/continuousTSLinksResults/' tissue ...
      'SpaceMat_pVals.mat']) % spaceMat
spaceMat = spaceMat + 1e-20;
pVals =  spaceMat .* bin1; 
[pva, pvb, pv] = find(pVals);
pv = pv - 10e-12;

sparsePv = sparse(pva, pvb, pv, 18494, 18494);
[pvs, pvsInds] = sort(pv);
pvas = pva(pvsInds);
pvbs = pvb(pvsInds);

% just checking >>
holu = pvb - nab;
sum(holu)
% <<

pvNASorted = pv(nasInds);

% the rank correlation 
sib = corr(sortedCloseness, pvNASorted, 'type', 'Spearman') 
% overall rank correlation is very high .84, but I want the
% threshold milestone overlap rather than rank correlaion. 

% let's get the overlaps at different levels:
myPartitions = [1:3000:length(na)];

overlap = zeros(1, length(myPartitions)-1);

for i = 1:length(myPartitions) - 1
    e = myPartitions(i + 1);
    [a , b] = ismember(nasInds(1:e), pvsInds(1:e));
    overlap(i) = sum(a);
end

overlapPortions = overlap ./ myPartitions(2:end);
pvNASorted = pv(nasInds);
naPVSorted = na(pvsInds);

% TODO: what is missing and why, which one you like better, look at
% the FP from each of the links. 
e = myPartitions(15)
% the pv to na FP: get the minimum failing pv for the na selected
% base. 
[a, b] = ismember(pvsInds(1:e), nasInds(1:e));
tempInds = find(~a);
fpInds = pvsInds(tempInds); % these are the lowest pvalues which
                            % were not included in the lowest 
k = 1
fpInds(k)
[pva(fpInds(k)) pvb(fpInds(k))]
% It is simple from here, just get the top ones from pvsInds which
% are failed, those are the best ones which got failed

% the na to pv FN: those which were in pv, but not in na
[a, b] = ismember(nasInds(1:e), pvsInds(1:e));
tempInds = find(~a);
sum(a)/length(a)
fpInds = nasInds(tempInds); % these are the lowest pvalues which
                            % were not included in the lowest 
k = 3
fpInds(k)
[pva(fpInds(k)) pvb(fpInds(k))]

% the pv to na FN: Those which were in NA, but where not in pv.
[a, b] = ismember(nasInds(1:e), pvsInds(1:e));

% the na to pv FP

% >>>> get the conclusion that what is being missed and why. 


sparsePv = sparse(pva, pvb, pv, 18494, 18494);
[pvs, pvsInd] = sort(pv);
pvas = pva(pvsInds);
pvbs = pvb(pvsInds);

clear tind

% 11.5 the final THR 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% >> The conclusion was that I needed to get the union. For each
% method I have the distribution of the random. I will make the
% things. 

clear 
% the tissues
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}

% the NA distribution of the tissues
load('~/resultsAndFigures/TSlinks/tissues_NA_distr.mat')
tissuesNA = distr;
clear distr
load('~/resultsAndFigures/TSlinks/tissues_NA_histValues.mat')
tissuesNAHist = myHists;
clear myHists

% the pvalue dist of the tissues
load('~/resultsAndFigures/TSlinks/tissues_pVals_distr.mat')
tissuesP = distr;
clear distr

% the NA distr of the 10, 30
load('~/resultsAndFigures/TSlinks/pseudoTissues_NA_distr.mat')
pseudo10NA = distr;
clear distr

% the Pvalue dist of the 10, 30
load('~/resultsAndFigures/TSlinks/pseudoTissues_pVals_distr.mat')
pseudo10p = distr
clear distr

% the NA distr of the 15, 15
load('~/resultsAndFigures/TSlinks/15pseudoTissues_NA_distr.mat')
distr(9) = [];
pseudo15NA = distr
clear distr
load(['~/resultsAndFigures/TSlinks/' ...
      '15pseudoTissues_NA_histValues.mat'])
pseudo15NAHists = hists;

h = figure;
randHist = pseudo15NAHists(1, :)/sum(pseudo15NAHists(1, :));
plot([0:.01:.99], randHist, '-')
hold on
tissueHist = tissuesNAHist(4, :)/sum(tissuesNAHist(4, :));
plot([0:.01:.99], tissueHist, '-', 'color', 'r')

% the Pvalue dist of the 15, 15
load(['~/resultsAndFigures/TSlinks/' ...
      '15pseudoTissues_pVals_distr.mat'])
distr(9) = [];
pseudo15p = distr
clear distr

% the NA distr of the 7, 14
load('~/resultsAndFigures/TSlinks/7pseudoTissues_NA_distr.mat')
pseudo7NA = distr
clear distr

% the Pvalue dist of the 7, 14
load('~/resultsAndFigures/TSlinks/7pseudoTissues_pVals_distr.mat')
pseudo7p = distr
clear distr

h = figure
for i = 1:14
    plot(pseudo7NA(i).xi, pseudo7NA(i).fi)
    hold on
end
hold on
for i = 1:15
    plot(pseudo15NA(i).xi, pseudo15NA(i).fi, 'color', 'k')
    hold on
end

hold on
for i = 1:10
    plot(pseudo10NA(i).xi, pseudo10NA(i).fi, 'color', 'g')
    hold on
end

smna = tissuesNA(5);
lgna = tissuesNA(4);
plot(smna.xi, smna.fi, 'color', 'r')

h = figure
for i = 1:14
    plot(pseudo7p(i).xi, pseudo7p(i).fi)
    hold on
end
hold on
for i = 1:8
    plot(pseudo15p(i).xi, pseudo15p(i).fi, 'color', 'k')
    i = i + 1
    hold on
end

hold on
for i = 1:10
    plot(pseudo10p(i).xi, pseudo10p(i).fi, 'color', 'g')
    hold on
end

smna = tissuesP(5);
plot(smna.xi, smna.fi, 'color', 'r')

% load the random one and get the random one
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% >>>>> The follwoing lines did this part using the already defined
% TSlinks only. I wanted to check everything from the beginning, so
% I did the above code. 

clear
load('~/resultsAndFigures/TSlinks/TSlinksExtended.mat')
TSlinksEx = TSlinks;
load('~/resultsAndFigures/TSlinks/TSlinks.mat')

TSlinks{1}
corrs = zeros(1 ,5)
for t = 1:5
    mod = 1 -  TSlinks{t}.na;
    corrs(t) = corr(TSlinks{t}.pv, mod, 'Type', 'Spearman');
end

corrs

% get the k best links for each tissue:
k = 5  % count of links to get
t = 1  % tissue
[ana, bna] = sort(TSlinks{t}.na, 'descend');
[ap, bp] = sort(TSlinks{t}.pv);

% sorted for the NegAdjusted values
[TSlinks{t}.a(bna(1:k)), TSlinks{t}.b(bna(1:k)), TSlinks{t}.a(bp(1:k)), TSlinks{t}.b(bp(1:k))]
[ana(1:k), ap(1:k)]
TSlinks{t}.na(bp(1:k))

plot(TSlinks{t}.na, (TSlinks{t}.pv), '.')

sum(TSlinksEx{t}.na > .71)
length(TSlinksEx{t}.na)

%12. Get the hist/distr for all the neg adjusted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

% getting the list of files to read 
% select:
rCounts = [7, 15, 12, 9 10]
rCounts = 9
for i = 1:length(rCounts)
    rCount = rCounts(i)

    files = dir(sprintf(['~/networks/continuousTSLinksResults/qnorms/' ...
                        'random%d/%drandom*.mat'], rCount, rCount))
    % files = dir(sprintf(['~/networks/continuousTSLinksResults/random%d/%drandom*_NegAdjusted' ...
    %                     '.mat'], rCount, rCount))
    
    % in each file I have matrix of values
    n = length(files)
    rs = rCount
    hists = zeros(n, 100);
    distr(n, 1) = struct('xi' , zeros(1, 100), 'fi', zeros(1, 100));
    for r = 1:n

        load(sprintf(['~/networks/continuousTSLinksResults/qnorms/' ...
                        'random%d/%s'], rCount, files(r).name))

        % load(sprintf('~/networks/continuousTSLinksResults/random%d/%s', ...
        %              rCount, files(r).name))

            trueVals = smallMat ./ ((53-rs)*rs*.83);
        
        [a, b, c] = find(trueVals);
        length(c) 
        hists(r,:) = hist(c, [0:.01:.99]);
        [distr(r).fi distr(r).xi] = ksdensity(c, 'Support', [0,1])
    end

    save(sprintf('~/resultsAndFigures/TSlinks/randomDistsAndHists/%dpseudoTissues_NegAdjusted_qnorms_histValues.mat', rs), ...
         'hists')
    save(sprintf('~/resultsAndFigures/TSlinks/randomDistsAndHists/%dpseudoTissues_NegAdjusted_qnorms_distr.mat', rs), ...
         'distr')

    % save(sprintf('~/resultsAndFigures/TSlinks/randomDistsAndHists/%dpseudoTissues_NegAdjusted_ExpAdjusted_histValues.mat', rs), ...
    %      'hists')
    % save(sprintf('~/resultsAndFigures/TSlinks/randomDistsAndHists/%dpseudoTissues_NegAdjusted_ExpAdjusted_distr.mat', rs), ...
    %      'distr')

    % save(sprintf('~/resultsAndFigures/TSlinks/randomDistsAndHists/%dpseudoTissues_NegAdjusted_histValues.mat', rs), ...
    %      'hists')
    % save(sprintf('~/resultsAndFigures/TSlinks/randomDistsAndHists/%dpseudoTissues_NegAdjusted_distr.mat', rs), ...
    %      'distr')
end

h = figure
i = 1
plot(distr(i).xi, distr
(i).fi)
hold on
i = i + 1

type = 'NegAdjusted'
rCount = 10
load(sprintf('~/resultsAndFigures/TSlinks/%dpseudoTissues_%s_distr.mat', ...
             count, type))
allDist = distr;
hold on
plot(distr(1).xi, distr(1).fi)

% 12.5 get the hist/dist for all the pVals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

% getting the list of files to read 
% select:
rCounts = [7, 15, 12, 9 10]
%rCounts = [12, 9]

for i = 1:length(rCounts)

    rCount = rCounts(i)
    files = dir(sprintf(['~/networks/continuousTSLinksResults/random%d/%drandom*_pVals_oneSided' ...
                        '.mat'], rCount, rCount))

    n = length(files)
    rs = rCount
    hists = zeros(n, 100);
    modHists = zeros(n, 12);
    distr(n, 1) = struct('xi' , zeros(1, 100), 'fi', zeros(1, ...
                                                      100));
    
    for r = 1:n
        load(sprintf('~/networks/continuousTSLinksResults/random%d/%s', ...
                     rCount, files(r).name))

        trueVals = smallMat;
        [a, b, c] = find(trueVals);
        length(c) 
        %        hists(r,:) = hist(c, [0:.01:.99]);
        modHists(r, :) = hist(c, [1e-6, 5e-6, 1e-5, 5e-5, 1e-4, ...
                            5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, ...
                            5e-1]);
        %        [distr(r).fi distr(r).xi] = ksdensity(c, 'Support', [0,1])
    end
    
    save(sprintf('~/resultsAndFigures/TSlinks/randomDistsAndHists/%dpseudoTissues_pVals_OS_ModHistValues.mat', rs), ...
         'modHists')

    % save(sprintf('~/resultsAndFigures/TSlinks/randomDistsAndHists/%dpseudoTissues_pVals_OS_histValues.mat', rs), ...
    %      'hists')

    % save(sprintf('~/resultsAndFigures/TSlinks/randomDistsAndHists/%dpseudoTissues_pVals_OS_distr.mat', rs), ...
    %      'distr')
end
h = figure
i = 1
plot(distr(i).xi, distr(i).fi)
hold on
i = i + 1

load(sprintf(['~/data/general/randomDSselections/' ...
              '9randomDSselection5.mat'])) % dsList
load('~/data/general/tissueExpGenes/allDSExprInfo.mat')

[a, b] = ismember({allDSExprInfo(:).GSEID}, dsList);
b(13:21)
{allDSExprInfo(:).tissue}

% get the pValues distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rs = 15
% n = 15
% distr(n, 1) = struct('xi' , zeros(1, 100), 'fi', zeros(1, 100));
% for r = 1:n
%     r
%     load(['~/networks/continuousTSLinksResults/' files(r).name])
%     upperSingle = smallMat(logical(triu(ones(size(smallMat)), 1)));
%     [distr(r).fi distr(r).xi] = ksdensity(upperSingle,'Support', [-.001,1.001]);
% end

% plot(distr(1).xi, distr(1).fi)
% save(sprintf('~/resultsAndFigures/TSlinks/%dpseudoTissues_pVals_distr.mat', rs), ...
%      'distr')
% load('~/resultsAndFigures/TSlinks/pseudoTissues_pVals_distr.mat')
% allDist = distr;

% % get the hist for the random pvalues
% clear
% files = dir(sprintf(['~/networks/continuousTSLinksResults/*random*_pVals' ...
%                     '.mat']))

% % in each file I have matrix of values
% hists = zeros(30, 12);
% hists = zeros(30, 101);
% for r = 1:30
%     r
%     load(['~/networks/continuousTSLinksResults/' files(r).name])
%     upperSingle = smallMat(logical(triu(ones(size(smallMat)), 1)));
%     hists(r, :) = hist(upperSingle, [0:.01:1]);
%     % hists(r, :) = hist(upperSingle, [1e-6, 5e-6, 1e-5, 5e-5, 1e-4, ...
%     %                     5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1]);
% end

% save('~/resultsAndFigures/TSlinks/pseudoTissues_pVals_ModHistValues.mat', ...
%      'hists')

% save('~/resultsAndFigures/TSlinks/pseudoTissues_pVals_histValues.mat', ...
%      'hists')

% load(['~/resultsAndFigures/TSlinks/' ...
%       'pseudoTissues_pVals_ModHistValues.mat'])
% load('~/resultsAndFigures/TSlinks/pseudoTissues_pVals_histValues.mat')

% 13. Get the tissue negAdjusted and pvalues. Also, modified pvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%13.0
clear
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}

distr(5, 1) = struct('xi' , zeros(1, 100), 'fi', zeros(1, 100));
for t = 1:5
    t
    tissue = tissues{t};

    load(['~/networks/continuousTSLinksResults/negAdjusted_ExpAdjusted/' tissue '/SpaceMat_NA_EA.mat'])
    % load(['~/networks/continuousTSLinksResults/' tissue '/' tissue ...
    %       'SpaceMat_negAdjusted.mat'])
    % load the tissue binary network
    load(['~/networks/tissues/' tissue ['/binaryNet_FDR5e-' ...
                        '5_0.8Expr_Ind0.10.mat']])
    kado = binNet .* spaceMat;
    [a, b, c] = find(kado);
    [distr(t).fi distr(t).xi] = ksdensity(c, 'Support', [0,1])
end


save('~/resultsAndFigures/TSlinks/tissueATN_NA_EA_distr.mat', ...
      'distr')

save('~/resultsAndFigures/TSlinks/tissueATN_NA_distr.mat', ...
      'distr')

% 13.1: tissue Dist for NA
clear
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}

distr(5, 1) = struct('xi' , zeros(1, 100), 'fi', zeros(1, 100));
for t = 1:5
    t
    tissue = tissues{t};

    load(['~/networks/continuousTSLinksResults/negAdjusted_ExpAdjusted/' tissue '/SpaceMat_NA_EA.mat'])

    % load(['~/networks/continuousTSLinksResults/' tissue '/' tissue ...
    %       'SpaceMat_negAdjusted.mat'])
    [a, b, c] = find(spaceMat);
    [distr(t).fi distr(t).xi] = ksdensity(c, 'Support', [0,1])
end

save('~/resultsAndFigures/TSlinks/tissue_NA_EA_distr.mat', ...
      'distr')
load('~/resultsAndFigures/TSlinks/tissue_NA_distr.mat')

% plot(distr(5).xi, distr(5).fi)
% hold on
% just trying out some plots
% i = 1
% plot(allDist(i).xi, allDist(i).fi)
% hold on
% t = 5
% plot(tissueDist(t).xi, tissueDist(t).fi, 'color', 'r')

% 13.2: tissue Dist for PV
distr(5, 1) = struct('xi' , zeros(1, 100), 'fi', zeros(1, 100));
for t = 1:5
    t
    tissue = tissues{t}
    load(['~/networks/continuousTSLinksResults/' tissue ...
          'SpaceMat_pVals.mat'])
    load(['~/data/general/tissueExpGenes/' tissue ...
          'ExpGenes0.8.mat']);
    sum(sum(spaceMat > 0))
    sum(expGenesInd)* (sum(expGenesInd) - 1)/2
    % kado1 = zeros(18494, 18494);
    % kado2 = zeros(18494, 18494);
    % kado1(expGenesInd, :) = 1;
    % kado2(:, expGenesInd) = 1;
    % sib = kado1 + kado2; 
    % temp = sib == 2;
    % temp = temp .* 10e-10;
    % clear kado1 kado2 sib
    % temp = triu(temp, 1);
    % spaceMat = spaceMat + temp;
    [a, b, c] = find(spaceMat);
    [distr(t).fi distr(t).xi] = ksdensity(c,'Support', [-.001,1.001]);
end

plot(distr(1).xi, distr(1).fi)
save('~/resultsAndFigures/TSlinks/tissues_pVals_distr.mat', ...
     'distr')

% 13.3 Hist for tissues NA 
clear
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}
myHists = zeros(5, 100);
for t = 1:5
    t
    tissue = tissues{t};
    
 load(['~/networks/continuousTSLinksResults/qnorms/' tissue '/SpaceMat_NA_qnorms.mat'])
    % load(['~/networks/continuousTSLinksResults/negAdjusted_ExpAdjusted/' tissue '/SpaceMat_NA_EA.mat'])

    % load(['~/networks/continuousTSLinksResults/' tissue '/' tissue ...
    %       'SpaceMat_negAdjusted.mat'])
    [aa, bb, cc] = find(spaceMat);
    kado = hist(cc, [0:0.01:.99]);
    myHists(t, :) = kado;
end


save('~/resultsAndFigures/TSlinks/tissue_NegAdjusted_qnorms_histValues.mat', ...
     'myHists')

save('~/resultsAndFigures/TSlinks/tissue_NegAdjusted_ExpAdjusted_histValues.mat', ...
     'myHists')

save('~/resultsAndFigures/TSlinks/tissue_NegAdjusted_histValues.mat', ...
     'myHists')

% 13.35 Get the hists distr for the tisse nets. 


% 13.4 get the hist for the tissue pvalues = Modified. 
clear
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}

myHists = zeros(5, 12);
myHists = zeros(5, 101);
for t = 1:5
    t
    tissue = tissues{t};
    load(['~/networks/continuousTSLinksResults/' tissue '/' tissue 'SpaceMat_pVals_oneSided.mat'])
    [aa, bb, cc] = find(spaceMat);
    length(cc)
    myHists(t, :) = hist(cc, [1e-6, 5e-6, 1e-5, 5e-5, 1e-4, ...
                        5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, ...
                        5e-1]);
    %    myHists(t, :) = hist(cc, [0:.01:1]);
end

save('~/resultsAndFigures/TSlinks/tissue_pVals_oneSided_ModHistValues.mat', ...
     'myHists')

save('~/resultsAndFigures/TSlinks/tissue_pVals_ModHistValues.mat', ...
     'myHists')
save('~/resultsAndFigures/TSlinks/tissue_pVals_histValues.mat', ...
     'myHists')

% 13.5 : random plottings
clear
type = 'NegAdjusted'
rCount = 12
load(sprintf('~/resultsAndFigures/TSlinks/pseudoTissueAndHists/%dpseudoTissues_%s_distr.mat', ...
             rCount, type))
rDist = distr;

load(sprintf('~/resultsAndFigures/TSlinks/%dpseudoTissues_%s_histValues.mat', ...
             rCount, type))
rHists = hists;

load(['~/resultsAndFigures/TSlinks/tissueDistsAndHists/tissues_' ...
      type '_histValues.mat'])
tissueDist = distr;
load(['~/resultsAndFigures/TSlinks/tissueDistsAndHists/tissues_' ...
      type '_histValues.mat'])
tissueHists = hists;
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}
rCounts = [9, 12, 10, 15, 7]

%just trying out some plots
h = figure
hold all
for i = 1:5
    plot(rDist(i).xi, rDist(i).fi, 'color', 'k')
end

t = 2
plot(tissueDist(t).xi, tissueDist(t).fi, 'color', 'r')
hold on


% 14. getting the Thr by combining the two methods. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% todo: check the hists for each of the tissues: NA and pVals

% for each tissue, load the relevant NA and relevant Pvalues _ I
% can have the counts for the FDR. 

clear
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}
rCounts = [9, 12, 10, 15, 7]

type = 'NegAdjusted'
type = 'pVals'
type = 'pVals_oneSided'

mod = 'ModHistValues'
mod = 'histValues'
mod = 'distr'

mod = 'distr'
load(sprintf('~/resultsAndFigures/TSlinks/tissueDistsAndHists/tissues_%s_%s.mat', ...
             type, mod))
tissueDistr = distr;

mod = 'histValues'
mod = 'ModHistValues'
load(sprintf('~/resultsAndFigures/TSlinks/tissueDistsAndHists/tissues_%s_%s.mat', ...
             type, mod))
tissueHists = hists;

% These are the four groups of plots for the random and tissues,
% distr and hists. 
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% h = figure
% hold all
% for i = 1:5
%     plot(tissueDistr(i).xi, tissueDistr(i).fi)
% end

% h = figure
% hold all
% for i = 1:length(rDists)
%     plot(rDists(i).xi, rDists(i).fi, 'color', 'k')
% end

% h = figure
% hold all
% for i = 1:5
%     plot([0:.01:.99], tissueHists(i, :))
% end

% h = figure
% hold all
% for i = 1:size(rHists, 1)
%     plot([0:.01:.99], rHists(i, :))
% end
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% 15. get the average values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is for the Neg Adjusted FDRS:
clear
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}
rCounts = [9, 12, 10, 15, 7]

type = 'NegAdjusted_ExpAdjusted'
mod = 'histValues'

type = 'NegAdjusted_qnorms'
mod = 'histValues'

load(sprintf('~/resultsAndFigures/TSlinks/tissueDistsAndHists/tissue_%s_%s.mat', ...
             type, mod))
tissueHists = myHists;
% get the FDR from hists from NegAdjusted:
values = [30:1:100]
tsDiscovery = zeros(5, length(values));
tFDR = zeros(size(tsDiscovery));
for v = 1:length(values)
    v
    value = values(v)
    for t = 1:5
        myHist = tissueHists(t, :);
        rCount = rCounts(t);
        load(sprintf('~/resultsAndFigures/TSlinks/randomDistsAndHists/%dpseudoTissues_%s_%s.mat', ...
                     rCount, type, mod))

        rHists = hists;
        rSetCount =  size(rHists, 1);

        % normalize the hists for count
        normRHists = zeros(size(rHists));
        for i = 1:rSetCount
            normRHists(i, :) = rHists(i,:)/sum(rHists(i,:));
        end    
        
        normMyHist = myHist / sum(myHist);
        tsDiscovery(t, v) = sum(normMyHist(value:end));

        randomDiscovery = zeros(1, rSetCount);
        for i = 1:rSetCount;
            randomDiscovery(i) = sum(normRHists(i, value:end));
        end
        % getting the mean FDR
        tFDR(t, v) = mean(randomDiscovery)/tsDiscovery(t, v);

    end
end

save('~/resultsAndFigures/TSlinks/NegAdjusted_ExpAdjusted_FDR.mat', ...
     'tFDR')

save('~/resultsAndFigures/TSlinks/NegAdjusted_FDR.mat', 'tFDR')
load('~/resultsAndFigures/TSlinks/NegAdjusted_FDR.mat')

%% This is for the pValus FDR
clear
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}
rCounts = [9, 12, 10, 15, 7]

type = 'pVals_oneSided'
mod = 'ModHistValues'
% This is for the Neg Adjusted FDRS:

load(sprintf('~/resultsAndFigures/TSlinks/tissueDistsAndHists/tissue_%s_%s.mat', ...
             type, mod))
tissueHists = myHists;

% get the FDR from hists from NegAdjusted:
values = [1:1:12]
tsDiscovery = zeros(5, length(values));
tFDR = zeros(size(tsDiscovery));
for v = 1:length(values)
    v
    value = values(v)
    for t = 1:5
        myHist = tissueHists(t, :);
        rCount = rCounts(t);
        load(sprintf('~/resultsAndFigures/TSlinks/randomDistsAndHists/%dpseudoTissues_%s_%s.mat', ...
                     rCount, type, mod))

        rHists = hists;
        rSetCount =  size(rHists, 1);
        % normalize the hists for count
        normRHists = zeros(size(rHists));
        for i = 1:rSetCount
            normRHists(i, :) = rHists(i,:)/sum(rHists(i,:));
        end    
        
        normMyHist = myHist / sum(myHist);
        tsDiscovery(t, v) = sum(normMyHist(1:value));

        randomDiscovery = zeros(1, rSetCount);
        for i = 1:rSetCount;
            randomDiscovery(i) = sum(normRHists(i, 1:value));
        end
        tFDR(t, v) = mean(randomDiscovery)/tsDiscovery(t, v);
    end
end
save('~/resultsAndFigures/TSlinks/pVals12Counts_FDR.mat', 'tFDR')

% Just pick a relevant FDR for each tissue based on the two. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



