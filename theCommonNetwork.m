% basically, the common network and other things. 
% here I find the common network amongst the genes that are
% expressed in 80% of the datasets for each tissue. I got these
% genes from the file theContiunousTSlinkSelection. 
% 1. get the common network from individual thr
% 2. get the common network from the overal thr for each tissue 
% ?. get the common network amongs all the datasets and with one
% thr only. 
% 3. the distribution of TS values amongst all the tissues
% 3.5 distance for the next TS of the links
% 4. the stat test for the common network and also based on the
% node degree. 
% 5. Comparison of links in different thresholds for individual
% networks
% 6. Genes shared between the networks

% 1. get the common network from individual thr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

gCount = 18494;
tsCount = 5;
tissues = {'blood' , 'brain' , 'skeletalMuscle', 'liver', 'lung'}
genesInd = zeros(gCount, tsCount);

summaryTable.tissues = tissues;
summaryTable.perc = [0.4 0.6 0.8 1];
summaryTable.TS = zeros(4,5);
summaryTable.TC = zeros(4,5);
summaryTable.AC = zeros(4,1);

% loading the genes which are not expressed in tissues
nonExpMat = zeros(gCount, 5);
for i = 1:5
    tissue = tissues{i}
    load(['~/data/general/tissueExpGenes/' tissue ...
          'NonExpGenes0.2.mat'])
    nonExpMat(:, i) = nonExpGenesInd;
end

expGenes = zeros(gCount, 5);;
%nonExpMat = zeros(gCount, 5);
for i = 1:5
    tissue = tissues{i}
    load(['~/data/general/tissueExpGenes/' tissue ...
          'ExpGenes0.8.mat'])
    expGenes(:, i) = expGenesInd;
end
expGenesInd = sum(expGenes');
expGenesInd = expGenesInd == 5;

% matrix for the genes which are not expressed in all the other
% tissues - the union of those expressed only in a tissue will be
% TS genes

negExpMat = zeros(gCount, 5);
expGenes = zeros(gCount, 5);
for j = 1:4
    p = summaryTable.perc(j);
    % first go through tissues, for loading the counts and filling
    % up the matrix

    for i = 1: tsCount
        tissue = tissues{i}
        load(sprintf(['~/data/general/tissueExpGenes/%s' ...
                      'ExpGenes%.1f.mat'], tissue,p ))
        genesInd(:, i) = expGenesInd;
        tt = nonExpMat; 
        tt(:, i) = [];
        negExpMat(:, i) = sum(tt, 2);
        book = negExpMat == 4;

    end
    % count of tissue common genes for each tissue for j'th thr
    summaryTable.TC(j, :) = sum(genesInd);
    sib = sum(genesInd, 2);
    summaryTable.AC(j) = sum(sib == 5);

    % second go through tissues, now getting the TS genes
    for i = 1:tsCount
        tissues(i)
        tsExp = sib .* genesInd(:,i);
        tsExp = tsExp == 1;
        sum(tsExp)
        tsExp = tsExp + book(:, i);
        summaryTable.TS(j, i) = sum(tsExp == 2);
    end
end

save('~/data/general/tissueExpGenes/summaryTable.mat', ...
     'summaryTable')

load('~/data/general/tissueExpGenes/summaryTable.mat')

% I want to base my thr on Ribosome and how well could my thr
% retrieve ribosome links and the list from Jesse's. One way is
% that, and the other way is to simply sum up and put the final
% thr. I want to see how well they overlap and how well do they
% perform on the same ribosome and everywhere links. parts 2 and 3
% will be about that. 
% >> This part will only asve the common net for all tissues, the
% rest of info is for comparison
% 2. get the common network from the overal thr for each tissue 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
gCount = 18494;
tsCount = 5;
tissues = {'blood' , 'brain' , 'skeletalMuscle', 'liver', 'lung'}

% common network for all of the tissues:
% for each tissue
expThr = '0.6'
aggNet = zeros(gCount, gCount);
for t = 1:tsCount
    tissue = tissues{t}

    %%%%%%%%%%%% Pick one - the top two are the old versions of
    %%%%%%%%%%%% file naming and should be fixed
    % load(['~/networks/tissues/' tissue '/' ...
    %                     'sumNetworkExpThrAllTissue' expThr ...
    %       '_0.10.mat']) % for expThr 0.8 and common genes

    % load(['~/networks/tissues/' tissue '/' tissue 'SumNetworkExpThr' ...
    %      expThr '_0.10.mat']) % for expThr 0.8 and tissueEXp genes

    % load(['~/networks/tissues/' tissue '/' tissue 'SumNetworkExpThr' ...
    %      expThr '_0.10.mat']) % for ExpThr 0.6 and tissueExp Genes

    load(['~/networks/tissues/' tissue '/' tissue 'SumNetworkExpThrAll' ...
         expThr '_0.10.mat']) % for expThr 0.6 and common genes only
    %%%%%%%%%%%%

    thr = ceil(max(max(sumNet))  * .5 )
    tempNet = sumNet >= thr;
    sum(sum(tempNet))/(gCount * (gCount -1)/2)
    book = sum(tempNet);
    sum(book > 0)
    aggNet = aggNet + tempNet;
end
sum(sum(aggNet == 5)) / (gCount * (gCount-1)/2)
aggNet = aggNet == 5;
aggNet = sparse(aggNet);

sum(sum(aggNet))

save(['~/networks/tissues/commonNetAllTS_' expThr '_.5.mat'], 'aggNet')
aggNet = aggNet';

[a, b, c]  = find(aggNet);

length(unique([a', b']))

commonNetGenes = (unique([a', b']));

load('~/data/general/ribosomeIDsGO.mat')

[ribA, ribB] = ismember(ribIDs, commonNetGenes);
cmNetRibID = ribIDs(ribA);

sum(ribA)

kado = aggNet(cmNetRibID, cmNetRibID);
kado = aggNet(cmNetRibID, :);
size(kado)
sum(sum(kado))

sum(sum(kado))

load('~/data/general/hkeNet.mat')
holu = sum(hkeNet);
hkeGeneInd = find(holu >0);

[a, b] = ismember(hkeGeneInd, commonNetGenes);
sum(a)

sib = aggNet + hkeNet';
sum(sum(sib ==2))

load('~/data/general/hkgInd.mat')
[a, b]  = ismember(hkgInd, commonNetGenes);

% get the median filtered networks 

% 2.2 the tissue networks shared between their genes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
gCount = 18494;
tsCount = 5;
tissues = {'blood' , 'brain' , 'skeletalMuscle', 'liver', 'lung'}
thrs = [7, 8 , 6, 7 9]

% common network for all of the tissues:
% for each tissue
expThr = '0.8'
for t = 3:4%1:tsCount
    tissue = tissues{t}
    load(['~/networks/tissues/' tissue '/' ...
                       tissue 'SumNetworkExpThr' expThr '_0.100.mat'])
    %    thr = ceil(max(max(sumNet))  * .5 )
    thr = thrs(t);
    tempNet = sumNet >= thr;
    sum(sum(tempNet))/(gCount * (gCount -1)/2)
    book = sum(tempNet) + sum(tempNet');
    sum(book > 0) 
    file = ['~/networks/tissues/' tissue '/' ...
                        'binaryNet_FDR5e-5_' expThr 'Expr_Ind0.10.mat'
    binNet = sparse(tempNet);
    save(file, 'binNet')
end

expThr = '0.8'
for t = 1:5
    tissue = tissues{t}
    load( ['~/networks/tissues/' tissue '/' ...
                        'binaryNet_0.5DSCountCorr_' expThr ...
            'Expr.mat'])
    gc = summaryTable.TC(3, i)
    sum(sum(tempNet)) / (gc * (gc -1))/2

end
% 3. study of the tissue networks and TS links: I am getting the
% distribution of the TS values for links between everywhere
% expressed links, 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

gCount = 18494;
tsCount = 5;
tissues = {'blood' , 'brain' , 'skeletalMuscle', 'liver', 'lung'}
load('~/networks/tissues/commonNetAllTS_.08_.5.mat') % aggNet
load('~/data/general/tissueExpGenes/allExp0.8.mat') % everywhere
                                                    % Exp Genes

totTissueNetTSHist = zeros(101, 5);
tissueComGsNetTSHist = zeros(101, 5);
tissueComNetTSHist = zeros(101, 5);
unComNetTSHist = zeros(101, 5);
counts = zeros(5, 5);
txtInfo = 'each column is a tissue'
for t = 1:5 

    tissue = tissues{t};
    load(['~/networks/tissues/' tissue '/' ...
          'binaryNet_0.5DSCountCorr_0.8Expr.mat'])
    %tempNet/tissue Networks
    
    load(['~/networks/continuousTSLinksResults/' tissue ...
          'SpaceMat.mat']) % spaceMat

    totalTissueNet = spaceMat .* tempNet;
    allExpIndNet = zeros(18494, 18494);
    allExpIndNet(allExp, allExp) = 1;
    tissueCommonGenesNet = spaceMat .* tempNet .* allExpIndNet;
    tissueCommonNet = spaceMat .* aggNet;
    unCommonNet = totalTissueNet - tissueCommonNet;
    unCommonGenesNet = totalTissueNet - tissueCommonGenesNet;

    bins = [0:0.01:1];

    % tissue specificity of the links for the genes exp in the tissue
    [a, b, c] = find(totalTissueNet);
    tempHist = hist(c, bins);
    tempHist = tempHist / length(c);
    max(c)
    totTissueNetTSHist(:, t) = tempHist;
    counts(1, t) = length(c);

    % tissue specificyt of the links which are between the common genes
    [a, b, c] = find(tissueCommonGenesNet);
    tempHist = hist(c, bins);
    tempHist = tempHist / length(c);
    max(c)
    tissueComGsNetTSHist(:, t) = tempHist;
        counts(2, t) = length(c);

    % tissue specificity of the links for the common links
    [a, b, c] = find(tissueCommonNet);
    tempHist = hist(c, bins);
    tempHist = tempHist / length(c);
    max(c)
    tissueComNetTSHist(:, t) = tempHist;
    counts(3, t) = length(c);

    % tissue specificity of the links for all but those in the
    % common links
    [a, b, c] = find(unCommonNet);
    tempHist = hist(c, bins);
    tempHist = tempHist / length(c);
    max(c)
    unComNetTSHist(:, t) = tempHist;
        counts(4, t) = length(c);

        % for the links which don't have BOTH genes in teh common
        % network, but only one of  them or none
    [a, b, c] = find(unCommonGenesNet);
    tempHist = hist(c, bins);
    tempHist = tempHist / length(c);
    max(c)
    unComGsNetTSHist(:, t) = tempHist;
        counts(5, t) = length(c);

end

TSstruct.totTissueNetTSHist = totTissueNetTSHist;
TSstruct.tissueComGsNetTSHist = tissueComGsNetTSHist;
TSstruct.tissueComNetTSHist = tissueComNetTSHist;
TSstruct.unComNetTSHist = unComNetTSHist;
TSstruct.unComGsNetTSHist = unComGsNetTSHist;
TSstruct.counts = counts;
TSstruct.tissues = tissues;
TSstruct.txtInfo = txtInfo;

save('~/networks/TissueSpecificityHistsForNetworks.mat', ...
     'TSstruct')

% 3.5 the stat for the next TS of the links
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for the links in the tissue, I get the TS value, then I get the
% TS value for those links from other tissues. Then I set a thr on
% the values and get like the distribution of the distance between
% them. 
clear

gCount = 18494;
tsCount = 5;
tissues = {'blood' , 'brain' , 'skeletalMuscle', 'liver', 'lung'}
load('~/networks/tissues/commonNetAllT_.S08_.5.mat') % aggNet
load('~/data/general/tissueExpGenes/allExp0.8.mat') % everywhere

clear                                                    % Exp Genes
thr = 0.01
for t = 1:5

    tissue = tissues{t};
    % load(['~/networks/tissues/' tissue '/' ...
    %       'binaryNet_0.5DSCountCorr_0.8Expr.mat'])
    %binNet/tissue Networks
    
    load(['~/networks/tissues/' tissue '/' ...
      'binaryNet_4QDSCounter_0.8Expr_Ind0.05.mat'])
    mainBinNet = binNet;
    clear binNet;

    load(['~/networks/continuousTSLinksResults/' tissue ...
          'SpaceMat.mat']) % spaceMat for the tissue

    totalTissueNet = spaceMat .* binNet;
    finalNet = (totalTissueNet > .3);
    myTSlinks = finalNet .* spaceMat;

    gCount = 18494;
    maxMat = zeros(gCount, gCount);
    maxTissueMat = zeros(gCount, gCount);
    book = zeros(2, gCount, gCount);

    % getting the maximum TS value for each link. 
    for i = 1:5 
        i
        if i ~= t % if not the current tissue
            i
            tissue = tissues{i};
            load(['~/networks/continuousTSLinksResults/' tissue ...
                  'SpaceMat.mat']) % spaceMat for the tissue
            
             load(['~/networks/tissues/' tissue '/' ...
                   'binaryNet_4QDSCounter_0.8Expr_Ind0.05.mat'])

            tMat = spaceMat .* binNet .* mainBinNet;
            book(1, : ,:) = tMat;
            book(2, : , :) = maxMat;
            holu = squeeze(max(book));
            diffMat = holu - maxMat;
            maxTissueMat(diffMat ~= 0) = i; 
            fmtm = reshape(maxTissueMat, gCount, gCount);
            maxMat = holu;
        end
    end
    
    nextBest = maxMat .* finalNet;
    diff = myTSlinks - nextBest;
    
    [k1, k2, k3] = find(diff);

    clear holu

    % how many of the links are presenet in any other tissues
    sum(sum(maxMat > 0)) / sum(sum(binNet))
    maxMat = sparse(maxMat);

    holu = c > thr;
    sum(holu)

    intA = a(holu);
    intB = b(holu);

    indices = [intA, intB];

    sib = zeros(sum(holu), 1);
    sib2 = zeros(sum(holu), 1);

    for i = 1:length(sib)
        sib(i) = maxMat(intA(i), intB(i));
        sib2(i)  = fmtm(intA(i), intB(i));
    end

    bins = [0:.01:1];
    kado = hist(sib, bins)
    sum(holu)

    TSlinks = [intA, intB, c(holu), sib, sib2];

    
    save(sprintf('~/networks/tissues/%s/%sTSlinksThr%f_TSVal_nextImmTSVal_4QDSCounter_0.8Expr_Ind0.05.mat', tissues{t},tissues{t}, thr), ...
         'TSlinks')

end

clear
tissue = 'blood'
load(['~/networks/tissues/' tissue '/' tissue 'TSlinksThr0.63_TSVal_nextImmTSVal.mat'])

% each link should select a fav tissue, but then, I also need to
% put a thr since I don't want them too close so they get mixed
% up. 
clear
tissues = {'blood' , 'brain' , 'skeletalMuscle', 'liver', 'lung'}
for i = 1:5
    i
    tissue = tissues{i}
    load(['~/networks/tissues/' tissue '/' tissue ...
          'TSlinksThr0.63_TSVal_nextImmTSVal.mat'])
    tsLinkStr(i).tissue = tissue;
    tsLinkStr(i).tsl = TSlinks;
    ung = unique([TSlinks(:,1)', TSlinks(:,2)']);
    tsLinkStr(i).uniqueGenes = ung;
end

% finding the genes which have different TS partners  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% getting the count of TS partners for a gene
genes = zeros(1, 18494);
for i = 1:5
    genes(tsLinkStr(i).uniqueGenes) = genes(tsLinkStr(i).uniqueGenes) ...
        + 1;
end
sum(genes > 1)
multiTSpartnerGenes = genes > 1;

% this is the net in which I have all the TS links, with the thr
% from the [tissue]TSlinksThr[thr]_TSVal_nextImmTSVal.mat
gCount = 18494
specificNet = zeros(gCount, gCount);
for i = 1:5
    i
    for j = 1:length(tsLinkStr(i).tsl)
        specificNet(tsLinkStr(i).tsl(j, 1), tsLinkStr(i).tsl(j,2)) ...
            = i;
    end
end
sum(sum(specificNet > 0)) % total TS links
specificNet = sparse(specificNet);

specificNetStr.net = specificNet;
specificNetStr.tissues = tissues;
specificNetStr.note = 'all the TS links with thr .63'

save('~/networks/SpecificNet_thr0.63.mat', 'specificNetStr')
load('~/networks/SpecificNet_thr0.63.mat')
gCount = 18494
template = zeros(gCount, gCount);
template(:, multiTSpartnerGenes) = 1;
template(multiTSpartnerGenes, :) = 1;
multiGeneSpecificNet = specificNet .* template;
mgsn = sparse(multiGeneSpecificNet);
sum(sum(mgsn>0))
mgsn(1:1000, 1:1000)

sib = sum(mgsn);
sum(sib > 0)

tt = find(multiTSpartnerGenes);
sib = zeros(1, length(tt));
for i = 1:length(tt)
    tt(i)
    mgsn(tt(i), :)
    book = find(mgsn(tt(i), :)>0)
    sib(i) = ((sum(mgsn(tt(i), book) == 1) > 0) && (sum(mgsn(tt(i),book) ...
                                                      ==2)>0) &&...
 (sum(mgsn(tt(i), book) == 5) > 0));
    tissues;
end

i = 4
tt(i)
mgsn(tt(i), :)
tissues
i = i + 1
% for a given set of genes, I want to get all the networks and links
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
holu = find(mgsn(tt(i),:) >0)
myNetGenes = [tt(i), holu] % my set of genes now


load('~/data/general/GPL570GemmaMapNEW.mat')
myGeneSymbols = gpl570.uniqueSymbols(myNetGenes) % just getting the
                                                 % list

template = sparse(zeros(gCount, gCount));
template(myNetGenes, myNetGenes) = 1;

load('~/networks/tissues/commonNetAllTS_.08_.5.mat') % aggNet
myNet = aggNet .* template;
[a, b, c] = find(myNet); 
sum(a)

load('~/networks/SpecificNet_thr0.63') % TS net
specificNetStr.net .* template

tissues = specificNetStr.tissues;

tissueNets06 = struct
tissueNet08 = struct
for t = 1:length(tissues)
    tissue = tissues{t}
    load(['~/networks/continuousTSLinksResults/' tissue ...
          'SpaceMat.mat']) % spaceMat for the tissue

    load(['~/networks/tissues/' tissue '/' ...
          'binaryNet_0.5DSCountCorr_0.8Expr.mat'])
    %tempNet/tissue Networks
    tissueNets08(t).net = tempNet .* spaceMat .* template;

    
    load(['~/networks/tissues/' tissue '/' ...
          'binaryNet_0.5DSCountCorr_0.6Expr.mat'])
    %tempNet/tissue Net
    tissueNets06(t).net = tempNet .* spaceMat .* template;
end

i = 1
tissueNets08(i).net
tissueNets06(i).net
i = i + 1

netMat = zeros(20, 10); % ID, ID2, commonNet, TSNet,
                        % tissuesNet(5).8, commonNet(5).6
c = 1
for i = 1:length(myNetGenes)
    ID1 = myNetGenes(i);
    for j = i+1: length(myNetGenes)
        ID2 = myNetGenes(j);
        netMat(c, 1) = myNetGenes(i);
        netMat(c, 2 ) = myNetGenes(j);
        netMat(c, 3) = myNet(ID1, ID2); % common net
        netMat(c, 4) = specificNetStr.net(ID1, ID2); % TS net
        for t = 1:5
            netMat(c, 4+t) = tissueNets08(t).net(ID1, ID2);
            netMat(c, 9+t) = tissueNets06(t).net(ID1, ID2);
        end
            c = c + 1;
    end
end

book = sum(netMat(:, 3:end), 2) > 0;
sum(book)

netMat = netMat(book, :);

% reading netMat
i = 1

[gpl570.uniqueSymbols(netMat(i, 1)), gpl570.uniqueSymbols(netMat(i, 2))]
netMat(i, 3:4)

[netMat(i,5:9);netMat(i, 10:end)]
tissues
[netMat(i, 1), netMat(i, 2)]
i = i + 1

% making the table for the links

% getting all the links for a given gene: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load('~/networks/tissues/commonNetAllTS_.08_.5.mat') % aggNet
    load('~/networks/SpecificNet_thr0.63.mat') % TS net
sumEdge = zeros(length(tt), 3);

IDs = randi(18494, 597, 1);
for j = 1:length(tt)
    j
      ID = tt(j)
      % ID = IDs(j)
    dummyVec = zeros(1, 18494);

    dummyVec = aggNet(ID, :) + dummyVec;
    dummyVec = aggNet(:, ID)' + dummyVec;
    sumEdge(j, 1) = sum(dummyVec > 0);

    dummyVec = specificNetStr.net(ID, :) + dummyVec;
    dummyVec = specificNetStr.net(:, ID)' + dummyVec;
    sumEdge(j, 2) = sum(dummyVec > 0) - sumEdge(j, 1);


    for i = 1:5
        %dummyVec = zeros(1, 18494);
        tissue = tissues{i};
        load(['~/networks/tissues/' tissue '/' ...
              'binaryNet_0.5DSCountCorr_0.8Expr.mat']);
        dummyVec = tempNet(ID, :) + dummyVec;
        dummyVec = tempNet(:, ID)' + dummyVec;
        %sum(dummyVec > 0) - 76
        % i = i + 1
    end
    sumEdge(j, 3) = sum(dummyVec > 0);
end

sumEdgeCopy = sumEdge;
book = dummyVec > 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% >>>>>>>>>> then compare it with Illumina


% and then, how many genes are changing partners. Those, I would
% like to see. How many of them are expressed in both tissues and
% are changing partner? 

% also, how do their partners change function. 

% lung is suspicious. Maybe I should check it with three batches of
% half/ count datasets? hum? That means re-running all this for
% lung. But I will have that anyway? 

% expression level of these genes and their partners. 

% another batch of experiment is that how many of the TS links
% shared between genes which have stable expression level amongst
% the tissues? I think for that I should first find the expression
% levels dist for the tissues. 

% what is the general message? 
% the general message is that if a link is presenet in one tissue,
% having both its genes presenet, there is %% chance you observe it
% in other tissues as well. 
% general functional prediction of those links

% another study is to find the genes which are expressed at similar
% levels here. How would I find these genes? And then, find the
% links in them. What is the model for that? 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tempInd = finalNet > 0;
allTSValues = finalNet(tempInd);
TSHist = hist(allTSValues, [0:.1:max(allTSValues)]);

thr = .5
TSLinks = finalNet > thr;
sum(sum(finalNet > thr))

book = TSLinks(allExp, allExp);
sum(sum(book))


% common network for all of the tissues: jackknife
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
gCount = 18494;
tsCount = 5;
tissues = {'blood' , 'brain' , 'skeletalMuscle', 'liver', 'lung'}

% for each tissue
aggNet = zeros(gCount, gCount);
jkT = 5;
for t = 1:tsCount
    if t ~= jkT
        tissue = tissues{t}
        load(['~/networks/tissues/' tissue '/' ...
              'sumNetworkExpThr5_0.05.mat'])
        thr = ceil(max(max(sumNet))  * .5 )
        tempNet = sumNet >= thr;
        sum(sum(tempNet))/(gCount * (gCount - 1)/2)
        book = sum(tempNet);
        sum(book > 0)
        aggNet = aggNet + tempNet;
    end
end

sum(sum(aggNet == 4)) / (gCount * (gCount-1)/2)
aggNet = aggNet == 4;
aggNet = sparse(aggNet);

[a, b, c]  = find(aggNet);
length(a)

length(unique([a', b']))

commonNetGenes = (unique([a', b']));

load('~/data/general/ribosomeIDsGO.mat')
load('~/data/general/linkExprInfo/wholeGeneExpr.mat')

[ribA, ribB] = ismember(ribIDs, commonNetGenes);
cmNetRibID = ribIDs(ribA);

ribNet = aggNet(cmNetRibID, cmNetRibID);
ribExNd = sum(ribNet);
kado = aggNet(cmNetRibID, :);
ribAllNd = sum(kado);

sum(sum(kado))

load('~/data/general/hkeNet.mat')

% how many of the hkeNet genes were expressed? 
hkeNetExp = sum(hkeNet) > 0;
temp1 = hkeNetExp + expGenesInd;

sib = aggNet + hkeNet;

load('~/data/general/hkgInd.mat')
[hkgA, hkgB] = ismember(hkgInd, commonNetGenes);


% 4. the stat test for the common network and also based on the
% node degree. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%
kado = 0
p = 0.1
n = 11
b = ceil(n/2)
for i = b:n
    kado = kado + binopdf(i, n, p)
end

counter = 0
count = (7835 * 7834)/2;

for i =2511818:count
    r = binornd(15, 0.1);
    if r >=8
        r = binornd(10, 0.1);
        if (r >=5) 
            r = binornd(9, 0.1);
            if(r >=5)
                r = binornd(11, 0.1);
                if(r >=6)
                    r = binornd(7, 0.1);
                    if(r >= 4)
                        counter = 1;
                    end
                end
            end
        end
    end
end
counter

% 5. comparison of the links for different individual thresholds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the idea is that, if I am more stringient on the coexpression
% thr, I must be less stringient on the repeats. my original thr is
% at .1 for individual and .5 for repeats, now I want to see with
% what repeat for a given tissue, I can go for smallin individual
% thr and still get the most similar network. 

clear

tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}
dsCounts = [9, 12, 10, 15, 7] % phew. I just know them. 
expThr = '0.8' % I have two of these
binThr = '0.05' % I have other thr of interest: 0.05, 0.01, 0.005

% loading the network and expressed genes
for t = 1:length(tissues)
    tissue = tissues{t}
    % sumThr = .50
    dsCount = dsCounts(t);

    load(['~/data/general/tissueExpGenes/' tissue 'ExpGenes' expThr ...
          '.mat'])
    geneCount = sum(expGenesInd)

    netFolder = '~/networks/tissues/'
    load([netFolder tissue '/' tissue 'SumNetworkExpThr' expThr '_' binThr ...
          '.mat']) 
    % thisSumThr = max(max(net10)) * sumThr
    % binNet10_5 = net10 >= thisSumThr;
    % sum(sum(binNet10_5))

    % getting the density of each edge count and the .pdf
    edgeCount = zeros(1, dsCount);
    netDensity = zeros(1, dsCount); 
    genesIn = zeros(1, dsCount); 
    genesIn
    totalEdgeCount = geneCount * (geneCount - 1) /2
    for i = 0:dsCount
        edgeCount(i + 1) = sum(sum(sumNet == i));
        netDensity(i + 1) = edgeCount(i + 1)/ totalEdgeCount;
        genesIn(i + 1) = sum((sum(sumNet == i) + sum(sumNet' == i)) > 0);
    end
    
    genesInPer = genesIn / geneCount;

    netPDF = zeros(1, dsCount); 
    for i = 1:(dsCount+1)
        netPDF(i)  = sum(netDensity(i:end));
    end
    
    
    result(t).tissue = tissue;
    result(t).edgeCount = edgeCount;
    result(t).netDensity = netDensity;
    result(t).netPDF = netPDF;
    result(t).genesIn = genesIn;
    result(t).genesInPer = genesInPer;
    
    % I will study three sumThr 30% %50 70% for making my final binary
    % networks. 
    
    % >> for each  tissue, for each singelNet binary thr, I will
    % have 3 ND distributions. That means 60 plots. I am guessing
    % many of them are already fail!
    
    % each of the sumThr give one binary network. For the binary
    % networks, I want to see the ND distribution and the number of
    % genes involved. 
end

% setting questions:
% 1. how many links do I get with my smallest binary thr? For this
% part, I just want to see this. 

% The percentage of the links I get from different cut threshold
% (it is defined by 7 quantiles) and different binary thresholds I
% have. I have four binary thresholds

h = figure
colors = [228, 26, 28;
          55, 126, 184;
          77, 175, 74;
          152, 78, 163;
          225, 127, 0
         ]/255

for t = 1:length(result)
    l = length(result(t).netPDF)
    x = result(t).netPDF(3:end)
    y = [2:1:(l-1)]
    plot(y, x,'*' ,'color', colors(t, : ...
                                                      ))
    % p1 = ceil((l-1) * .3)
    % p2 = ceil((l-1) * .5)
    % p3 = ceil((l-1) * .7)
    % line([p1, p1], [0 .1], 'Color', colors(t, :))
    % line([p2, p2], [0 .1], 'Color', colors(t, :))
    % line([p3, p3], [0 .1], 'Color', colors(t, :))
    
    hold on
end
xlim([0,16])
ylim([0, .03])
title('binThr .10')

legend(tissues)
    
result005 = result;
result01 = result;
result05 = result;
result10 = result;

results{1} = result005;
results{2} = result01;
results{3} = result05;
results{4} = result10;

res = results{2}
t = 1
l = length(result(t).netPDF)
y = result(t).netDensity(3:end)
x = [2:1:(l-1)]
% fitObj = fit(y, x, 'fitType', poly2)

h = figure
plot(x, y, '*')
hold on
esy = 1 ./ (x .^ 3)
plot(x, esy, '*', 'color', 'r')


% reference
% results{1} = result005;
% results{2} = result01;
% results{3} = result05;
% results{4} = result10;

% Here I set quantiles for cutting the sum, and am looking for
% different quantiles with similar edge count in the networks.

% the commented ones are those I tried for a kind of normalization
binaryThrList = [.005 .01 .05 .1]
for r = 1:4
    myQs = [.1 .2 .35 .5 .65 .7 .8]
    res = results{r}
    netCDFQ = zeros(5, length(myQs));
    netPDFQ = zeros(5, length(myQs));
    ts = [5, 1, 3, 2, 4]
    for t = 1:5
        sib = res(ts(t)).netPDF(2:end);
        n = length(sib);
        book = quantile(1:n, myQs);
        netCDFQ(t, :) = sib(round(book));
        sib = res(ts(t)).netDensity(2:end);
        netPDFQ(t, :) = sib(round(book));
        % normNetCDFQ_p23(t, :) = sib(round(book))*(power(n, 2/3))
        % normNetCDFQ_p1(t, :) = sib(round(book))*(power(n, 1))
        % normNetCDFQ_p32(t, :) = sib(round(book))*(power(n, 3/ ...
        %                                                 2))
        % normNetCDFQ_p2(t, :) = sib(round(book))*(power(n, 2))
    end
    finalRes(r).CDFMat = netCDFQ;
    finalRes(r).PDFMat = netPDFQ;
    finalRes(r).binThr = binaryThrList(r)

    % addpath '~/codes/MATLAB/myCodes/general/'
    % h= figure

    % subplot(3, 2, 1)
    % heatmap(log10(netPDFQ), [], [], [], 'colorbar', true, 'colormap', ...
    %         'bone')
    % title(sprintf('result%d', r))

    % subplot(3, 2, 2)
    % heatmap(log10(netCDFQ), [], [], [], 'colorbar', true, 'colormap', ...
    %         'bone')

    % subplot(3, 2, 3)
    % heatmap(log10(normNetCDFQ_p23), [], [], [], 'colorbar', true, 'colormap', ...
    %         'bone')

    % subplot(3, 2, 4)
    % heatmap(log10(normNetCDFQ_p1), [], [], [], 'colorbar', true, 'colormap', ...
    %         'bone')

    % subplot(3, 2, 5)
    % heatmap(log10(normNetCDFQ_p32), [], [], [], 'colorbar', true, 'colormap', ...
    %         'bone')

    % subplot(3, 2, 6)
    % heatmap(log10(normNetCDFQ_p2), [], [], [], 'colorbar', true, 'colormap', ...
    %         'bone')
end

% so, my conclusion is that, with binThr .01, I go for secondQ, for
% binThr .05 I go for the fourthQ, for binThr .1 I go for the
% 6thQ. But I don't like this last one. 

summResult.tissues = {'sm', 'blood', 'liver', 'brain', 'lung'}
summResult.bin10Q6 = finalRes(4).CDFMat(:, 6);
summResult.bin05Q4 = finalRes(3).CDFMat(:, 4);
summResult.bin01Q2 = finalRes(2).CDFMat(:, 2)

sib = [finalRes(4).CDFMat(:, 5)';
finalRes(3).CDFMat(:, 4)';
finalRes(2).CDFMat(:, 2)']

% I have three final networks for each tissue, get the ND dist and
% gene count for these networks and see if anything fails there!
% for the remaining, just keep the network. 

% TODO: 
% for these sets of THR, for each tissue, build the networks, get
% the ND dist, get the gene count, get the overlap and decide on
% the final Networks.
clear
tissues = {'skeletalMuscle', 'blood', 'liver', 'brain', 'lung'}
dsCounts = [7, 9, 10, 12, 15]
myQs = [.1 .2 .35 .5 .65 .7 .8]
binThrList = {'0.10', '0.05', '0.01'}
sumCutQ = [5, 4, 2]
% for each (binaryThr, sumCutQ) - only 3 cause I use them in pair, for each tissue

% TODO: you wanna save teh figure here?
expThr = '0.8'
for i = 1:3 % for the pair
        binThr = binThrList{i}
    for t = 1:5 % for the tissue
        tissue = tissues{t}

        % 1. load the sumNet
        load(['~/networks/tissues/' tissue '/' tissue 'SumNetworkExpThr0.8_' ...
              binThr '.mat'])

        % 2. get the sumThr
        dsCount = dsCounts(t)
        book = quantile(1:dsCount, myQs);
        sib = round(book - .1) % the round qunatile values
        sumThr = sib(sumCutQ(i)) % getting the sum value for this  tissue

        % 3. get the networks
        binNet = sumNet >= sumThr;
        load(['~/data/general/tissueExpGenes/' tissue 'ExpGenes' expThr ...
              '.mat'])
        geneCount = sum(expGenesInd)

        % 4. get the net density
        netDensity = sum(sum(binNet)) / (geneCount * (geneCount -1)/2)
    
        % 5. get the node drege dist, (save it). 
        nd = sum(binNet) + sum(binNet');

        % 6. get the genes involved 
        genesInNet = sum(nd >0)
        genesInNet/geneCount

        kado = nd >=1;
        h = figure
        hist(nd(kado), 100)
        
        title(sprintf('tissue: %s, DSCounter: %d, binThr: %s', tissue, ...
                      sumThr, binThr))

        % 7. save the network. 
        save(sprintf(['~/networks/tissues/%s/' ...
                      'binaryNet_%dQDSCounter_0.8Expr_Ind%s.mat'], tissue, ...
                     sumCutQ(i), binThr), 'binNet')
    end
end

% 6. genes shared between the networks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}
file = ['~/networks/tissues/' tissue '/' ...
        'binaryNet_FDR5e-5_' expThr 'Expr_Ind0.10.mat']


