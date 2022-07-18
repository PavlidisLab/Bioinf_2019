% load the functional terms for each of them 
% load the pure edges: common functionality of the two genes in one
% tissue, versus their functionality in other tissues: How many of
% them have list of functions in other tissues, and they should not
% be overlapped. 
% This code is for the last part of the results and studies three
% things:
% 1. shift of functioanlity for the confidently highly expressed
% genes in ALL tissues : subset of the 3.
% 1.5 shift of functionality for confidently highly expressed in
% more than one tissue. subset of the 3. 
% 2. shift of functionality for the partners of the pure links. 
% 3. shift of functionality for any of the genes expressed in more
% than one tisssue. 
% 3.2 using JCC distance: here is to get genes whith little SOF
% 4. get the histograms. 
% 6. I check the pureLinkMat (with shift of functionality) in
% GTEx. The code is in GTEx, part 9. 
% 7. To test the clusters
% 8. the plots for the functional identity of the genes >> those
% volcano like plots are in this part, 1, 2, 3.
% 8.35 studying the selected instances: 
% 9. examining the moderate JCC instances for the cases
% I can use almost all the above code
% 10. getting the tsn for all vs things : also, this is for JCC <
% 0.01, I should do it for JCC < .1 really and then it is done. 
% 11. that tiny bar plot of count of genes expressed in multiple
% tissues 
% 12. The violin plot of the JCC and FI# for the pure, semi-pure and
% null
% 13. getting the heatmap of XYWZ
% 14. Individual examples in the last two result parts (SOF and
% pure)
% 15. links that represent functional similarity in ATN
% 16. links that represent functional similarity in TSN
% 17. pure links that represenet functional similarity 
% 18. functional implication of different groups (networks) of
%links: I find how well functions are represented in the networks
% 19. Which functions are represented in WHICh networks
% 20. The network numbers: how the TS link portion decrease 
% 21. The node degree distribution for genes which have pure link
%     vs those who don't, between the TS links
% 22. The functional enrichment for genes which have pure links,
%      versus those who do not
% 23. TANs and TSNs node degree distributions
% 24. STD of the expression for partners of pure links between
% tissues: examples of low STD with examples of high STD. 
% 25.
% 26. 
% 27.
% 28. the functions enriched in tissues 
% 29. 
% 30. response to Bioinf Reviewer: random links for TSN and pure
% (copy of 19)

% 0. get the count of represented functions for each of the genes
% in the network, in all the tissue networks. 
load('~/data/general/GOdata_GPL570_07.mat')
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
inGOIDs = GOdata.GOID(inFInds);
geneCounts = geneCounts(inFInds);

%load('~/data/clusterings/inTerms8611.mat')
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'};

FDR = '0012'
FC = '3'

load(sprintf(['~/resultsAndFigures/TSlinks/finalTables/' ...
              'finalTable_CG13_FC%slog_FDR%s.mat'], FC, FDR))

myGenes = [1:18494];
alpha = .05
tsfCounts = zeros(length(myGenes), 5);
for g = 1: length(myGenes)
    g
    geneID = myGenes(g);
    % 70. get the partners in the two networks
    atpartners = zeros(5, 18494);
    tspartners = zeros(5, 18494);
    for t = 1 : 5
        tissue = tissues{t};
        load( ['~/networks/tissues/' tissue '/' ...
               'binaryNet_FDR5e-5_0.8Expr_Ind0.10.mat'])
        atn{t} = binNet;
        atpartners(t, :) = (binNet(geneID, :)) + (binNet(:, geneID)');
        tsnet = finalTable(t).wholeNet;
        tspartners(t, :) = (tsnet(geneID, :)) + (tsnet(:, geneID)');
    end

    % first try: just by the count of the ts vs atn links for the
    % for a given gene, get its neighbours in the three ATNs 
    % 71. get the functional enrichment
    type = 'ATN';
    if strcmp(type, 'ATN')
        partners = atpartners; 
    else
        partners = tspartners;
    end

    enrichedF = zeros(5, length(inTerms));
    enrichedfp = zeros(5, length(inTerms));
    for t = 1 : 5
        tissue = tissues{t};
        geneSet = logical(partners(t,:));

        myMatP = GOdata.matP(geneSet, inFunctions);
        myMatC = GOdata.matC(geneSet, inFunctions);
        sum(geneSet);

        % sum of the functional mats (so all of them are together) 
        myTotalMat = myMatP + myMatC; %+myMatF;

        % count of genes for each of the functions
        mySumF = sum(myTotalMat);

        % getting p-value for each function
        geneSetCount = size(myTotalMat, 1);
        %            pVals = 1 - binocdf(mySumF, geneSetCount,
        %            totalSpec);
        
        mySumF(mySumF < 5) = 0;
        pVals = 1 - hygecdf(mySumF, 18494, geneCounts, sum(geneSet));

        % h = figure;
        % hist(pVals)
        % getting the FDR thr
        myTermInds = find(mySumF > 0);
        myFCount = length(myTermInds);
        miniMat = myTotalMat(:, myTermInds);
        miniPVals = pVals(myTermInds);

        [a, b] = sort(miniPVals);
        sortedPVals = miniPVals(b);
        sortedTermIDs = myTermInds(b);
        counter = 0;

        % Count of observations
        for i = 1:length(miniPVals)
            if (sortedPVals(i) <= alpha / fCount)                
                %                if (sortedPVals(i) <= i * alpha / fCount)
                %if (sortedPVals(i) <= i * alpha / length(myTermInds)) 
                counter = counter + 1; 
            end
        end
        length(miniPVals);
        counter;

        % my functions: 
        hitTermIDs = sortedTermIDs(1:counter);
        hitTerms = inTerms(hitTermIDs);
        enrichedF(t, hitTermIDs) = 1;
        enrichedfp(t, hitTermIDs) = pVals(hitTermIDs);
        tsfCount(g, t) = counter;
    end
end

save('~/resultsAndFigures/geneEnrichments/functionCount.mat', ...
     'tsfCount')

load('~/resultsAndFigures/geneEnrichments/functionCount.mat')

% plot the distribution
for t = 1 : 5
    tissue = tissues{t};
    % get the genes with any link in the network
    inNetGenes = atnds(t, :) > 0;
    netGeneCount = sum(inNetGenes);
    sib = tsfCount(inNetGenes, t);
    h = figure;
    hist(sib, [3, 10, 15:15:max(sib)])
    title(sprintf('gene count = %d', netGeneCount))
    file = ['~/resultsAndFigures/geneEnrichments/figures/' tissue 'functionRep_Hist.pdf']
    print(h, '-dpdf', file)
end

% the cGenes
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}
load('~/resultsAndFigures/cGenes.mat')
for t = 1 : 5
    tissue = tissues{t};
    % get the genes with any link in the network
    nds = atnds(t, cGenes);
    netGeneCount = sum(nds > 0);
    sib = tsfCount(cGenes, t);
    h = figure;
    hist(sib, [3, 10, 15:15:max(sib)]);
    fccGenes(t).hist = hist(sib, [3, 10, 15:15:max(sib)])/length(sib);
    fccGenes(t).binCenter = [3, 10, 15:15:max(sib)];
    title(sprintf('count of cgenes with link = %d (out of 1799)', netGeneCount))
    file = ['~/resultsAndFigures/geneEnrichments/figures/' tissue '_1799cGenes_functionRep_Hist.pdf']
    print(h, '-dpdf', file)
end

% the NOT cGenes
for t = 1 : 5
    tissue = tissues{t};
    % get the genes with any link in the network
    inGenes = atnds(t, :) > 0;
    selectedGenes = (inGenes + ~cGenes) == 2;
    nds = atnds(t, selectedGenes);
    netGeneCount = sum(selectedGenes);
    sib = tsfCount(selectedGenes, t);
    h = figure;
    hist(sib, [3, 10, 15:15:max(sib)]);
    fcNOTcGenes(t).hist = hist(sib, [3, 10, 15:15:max(sib)]) / length(sib);
    fcNOTcGenes(t).binCenter = [3, 10, 15:15:max(sib)];
    title(sprintf('count of genes with link = %d, %s', netGeneCount, ...
                  tissue))
    file = ['~/resultsAndFigures/geneEnrichments/figures/' tissue '_NOT1799cGenes_functionRep_Hist.pdf']
    print(h, '-dpdf', file)
    % h = figure
    % plot(atnds(t, selectedGenes), tsfCount(selectedGenes, t), 'o')
    % hold on
    % plot(atnds(t, cGenes), tsfCount(cGenes, t), 'o', 'color', 'r')
    % xlim([0 max(atnds(t, cGenes))])
    % ylim([0 max(atnds(t, cGenes))])
end
d
% genes which have profile in blood and something else:
for t = 1 : 5
    tissue = tissues{t};
    % get the genes with any link in the network
    fCounts = tsfCount(shiftFlists.third, t);
    sib = ;
    h = figure;
    hist(sib, [3, 10, 15:15:max(sib)]);
    sofThird(t).hist = hist(sib, [3, 10, 15:15:max(sib)]) / sum(sib);
    sofThird(t).binCenter = [3, 10, 15:15:max(sib)];
    title(sprintf('count of genes with link = %d, %s', , ...
                  tissue))
    file = ['~/resultsAndFigures/geneEnrichments/figures/' tissue '_NOT1799cGenes_functionRep_Hist.pdf']
    print(h, '-dpdf', file)
end

for t = 1:5
    tissue = tissues{t};
    notC = fcNOTcGenes(t).hist;
    C = fccGenes(t).hist;
    [a, b] = max([length(notC), length(C)])
    if b == 1
        binCenter = fcNOTcGenes(t).binCenter;
        temp = C;
        C = zeros(size(notC));
        C(1:length(temp)) = temp;
    else
        binCenter = fccGenes(t).binCenter;
        temp = notC;
        notC = zeros(size(C));
        notC(1:length(temp)) = temp;
    end

    barMat = [notC; C];
    h = figure
    bar(barMat')
    set(gca, 'Xtick', [1:length(binCenter)], 'XTickLabel', binCenter ...
             )
    title(tissue)
    file = ['~/resultsAndFigures/geneEnrichments/figures/' tissue '_NOT1799cGenes_cGenes_functionRep_Hist.pdf']
    print(h, '-dpdf', file)

end

% What percentage of these 

% getting the data for the violin plot:
violinMat = zeros(60000, 3); % value, tissue, c or not 
counter = 1;
for t = 1:5
     % get the tsfn of the rest of the genes in the network
    inNet = ((atnds(t, :) > 0) + ~cGenes) == 2;
    l = sum(inNet);
    violinMat(counter: (counter + l - 1), 1) = tsfCount(inNet, t);
    violinMat(counter: (counter + l - 1), 2) = t; 
    counter = counter + l;
    
    % get the tsfn of the common genes
    v = tsfCount(cGenes, t);
    l = length(v);
    violinMat(counter:(counter + l - 1), 1) = tsfCount(cGenes, t);
    violinMat(counter:(counter + l - 1), 2) = t;
    violinMat(counter:(counter + l - 1), 3) = 1;
    counter = counter + l;
end

violinMat = violinMat(1:(counter-1), :);

file = ['~/resultsAndFigures/commonAndPure/' ...
            'fCounts.csv']
fid = fopen(file, 'w')
fprintf(fid, ['fc,', 'tissue,', 'type\n'])
for i = 1:length(violinMat)
    fprintf(fid, ['%d,%d,%d\n'], violinMat(i, 1), violinMat(i, 2), ...
            violinMat(i, 3));
end

% 1. get the highly expressed genes with ... 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% list of genes highly expressed in more than one tissue : 2/3,
% more than one tissue: I got it from the pureTSlinkStudy.m
save('~/resultsAndFigures/commonAndPure/GeneList_23MoreThanOneTissue.mat', ...
     'myGenes')

% load the common genes : change of their functionality through the
% tissues: common genes with significant TS pratners in one tissue
% or another : the current list is fine. 

% gpl570, cGenes, tsnds and atnds are loaded in the 18. part of
% purelinkstudy 
cgIDs = find(cGenes);
sib = tsnds ./ (atnds + 1);
cGenesNdsR = sib(:, cGenes);

selectedGenes = cGenesNdsR > .2;

hitList = sum(selectedGenes > 0) > 0;
hitIDs = cgIDs(logical(hitList));

% load the function stuff from the functional Enrichment
% get the function of these genes in each of the tissues. 
myGenes = hitIDs;
myGenes = find(cGenes);
alpha = .05
tsfCounts = zeros(length(myGenes), 5);
tsfs = zeros(length(myGenes), 5, length(inTerms));
for g = 1: length(myGenes)
    g
    geneID = myGenes(g);
    % 70. get the partners in the two networks
    atpartners = zeros(5, 18494);
    tspartners = zeros(5, 18494);
    for t = 1 : 5
        tissue = tissues{t};
        load( ['~/networks/tissues/' tissue '/' ...
               'binaryNet_FDR5e-5_0.8Expr_Ind0.10.mat'])
        atn{t} = binNet;
        atpartners(t, :) = (binNet(geneID, :)) + (binNet(:, geneID)');
        tsnet = finalTable(t).wholeNet;
        tspartners(t, :) = (tsnet(geneID, :)) + (tsnet(:, geneID)');
    end

    % first try: just by the count of the ts vs atn links for the
    % for a given gene, get its neighbours in the three ATNs 
    % 71. get the functional enrichment
    type = 'ATN';
    if strcmp(type, 'ATN')
        partners = atpartners; 
    else
        partners = tspartners;
    end

    enrichedF = zeros(5, length(inTerms));
    enrichedfp = zeros(5, length(inTerms));
    for t = 1 : 5
        tissue = tissues{t};
        geneSet = logical(partners(t,:));

        myMatP = GOdata.matP(geneSet, inFunctions);
        myMatC = GOdata.matC(geneSet, inFunctions);
        sum(geneSet);

        % sum of the functional mats (so all of them are together) 
        myTotalMat = myMatP + myMatC; %+myMatF;

        % count of genes for each of the functions
        mySumF = sum(myTotalMat);

        % getting p-value for each function
        geneSetCount = size(myTotalMat, 1);
        %            pVals = 1 - binocdf(mySumF, geneSetCount,
        %            totalSpec);p
        
        mySumF(mySumF < 5) = 0;
        pVals = 1 - hygecdf(mySumF, 18494, geneCounts, sum(geneSet));

        % h = figure;
        % hist(pVals)
        % getting the FDR thr
        myTermInds = find(mySumF > 0);
        myFCount = length(myTermInds);
        miniMat = myTotalMat(:, myTermInds);
        miniPVals = pVals(myTermInds);

        [a, b] = sort(miniPVals);
        sortedPVals = miniPVals(b);
        sortedTermIDs = myTermInds(b);
        counter = 0;

        % Count of observations
        for i = 1:length(miniPVals)
            if (sortedPVals(i) <= alpha / fCount)                
                %                if (sortedPVals(i) <= i * alpha / fCount)
                %if (sortedPVals(i) <= i * alpha / length(myTermInds)) 
                counter = counter + 1; 
            end
        end
        length(miniPVals);
        counter;

        % my functions: 
        hitTermIDs = sortedTermIDs(1:counter);
        hitTerms = inTerms(hitTermIDs);
        enrichedF(t, hitTermIDs) = 1;
        enrichedfp(t, hitTermIDs) = pVals(hitTermIDs);
        tsfCount(g, t) = counter;
        tsfs(g, t, hitTermIDs) = 1;
    end
end

save('~/resultsAndFigures/commonAndPure/cGenesFunctions_1799.mat', 'tsfs', ...
     '-v7.3');

load('~/resultsAndFigures/commonAndPure/cGenesFunctions_1799.mat')
% finding genes with 10f in more than one tissue. 
highf = tsfCount > 10;
highFCount = sum(highf, 2);
IDlevel3 = find(highFCount > 2);
myGenes3 = myGenes(IDlevel3); % real IDs 
finalTsFs = tsfs(IDlevel3, :, :);
% I have 274 of such genes: looking for exterem shift of functions
% get the mutually excluded functions for each of these genes : I
% don't need ME, I just need it if they are different in even one
% group. 

g = 7
myMat = reshape(finalTsFs(g, :, :), [5, 8611]);
overlap = myMat * myMat'
heatmap(overlap, [],  [],  [], 'Colorbar', true, 'Colormap', ...
         'summer')
myGenes3(g)
gpl570.uniqueSymbols(myGenes3(g))

% I want to select the genes which have different functions in
% different tissues, and have enough functions in each of the
% tissus. 
shiftFunctionGenes = zeros(1, size(tsfs,1));
for g = 1:size(tsfs, 1)
    g
    myMat = reshape(tsfs(g, :, :), [5, 8611]);
    overlap = myMat * myMat'; % overlap of functions in the tissues
    d = diag(overlap);
    selectedT = find(d > 10); % finding the tissues which have more
                              % functions for this gene
    hitMat = zeros(5,5);
    for i = 1:5
        for j = i:5
            hitMat(i,j) = overlap(i, j) /(overlap(i,i) + .5);
            hitMat(j,i) = overlap(j, i)/(overlap(j,j) + .5);
            if ((hitMat(i,j) < .05) + (hitMat(j, i) < .05) + (d(i) > ...
                10) + (d(j) > 10) == 4)
                shiftFunctionGenes(g) = 1;
            end
        end
    end
end
sum(shiftFunctionGenes)
gs = find(shiftFunctionGenes);

firstList = find(cGenes);
firstList = firstList(logical(shiftFunctionGenes));

% print the functions for those genes 
saveFolder = '~/resultsAndFigures/geneEnrichments/shiftOfFunction/the1799Genes/'
for i = 1: length(firstList)
    i
    id = firstList(i);
    functionalEnrichment_function(id, 'ATN', saveFolder)
end

% whatendbookhthesehist(book)
% give me the count of different functionality: back up examples
% with previous publications.
addpath('~/codes/MATLAB/myCodes/general/')
heatmap(tsfCount(:, [1, 3, 4 ,5]), [],  [],  [], 'Colorbar', true, 'Colormap', ...
         'summer')

% get the partners of the pure links and do the same thing for
% them: there must be an overlap

% 2. get the pure links. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}
% get the genes expressed in more than one tissue. 
expMat = zeros(18494, 5);
for t = 1:5
    tissue = tissues{t};
    load(['~/data/general/tissueExpGenes/' tissue ...
          'ExpGenes0.8.mat'])
    expMat(:, t) = expGenesInd;
end

teCount = sum(expMat');
% for those genes, get this whole exp profile. 

% making the matrix exp template for links which could be present
% in more than one tissue
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}
% get the genes expressed in more than one tissue. 
expNetTemp = zeros(18494, 18494);
for t = 1:5
    tissue = tissues{t};
    load(['~/data/general/tissueExpGenes/' tissue ...
          'ExpGenes0.8.mat'])
    expNetTemp(expGenesInd, expGenesInd) = expNetTemp(expGenesInd, expGenesInd) + 1;
end

expNetTempMTO = expNetTemp>1; % the link is expressed in more than
                              % one tissue

% get the function of these genes in each of the tissues. 
myGenes = find(teCount > 0);

load(['~/resultsAndFigures/commonAndPure/' ...
     'multiTissueExpGenesFunctions_GO07_nonComp_expInOneT.mat'])

FDR = '0012'
FC = '3'

load(sprintf(['~/resultsAndFigures/TSlinks/finalTables/' ...
              'finalTable_CG13_FC%slog_FDR%s.mat'], FC, FDR))

% >>>> making an alternative femplate for the cgnet, where I have
% the in more than one tissue
expTemplate = zeros(18494, 18494);
moreThanOneExp = teCount > 1;
expTemplate(moreThanOneExp, moreThanOneExp) = 1;
% get all the genes in the pure net.
sumNet = zeros(18494, 18494);
for t = 1: 5
    
    % links with genes expressed in more than one tissue
    tempPure = ((expNetTempMTO + finalTable(t).noTSNet) == 2) ...
        .* t;
     
    %links with genes expressed in all the tissuesn
    % tempPure = ((finalTable(t).cgNet + finalTable(t).noTSNet) == 2) ...
    %      .* t;

    sumNet = sumNet + tempPure;
end
 
book = sum(sumNet) + sum(sumNet');

% NOTE: I dont need this part as I already have the tsfs I need for
% all the genes expressed in more than one tissue. 
% load things from functionalEnrichment.m 1690. 
% for each of the genes, get their function matrix. 
% tissues = {'blood' , 'brain', 'liver', 'lung', 'skeletalMuscle'}
% myGenes = find(book>0);
% alpha = .05
% tsfCounts = zeros(length(myGenes), 5);
% tsfs = zeros(length(myGenes), 5, length(inTerms));
% for g = 1: length(myGenes)
%     geneID = myGenes(g);
%     % 70. get the partners in the two networks
%     atpartners = zeros(5, 18494);
%     tspartners = zeros(5, 18494);
%     for t = 1 : 5
%         tissue = tissues{t};
%         load( ['~/networks/tissues/' tissue '/' ...
%                'binaryNet_FDR5e-5_0.8Expr_Ind0.10.mat'])
%         atn{t} = binNet;
%         atpartners(t, :) = (binNet(geneID, :)) + (binNet(:, geneID)');
%         tsnet = finalTable(t).wholeNet;
%         tspartners(t, :) = (tsnet(geneID, :)) + (tsnet(:, geneID)');
%     end

%     % first try: just by the count of the ts vs atn links for the
%     % for a given gene, get its neighbours in the three ATNs 
%     % 71. get the functional enrichment
%     type = 'ATN'
%     if strcmp(type, 'ATN')
%         partners = atpartners; 
%     else
%         partners = tspartners;
%     end

%     enrichedF = zeros(5, length(inTerms));
%     enrichedfp = zeros(5, length(inTerms));
%     for t = 1 : 5
%         tissue = tissues{t};
%         geneSet = logical(partners(t,:));

%         myMatP = GOdata.matP(geneSet, inFunctions);
%         myMatC = GOdata.matC(geneSet, inFunctions);
%         sum(geneSet)

%         % sum of the functional mats (so all of them are together) 
%         myTotalMat = myMatP + myMatC; %+myMatF;

%         % count of genes for each of the functions
%         mySumF = sum(myTotalMat);

%         % getting p-value for each function
%         geneSetCount = size(myTotalMat, 1);
%         %            pVals = 1 - binocdf(mySumF, geneSetCount,
%         %            totalSpec);p
        
%         mySumF(mySumF < 5) = 0;
%         pVals = 1 - hygecdf(mySumF, 18494, geneCounts, sum(geneSet));

%         % h = figure;
%         % hist(pVals)
%         % getting the FDR thr
%         myTermInds = find(mySumF > 0);
%         myFCount = length(myTermInds);
%         miniMat = myTotalMat(:, myTermInds);
%         miniPVals = pVals(myTermInds);

%         [a, b] = sort(miniPVals);
%         sortedPVals = miniPVals(b);
%         sortedTermIDs = myTermInds(b);
%         counter = 0;

%         % Count of observations
%         for i = 1:length(miniPVals)
%             if (sortedPVals(i) <= alpha / fCount)                
%                 %                if (sortedPVals(i) <= i * alpha / fCount)
%                 %if (sortedPVals(i) <= i * alpha / length(myTermInds)) 
%                 counter = counter + 1; 
%             end
%         end
%         length(miniPVals)
%         counter

%         % my functions: 
%         hitTermIDs = sortedTermIDs(1:counter);
%         hitTermps = inTerms(hitTermIDs);
%         enrichedF(t, hitTermIDs) = 1;
%         enrichedfp(t, hitTermIDs) = pVals(hitTermIDs);
%         tsfCount(g, t) = counter;
%         tsfs(g, t, hitTermIDs) = 1;
%     end
% end

% save('~/resultsAndFigures/commonAndPure/cGenesFunctions.mat', 'tsfs', ...
%      '-v7.3');
% load('~/resultsAndFigures/commonAndPure/cGenesFunctions.mat')

% getting the ATND and TSND
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

% total count of pure links
[a, b, c] = find(sumNet);
% for each of these links:
%mytslGenes = unique([a, b]);
allol = zeros(length(a), 4);
allmin = zeros(length(a), 4);
shareFunction = zeros(length(a),1);
shareFunctionJK = zeros(length(a),1);
JKDotherTs = zeros(length(a), 4);
for i = 1:length(a)
    i
    g1 = a(i);
    g2 = b(i);
    tissueID = c(i);
    g1Ind = find(myGenes == g1);
    g2Ind = find(myGenes == g2);
    g1fs = reshape(tsfs(g1Ind, :, :), [5, 3507]);
    g2fs = reshape(tsfs(g2Ind, :, :), [5, 3507]);
    
    % Identification of the change of function for true pure link
    % partners: pure link partners which share a lot of functions:
    % > 20
    
    diag1 = sum(g1fs');
    diag2 = sum(g2fs');
    
    overlap1 = g1fs * g1fs';
    overlap2 = g2fs * g2fs';
    overlaptwo = g1fs * g2fs';
    diag3 = diag(overlaptwo);
    shareFunction(i) = diag3(tissueID);
    
    shareFunctionJK(i) = shareFunction(i) / (diag1(tissueID) + ...
        diag2(tissueID) - shareFunction(i) + 1);

    % diag3 = diag(overlaptwo);
    
    % nd1 = atnds(:, g1);
    % nd2 = atnds(:, g2);
    
    % g1net = zeros(5, 18494);
    % g2net = zeros(5, 18494);
    % for t =1 :5
    %     fullNet = atn{t} + atn{t}';
    %     g1net(t, :) = fullNet(g1, :);
    %     g2net(t, :) = fullNet(g2, :);
    % end
    
    % netOverlap = g1net * g2net';
    % brainF = (g1fs(2, :) + g2fs(2,:)) ==2;    
    % liverF = (g1fs(3, :) + g2fs(3,:)) ==2;
    
    % similarity = liverF + brainF;
    % sum(similarity == 2);
    
    
    % 1. these links are pure, so the genes should not have similar
    % functions in other tissues. 
    
    % 2. The links, potentially, have other functions in other
    % tissues. 
    
    % 3. do they? Do they have any measureable difference in
    % functionality in other tissues?
    
    % 4. Both of them in different tissues?
    
    % If the two two genes have a lot in a tissue, but nothing in
    % common: any of them have more than 10 but none sharing.
    
    % let's see if they share < 1% 
    
    % compare the minimum with the shared. 
    % get the diag - 
        
    % checking what is the minimum fcount for each of the genes in
    % the other tissues.
    diags = [diag1; diag2];
    diags(:, c(i)) = [];
    minDs = min(diags);
    allmin(i, :) = minDs;
    
    % getting the JK distance for the other tissues 
    diag3Rep = diag3;
    diag3Rep(tissueID) = [];
    diag1Rep = diag1;
    diag1Rep(tissueID) = [];
    diag2Rep = diag2;
    diag2Rep(tissueID) = [];
    JKDotherTs(i, :) = diag3Rep' ./ ((diag1Rep + diag2Rep - diag3Rep') + .5);
    
    % checking what is the overlap of the genes in the other tissues
    otherOverlaps = diag(overlaptwo);
    otherOverlaps(c(i)) = [];
    allol(i, :) = otherOverlaps;
end

save('~/resultsAndFigures/tempResults/otherJKD_GO07_nonComp.mat', 'JKDotherTs')
save('~/resultsAndFigures/tempResults/mainOverlapJK_GO07_nonComp.mat', 'shareFunctionJK')
save('~/resultsAndFigures/tempResults/mainOverlap_GO07_nonComp.mat', 'shareFunction')
save('~/resultsAndFigures/tempResults/allOverlap_GO07_nonComp.mat', 'allol')
save('~/resultsAndFigures/tempResults/allmin_GO07_nonComp.mat', ...
     'allmin')

% end of saving zone!
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 

clear
load('~/resultsAndFigures/tempResults/mainOverlap.mat')
load('~/resultsAndFigures/tempResults/allOverlap.mat')
load('~/resultsAndFigures/tempResults/allmin.mat')

clear
load('~/resultsAndFigures/tempResults/otherJKD_GO07_nonComp.mat')
load('~/resultsAndFigures/tempResults/mainOverlapJK_GO07_nonComp.mat')
load('~/resultsAndFigures/tempResults/mainOverlap_GO07_nonComp.mat')
load('~/resultsAndFigures/tempResults/allOverlap_GO07_nonComp.mat')
load('~/resultsAndFigures/tempResults/allmin_GO07_nonComp.mat')
load('~/resultsAndFigures/tempResults/semiPure_otherJKD_GO07_nonComp.mat')
load('~/resultsAndFigures/tempResults/semiPure_mainOverlapJK_GO07_nonComp.mat')
load('~/resultsAndFigures/tempResults/semiPure_mainOverlap_GO07_nonComp.mat')
load('~/resultsAndFigures/tempResults/semiPure_allOverlap_GO07_nonComp.mat')
load(['~/resultsAndFigures/tempResults/' ...
      'semiPure_allmin_GO07_nonComp.mat'])

% >>> getting counpuret of genes in the tissues 
f1s = zeros(length(a), 5);
f2s = zerops(length(b), 5);
for i = 1:length(a)
    ind = find(myGenes == a(i));
    f1s(i, :) = sum(tsfs(ind, :, :) ,3);
    
    ind = find(myGenes == b(i));
    f2s(i, :) = sum(tsfs(ind, :, :), 3);
end 

semiPureData.JCothers = JKDotherTs;
semiPureData.allmin = allmin;
semiPureData.allol = allol;
semiPureData.sf = shareFunction;
semiPureData.sfJC = shareFunctionJK;
semiPureData.links = [a b c f1s f2s];
save('~/resultsAndFigures/tempResults/semiPureData_allThings.mat', 'semiPureData')
load('~/resultsAndFigures/tempResults/semiPureData_allThings.mat')

load('~/resultsAndFigures/tempResults/pureData_allThings.mat')

% 2.1 get some general stuff. 
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

JKDotherTs;
allmin;
allol;
shareFunction;
shareFunctionJK;

pureData.JCothers = JKDotherTs;
pureData.allmin = allmin;
pureData.allol = allol;
pureData.sf = shareFunction;
pureData.sfJC = shareFunctionJK;
pureData.links = [a b f1s f2s];

save('~/resultsAndFigures/tempResults/pureData_allThings.mat', ...
     'pureData')
load('~/resultsAndFigures/tempResults/pureData_allThings.mat')

% >>> getting count of genes in the tissues 
f1s = zeros(size(a));
f2s = zeros(size(b));
for i = 1:length(a)
    ind = find(myGenes == a(i));
    t = c(i);
    f1s(i) = sum(tsfs(ind, t, :));
    
    ind = find(myGenes == b(i));
    t = c(i);
    f2s(i) = sum(tsfs(ind, t, :));
end 
hist(f1s)
book = min([f1s'; f2s']);
h = figure
hist(book)
myb = hist(book);
bar(log2(myb))

% getting the final plots for the pure
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pureData
% >>>>> NUMBERS:

% >> How many links have > 5 FI in the pure tissue. 
book = pureData.links(:, [3,4]);
minPureFI = min(book');
hist(minPureFI)
sum(minPureFI > 5)
pureMinFIpass = minPureFI > 5;

% >> How many links have > 5 FI in one of the other tissues. 
book = pureData.allmin;
minPass = book > 5;
sib = sum(minPass');
sum(sib > 0)
pureMinOtherPass = sib > 0;

% >> How many links have both. 
kado = pureMinOtherPass + pureMinFIpass;
sum(kado == 2)
finalLinks = kado == 2;

% >>>>>> PLOT3: Dist of minF for the pure tissue. For both pass 5
% and not pass fieve. 
h = figure
hist(minPureFI(finalLinks), 20)
title('Distribution of pure min FI')
figFolder = '~/resultsAndFigures/functionalIdentityPureLinks/figures/'
fileName = ['PureMinFICount_hist'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);

% >>>>>> PLOT1: The dist of JCC in the pure tissue. 
% we are looking at the final links. 
h = figure
hist(pureData.sfJC(finalLinks), 20)
title('Distribution of JCC for the finalHit : the 23,825')
figFolder = '~/resultsAndFigures/functionalIdentityPureLinks/figures/'
fileName = ['pureLinks_pureTissueJCC_finalHit__hist'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);

h = figure
hist(pureData.sfJC(pureMinFIpass), 20)
title('Distribution of JCC for the pureLinks passFI > 5')
figFolder = '~/resultsAndFigures/functionalIdentityPureLinks/figures/'
fileName = ['pureLinks_pureTissueJCC_pureMinFIpass__hist'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);
 
% >>>>>> PLOT2: Distribution of JCC in other links, for each linkMin
copyJCothers = pureData.JCothers;
copyJCothers(~minPass) = 1;
[minJCothers, minInd] = min(copyJCothers');
h = figure
hist(minJCothers(finalLinks), 20)
title('Distribution of JCC in the other links, for eachlinkMin')
figFolder = '~/resultsAndFigures/functionalIdentityPureLinks/figures/'
fileName = ['pureLinks_otherTissueMinJCC_finalLinks_hist'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);
 
% >>>>>> PLOT4: Distribution of the FIcount for links from plot2
minFcounts = zeros(size(minInd));
for i = 1:length(minInd)
    minFcounts(i) = pureData.allmin(i, minInd(i));
end 
h = figure
hist(minFcounts(finalLinks), 20)
title('Distribution of FIcount in the other links, for eachlinkMin')
figFolder = '~/resultsAndFigures/functionalIdentityPureLinks/figures/'
fileName = ['pureLinks_otherTissueMinJCC_FIcount_finalLinks_hist'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);

% >>>>>> PLOT5: scatter of 3 and 4

% not correlated. It is 0.03
x = minPureFI(finalLinks);
y = minFcounts(finalLinks);
plot(x, y, '.')
corr(x', y') % the correlation is about zeros. that means that
             % there is no pattern, genes just tend to have more
             % functions annotated to them in the tissue that they
             % have the pure link in. same is for the null? 

% >>>>>> PLOT6: scatter of plots 1 and 2
w = pureData.sfJC(finalLinks);
z = (minJCothers(finalLinks));
plot(z, w, '.')
corr(z, w') % They are FAR less Correlated in this case too... it
             % is 0.03

% The question was that, if the small count of f in the min tissue
% is contributing to the small JCC . The answer is not. These genes
% just tend to be less active in the other tissues. and I want to
% argue that the genes with coexpression links in multiple tissues,
% tend to have more annotation in their connected groups: take
% ribosome and protein synthesis stuff. So No. They are not
% correlated. However, they are HIGHLY correlated in the null case:
% if 

corr(x', w') % correlation of 2 and 4: if they are less similar,
             % they also have less function? this corr is
             % important. it is 0.1 : little. very, little. So I
             % think it contributes, but not so much! 

myMat = [x' y' w z'];
sib = corr(myMat, 'type', 'Pearson')
addpath('~/codes/MATLAB/myCodes/general/')
h = figure
heatmap(sib, [],  [],  [], 'Colorbar', true, 'Colormap', ...
         'summer')

pureXYWZmat = myMat;
save('~/resultsAndFigures/commonAndPure/pureXYXZmat.mat', 'pureXYWZmat')
sib = corr(myMat, 'type', 'Spearman')
addpath('~/codes/MATLAB/myCodes/general/')
h = figure
heatmap(sib, [],  [],  [], 'Colorbar', true, 'Colormap', ...
         'summer')
title('null')

h = figure
subplot(2, 2, 1)
hist(x, 20)
title('X')

subplot(2, 2, 2)
hist(y, 20)
title('Y')


subplot(2, 2, 3)
hist(w, 20)
title('W')

subplot(2, 2, 4)
hist(z, 20)
title('Z')


% end of plotting for the pure
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% >>> cases which have overlap in the tissue
hist(shareFunction)

sum(shareFunction > 5) 

hist(shareFunctionJK)

sum(shareFunctionJK > .5)

% >>> cases which have potential diverge in other tisssues: 
% (similar to above), what is the distribution of minimum and max overlap amongst the
% other tissues? 
maxol = max(allol');
hist(maxol) % >>> the overlap in other tissues is rarely big, as
            % oppose to the tissue they have a pure link in. 

% (similar to above), what is the distribution of minimum and max JCC
maxJK = max(JKDotherTs');
hist(maxJK) % >> the small overlap is also clear in the small JCC
            % distance. But, are they functional in other tissues?

% but, how many functions do they have in the other tissues?
maxallmin = max(pureData.allmin');
h = figure
hist(maxallmin)
figure
hist(shareFunction)
% >>> in the other tissues, it is VERY rare that both genes SOF to
% high count of functions

scatter(maxallmin, shareFunctionJK', '.')
corr(maxallmin', shareFunctionJK)


% 2.2 the rest (this is just to find things)
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% check the overlap too. haha. 

ol = allol(:);
mins = allmin(:);

plot(ol, mins, 'o')
xlabel('overlap')
xlim([1 400])
ylabel('minimum')

passMins = allmin > 5; % min count of function in the other tissues

%passol = (allol ./ (allmin + .05)) < .05; % their overlap in the
                                          % other tissues
passJKD = JKDotherTs < .1;

sib = (passJKD + passMins) >1;
sib = (sum(sib')) > 0; % those with soome funcion in the other
                       % tissues and
maxotherJK = max((JKDotherTs'));
hist(maxotherJK, 10)
hist(shareFunctionJK)

sib = passMins + passol;
sib = sib == 2; % those with some funcion in the other tissues and
                % smaller overlap
halva = sum(sib') > 0; 
halim = shareFunction > 5; % those who share sum function in the
                           % selected tissue
pureJKD = shareFunctionJK > .5;

% >>> ignore this part later
kad = sib(halim, :); % those who share function ... 
sum(sum(kad))  % we have 32128 incidents from the total of 80637*4
               % incidents : .0996

% get the count of possible instances for each of the links:
mlls = sparse(a, b, ones(length(a), 1), 18494, 18494);
pt = mlls .* expNetTemp;

[d e f] = find(pt);

sum(sum(pt)) - sum(sum(mlls))

% <<< ignore this part later

% getting the links with two incidents: 
holu = sum(sib');
moreThanOne = holu > 0;
hit = find((moreThanOne + halim') == 2); 

mynet = [a(hit), b(hit), c(hit)];

% OR we could go just by links and have length(Inds) / length(sib)
% : .28

Inds = find((halva + halim') == 2);

SOFPureNet = [a(Inds), b(Inds), c(Inds)];

save('~/resultsAndFigures/commonAndPure/SOFPureNet_GO07_nonComp.mat', ...
     'SOFPureNet')

load('~/resultsAndFigures/commonAndPure/SOFPureNet.mat')
load('~/resultsAndFigures/commonAndPure/SOFPureNet_GO07_nonComp.mat')

length(unique([a(shareFunction > 10) b(shareFunction > 10)]))
length(unique([a(shareFunction> 5) b(shareFunction > 5)]))

hitGenes = unique([a(Inds) b(Inds)]);
secondList = hitGenes;

% PLOTS: 
% p2.1: function shared between the pure link genes for each
% tissue. 
for t = 1 : 5
    tissue = tissues{t};
    tlinks = find(c == t);
    sib = shareFunction(tlinks);
    tGCount = length(unique([a(tlinks) b(tlinks)]));
    h = figure;
    hist(sib, [0,5, 10, 15:15:max(sib)])
    title(sprintf('gene count = %d, link count = %d', tGCount, length(sib)))
    file = ['~/resultsAndFigures/commonAndPure/figures/' tissue 'pureLinkFOverlap_Hist_GO07_nonComp.pdf']
    print(h, '-dpdf', file)
end

% for each partner link, see if the gene has different function in
% any of the other tissues. How did we get the links?

% 2.5 semi pure link stuff
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
clear
load('~/resultsAndFigures/tempResults/semiPureData_allThings.mat')
semiPureData

minPureFC = zeros(165844, 1); 
for i = 1:length(minPureFC)
    pti = semiPureData.links(i, 3);
    minPureFC(i) = min(semiPureData.links(i, (3 + pti)), semiPureData.links(i, (8 + pti)));
end

% getting the final plots for the semi pure
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% >>>>> NUMBERS:
semiPureData

% >> How many links have > 5 FI in the pure tissue. 
sum(minPureFC >  5)
pureMinFIpass = minPureFC > 5;

% >> How many links have > 5 FI in one of the other tissues. 
book = semiPureData.allmin;
minPass = book > 5;
sib = sum(minPass');
sum(sib > 0)
pureMinOtherPass = sib > 0;

% >> How many links have both. 
kado = pureMinOtherPass + pureMinFIpass';
sum(kado == 2)
finalLinks = kado == 2;

% >>>>>> PLOT3: Dist of minF for the pure tissue. For both pass 5
% and not pass five. 
h = figure
hist(minPureFC(finalLinks), 20)
title('Distribution of semi pure min FI')
figFolder = '~/resultsAndFigures/functionalIdentityPureLinks/figures/'
fileName = ['SemiPureMinFICount_hist'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);

% >>>>>> PLOT1: The dist of JCC in the pure tissue. 
% we are looking at the final links. 
h = figure
hist(semiPureData.sfJC(finalLinks), 20)
title('Distribution of JCC for the finalHit in semi pure : the 33,714')
figFolder = '~/resultsAndFigures/functionalIdentityPureLinks/figures/'
fileName = ['semiPureLinks_pureTissueJCC_finalHit__hist'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);

h = figure
hist(semiPureData.sfJC(pureMinFIpass), 20)
title('Distribution of JCC for the semiPureLinks passFI > 5')
figFolder = '~/resultsAndFigures/functionalIdentityPureLinks/figures/'
fileName = ['semiPureLinks_pureTissueJCC_pureMinFIpass__hist'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);
 
% >>>>>> PLOT2: Distribution of JCC in other links, for each linkMin
copyJCothers = semiPureData.JCothers;
copyJCothers(~minPass) = 1;
[minJCothers, minInd] = min(copyJCothers');
h = figure
hist(minJCothers(finalLinks), 20)
title('SemiPure Distribution of JCC in the other links, for eachlinkMin')
figFolder = '~/resultsAndFigures/functionalIdentityPureLinks/figures/'
fileName = ['semiPureLinks_otherTissueMinJCC_finalLinks_hist'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);
 
% >>>>>> PLOT4: Distribution of the FIcount for links from plot2
minFcounts = zeros(size(minInd));
for i = 1:length(minInd)
    minFcounts(i) = semiPureData.allmin(i, minInd(i));
end 
h = figure
hist(minFcounts(finalLinks), 20)
title('semiPure Distribution of FIcount in the other links, for eachlinkMin')
figFolder = '~/resultsAndFigures/functionalIdentityPureLinks/figures/'
fileName = ['semiPureLinks_otherTissueMinJCC_FIcount_finalLinks_hist'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);

% >>>>>> PLOT5: JCC max for the other tissues
copyJCothers = semiPureData.JCothers;
copyJCothers(~minPass) = -1;
[maxJCothers, maxInd] = max(copyJCothers');
h = figure
hist(maxJCothers(finalLinks), 20)
title('SemiPure Distribution of JCC in the other links, for eachlinkMin')
figFolder = '~/resultsAndFigures/functionalIdentityPureLinks/figures/'
fileName = ['semiPureLinks_otherTissueMinJCC_finalLinks_hist'];
print(h, '-

depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);

% >>>>>> PLOT5: scatter of plots 3 and 4: the count of min FI is
% not correlated. It is zero
x = minPureFC(finalLinks);
y = minFcounts(finalLinks);
plot(x, y, '.')
corr(x', y) % they are NOT correlated (just like pure). The higher the first one, the
             % higher other one. How is it for the pure and semi
             % pure ...

% >>>>>> PLOT6: scatter of plots 1 and 2
w = semiPureData.sfJC(finalLinks);
z = (minJCothers(finalLinks));
plot(x, y, '.')
corr(x, y') % They are FAR less Correlated in this case too... pure

% >>>>>> PLOT7: 2 and 4: is the lack of similarity in the other tissue casued by lack
% of functionality? slightly, but not completely: the corrleation
% is low 0.1, very close to pure case.
x = minFcounts(finalLinks);
y = (minJCothers(finalLinks)); % similarity JCC
plot(x, y, '.')
corr(x', y', 'type', 'Spearman') % They are FAR less
                                 % correlated. This is also like
                                 % pure 
myMat = [x y' w, z'];
semiXYZWmat = myMat;
save('~/resultsAndFigures/commonAndPure/semiXYZWmat.mat', 'semiXYZWmat')
sib = corr(myMat, 'type', 'Spearman');
addpath('~/codes/MATLAB/myCodes/general/')
h = figure
heatmap(sib,[] ,[], [], 'Colorbar', true, 'Colormap', ...
         'summer')
title('semi pure')

h = figure
subplot(2, 2, 1)
hist(x, 20)
title('X')

subplot(2, 2, 2)
hist(y, 20)
title('Y')


subplot(2, 2, 3)
hist(w, 20)
title('W')

subplot(2, 2, 4)
hist(z, 20)
title('Z')
% end of plotting for the semi pure
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

whos minPureFC
passMinPure = (minPureFC > 5);
hist(semiPureData.sfJC(passMinPure)) % exploratory 

% those who have both genes with > FI in more than one tissue : 
passMins = semiPureData.allmin > 5;
sib = sum(passMins');
passMinOthers = (sib > 0);

copyJCothers = semiPureData.JCothers;
copyJCothers(~passMins) = 1; 
[minOtherJCC, minInd] = min(copyJCothers'); % I got this, so now I
                                            % can compare it to JCC 

% getting the distribution of the min, for those who pass the two...
sib = passMinOthers + passMinPure';
sum(sib == 2)
final = (sib == 2);

hist(semiPureData.sfJC(final)) % FIRST 
hist(minOtherJCC(final)) % SECOND >>>> so, they diverge pretty
                         % much But alos, what is the min of their
                         % FCount in the other tissues...
MOJCC = minOtherJCC(final);
% min FI count in the other tissue of min JC (it should pass the 5)
for i = 1:165844
    minFICount(i) = semiPureData.allmin(i, minInd(i));
end

hist(minFICount(final), 30) 
MOFICount = minFICount(final);

scatter(MOJCC, MOFICount)

% the distribution of min JCC for links who pass minpure and pass min others,
% how many of them pass the minimum in the pure tissue, NOT the
% overlap 
 
% 3. shift of functionality for any of the genes expressed in more
% than one tissue. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
% >>>> loading and preparing data

% Tissues And final table
FC = '3'
FDR = '0012'
load(sprintf(['~/resultsAndFigures/TSlinks/finalTables/' ...
              'finalTable_CG13_FC%slog_FDR%s.mat'], FC, FDR))
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'};

% Function matrix
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
inGOIDs = GOdata.GOID(inFInds);
geneCounts = geneCounts(inFInds);
%save('~/data/clusterings/inTerms8611.mat', 'inTerms')

% function specificity:
%sumF = sum(GOdata.matF(:, inFunctions));
sumP = sum(GOdata.matP(:, inFunctions));
%sumC = sum(GOdata.matC(:, inFunctions));

%specF = sumF / gCount;
specP = sumP / gCount;
%specC = sumC / gCount;

totalSpec =  specP;% + specC;% + specF;

% FDR
alpha = .05
 
% get the genes expressed in more than one tissue. 
expMat = zeros(18494, 5);
for t = 1:5
    tissue = tissues{t};
    load(['~/data/general/tissueExpGenes/' tissue ...
          'ExpGenes0.8.mat'])
    expMat(:, t) = expGenesInd;
end

teCount = sum(expMat');
% for those genes, get this whole exp profile. 

% get the function of these genes in each of the tissues. 
myGenes = find(teCount > 0);

alpha = .05
tsfCounts = zeros(length(myGenes), 5);
tsfs = zeros(length(myGenes), 5, length(inTerms));
counter = 0;
for g = 1: length(myGenes)
    
    g
    geneID = myGenes(g);
    % 70. get the partners in the two networks
    atpartners = zeros(5, 18494);
    tspartners = zeros(5, 18494);
    for t = 1 : 5
        tissue = tissues{t};
        load( ['~/networks/tissues/' tissue '/' ...
               'binaryNet_FDR5e-5_0.8Expr_Ind0.10.mat'])
        atn{t} = binNet;
        atpartners(t, :) = (binNet(geneID, :)) + (binNet(:, geneID)');
        tsnet = finalTable(t).wholeNet;
        tspartners(t, :) = (tsnet(geneID, :)) + (tsnet(:, geneID)');
    end

    % first try: just by the count of the ts vs atn links for the
    % for a given gene, get its neighbours in the three ATNs 
    % 71. get the functional enrichment
    type = 'ATN';
    if strcmp(type, 'ATN')
        partners = atpartners; 
    else
        partners = tspartners;
    end

    enrichedF = zeros(5, length(inTerms));
    enrichedfp = zeros(5, length(inTerms));
    for t = 1 : 5
        tissue = tissues{t};
        geneSet = logical(partners(t,:));

        myMatP = GOdata.matP(geneSet, inFunctions);
        %        myMatC = GOdata.matC(geneSet, inFunctions);
        sum(geneSet);

        % sum of the functional mats (so all of them are together) 
        myTotalMat = myMatP; % + myMatC; %+myMatF;

        % count of genes for each of the functions
        mySumF = sum(myTotalMat);

        % getting p-value for each function
        geneSetCount = size(myTotalMat, 1);
        %            pVals = 1 - binocdf(mySumF, geneSetCount,
        %            totalSpec);p
        
        % mySumF(mySumF < 5) = 0;
        pVals = 1 - hygecdf(mySumF, 18494, geneCounts, sum(geneSet));

        % h = figure;
        % hist(pVals)
        % getting the FDR thr
        myTermInds = find(mySumF > 3); % functions with more than
                                       % one gene?
        myFCount = length(myTermInds);
        miniMat = myTotalMat(:, myTermInds);
        miniPVals = pVals(myTermInds);

        [a, b] = sort(miniPVals);
        sortedPVals = miniPVals(b);
        sortedTermIDs = myTermInds(b);
        counter = 0;

        % Count of observations
        for i = 1:length(miniPVals)
            if (sortedPVals(i) <= (alpha / fCount))                
                %                if (sortedPVals(i) <= i * alpha / fCount)
                %if (sortedPVals(i) <= i * alpha / length(myTermInds)) 
                counter = counter + 1; 
            end
        end
        length(miniPVals);
        counter;

        % my functions: 
        hitTermIDs = sortedTermIDs(1:counter);
        hitTerms = inTerms(hitTermIDs);
        enrichedF(t, hitTermIDs) = 1;
        enrichedfp(t, hitTermIDs) = pVals(hitTermIDs);
        tsfCount(g, t) = counter;
        tsfs(g, t, hitTermIDs) = 1;
    end
end

% report the extreme cases. :P 
%save(['~/resultsAndFigures/commonAndPure/' ...
%     'plusZeroTissueExpGenesFunctions.mat'], 'tsfs', '-v7.3')

save('~/resultsAndFigures/commonAndPure/multiTissueExpGenesFunctions_GO07_nonComp_expInOneT.mat', 'tsfs', ...
     '-v7.3');

load(['~/resultsAndFigures/commonAndPure/' ...
      'multiTissueExpGenesFunctions_GO07_nonComp_expInOneT.mat'])

% I should do the same thing for these genes: if they have big
% shift of function. 
shiftFunctionGenes = zeros(1, size(tsfs,1));
% for g = 1:size(tsfs, 1)
%     g
%     myMat = reshape(tsfs(g, :, :), [5, 8611]);
%     overlap = myMat * myMat'; % overlap of functions in the tissues
%     d = diag(overlap);
%     selectedT = find(d > 10); % finding the tissues which have more
%                               % functions for this gene
%     hitMat = zeros(5,5);
%     for i = 1:5
%         for j = i:5
%             hitMat(i,j) = overlap(i, j) /(overlap(i,i) + .05);
%             hitMat(j,i) = overlap(j, i)/(overlap(j,j) + .05);
%             if ((hitMat(i,j) < .05) + (hitMat(j, i) < .05) + (d(i) > ...
%                 10) + (d(j) > 10) == 4)
%                 shiftFunctionGenes(g) = 1;
%             end
%         end
%     end
% end
shiftFunctionGenes = zeros(1, size(tsfs,1));
sfInstances = zeros(2728, 6); % the count number here comes from
                              % the first time I run it (ran it
                              % once to get the count, ran it twice
                              % to fill this). Fields: gene
                              % secondary ID (maps to myGenes
                              % list), first tissue, second tissue,
                              % #FI in first, #FI in second, #overlap
count = 0;
for g = 1:size(tsfs, 1)
    g
    myMat = reshape(tsfs(g, :, :), [5, size(tsfs, 3)]);
    overlap = myMat * myMat'; % overlap of functions in the tissues
    d = diag(overlap);
    selectedT = find(d > 10); % finding the tissues which have more
                              % functions for this gene
    hitMat = zeros(5,5);
    for i = 1:5
        for j = i:5
            hitMat(i,j) = overlap(i, j) /(overlap(i,i) + .5);
            hitMat(j,i) = overlap(j, i)/(overlap(j,j) + .5);
            if ((hitMat(i,j) < .05) + (hitMat(j, i) < .05) + (d(i) > ...
                5) + (d(j) > 5) == 4)
                count = count + 1;
                shiftFunctionGenes(g) = shiftFunctionGenes(g) + 1;
                sfInstances(count, 1) =  g;
                sfInstances(count, 2) = i;
                sfInstances(count, 3) = j;
                sfInstances(count, 4) = d(i);
                sfInstances(count, 5) = d(j);
                sfInstances(count, 6) = overlap(i,j);
            end
        end
    end
end
sum(shiftFunctionGenes)
gs = find(shiftFunctionGenes);

% getting it with the JACCARD distance
% 3.2 using JCC distance: here is to get genes whith little SOF
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
count = 0;
jccSFInstances = zeros(12350, 7);
for g = 1:size(tsfs, 1)
    g
    myMat = reshape(tsfs(g, :, :), [5, size(tsfs, 3)]);
    overlap = myMat * myMat'; % overlap of functions in the tissues
    d = diag(overlap);
    hitMat = zeros(5,5);
    for i = 1:5
        for j = (i+1):5
            jcc = overlap(i,j)/ (d(i) + d(j) - overlap(i, j)+1);
            if(((d(j)>5) + (d(i) > 5)) ==2)
                count = count + 1;
                jccSFInstances(count, 1) =  g;
                jccSFInstances(count, 2) = i;
                jccSFInstances(count, 3) = j;
                jccSFInstances(count, 4) = d(i);
                jccSFInstances(count, 5) = d(j);
                jccSFInstances(count, 6) = overlap(i,j);
                jccSFInstances(count, 7) = jcc;
            end
        end
    end
end

bigJCC = find(jccSFInstances(:,7) < .1);
bigJCCInstances = jccSFInstances(bigJCC, :);

sib = unique((bigJCCInstances(:, 1)));
gpl570.uniqueSymbols(myGenes(sib))

% For each gene: get the cases where they could shift, get the
% cases where they don't. 
% their potential count of SHIFTS: 
book = sum(expMat');
realIDs = myGenes(sib);
theirExp = book(realIDs);

bars = hist(jccSFInstances(:,7), [0:.05:1]);
bar((bars))
set(gca, 'XTick', [1:21])
set(gca, 'XTickLabel',[0:.05:1])
title('Distribution of the jaccard distance for SOFI')
figFolder = '~/resultsAndFigures/functionalIdentityForGenes/figures/'
fileName = ['SOFI_JCC_hist'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);

% addpath('~/codes/MATLAB/myCodes/general/')
% heatmap(bigJCCInstances(1:100, 4:6), [],  [],  [], 'Colorbar', true, 'Colormap', ...
%          'summer')


% counting the potential cases of SOF
book = sum(expMat');
sum(book > 1) % count of genes expressed in more than one tissue
sum(book == 2) + sum(book == 3) * 3 + sum(book == 4) * 6 + sum(book ...
                                                  == 5) * 10 

thirdList = myGenes(logical(shiftFunctionGenes));

tlid = zeros(18494, 1);
tlid(thirdList) = 1;
sum(tlid)

slid = zeros(18494, 1);
slid(secondList) = 1;
sum(slid)

flid = zeros(18494, 1);
flid(firstList) = 1;
sum(flid)

sib = slid + tlid;
sum(sib == 2)

book = [flid slid tlid];
total = sum(book');
hits = find(total > 0);

% 3.5 Looking at the genes with high JCC in their SOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sutyding the SOF-JCC cases with similar Identities ... First look
% at the > .9, then > .8, then greater than >.7 and > .5 

% > .9
syms = GOdata.geneSymbols;
whos jccSFInstances
largerJCC = jccSFInstances(:, 7) > .4;
sum(largerJCC)
largerJCCIn = jccSFInstances(largerJCC, :);

largerJCCIn(1,:)
length(unique(syms(myGenes(largerJCCIn(:, 1)))))
length(unique(jccSFInstances(:,1)))

% >>> get the list of genes in each JCC bin, mark the ATP and RIB
% and other genes ... 
sofBinGenes = zeros(14301, 10);
s = 0 
for i = 1:10
    bt = (i - 1) / 10
    st = (i) / 10
    bigger = jccSFInstances(:, 7) >= bt;
    smaller = jccSFInstances(:, 7) < st;

    selected = (bigger + smaller)  == 2;
    % book = unique(jccSFInstances(selected, 1));
    % sofBinGenes(book, i) = 1; 
    s = s + sum(selected);
    selind = find(selected);
    for j = 1:length(selind)
        g = jccSFInstances(selind(j), 1);
        sofBinGenes(g, i) = sofBinGenes(g, i)  + 1;
    end
end

% get the genes with some SOF 
kado = sum(sofBinGenes') > 0;
sofInSBG = sofBinGenes(kado, :);

[a, b, c] = find(sofInSBG);
rs = (-1 + (2)* rand(length(b), 1)) / 3;

h = figure();
plot(b, c, '.')
plot((b + rs), c, '.') % there are more repeats at zero or smaller
                       % it is MORE common for genes to have big
                       % SOF between multiple tissues. 
                       
%getting the cumdist
cumDist = zeros(1, 10);
cumDist(1) = sum(sofBinGenes(:, 1) > 0);
for i = 2:10
    lowSOFSum = sum(sofBinGenes(:, [1:i])');
    loglss = lowSOFSum > 0;
    cumDist(i) = sum(loglss);
end

h = figure;
bar(cumDist)
xlim([0, 11])
set(gca, 'XTick', [1:10])
set(gca, 'XTickLabel',[.1:.1:1])
title('cum-dist of the count of genes with instances of minimum JCC')
figFolder = '~/resultsAndFigures/functionalIdentityForGenes/figures/'
fileName = ['cum_gene_SOFJCC_bar'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);
 
lowSOFSum = sum(sofBinGenes(:, [1])');
loglss = lowSOFSum > 0;
medSOFSum = sum(sofBinGenes(:, [3:7])');
logmss = medSOFSum > 0;
highSOFSum = sum(sofBinGenes(:, [6:10])');
loghss = highSOFSum > 0;


plotData.jccSFInstances = jccSFInstances;
plotData.sofBinGenes = sofBinGenes;
plotData.myGenes = myGenes
plotData.syms = syms;

save('~/resultsAndFigures/functionalIdentityForGenes/plotData_SOF.mat', ...
     'plotData')
load('~/resultsAndFigures/functionalIdentityForGenes/plotData_SOF.mat')


sib = loglss + loghss; % 512: very few genes change their SOF
                       % behavior : they either stay the same, or
                       % switch, it is not like the same gene would
                       % stay the same bewteen two tissues and then
                       % switch in the other one ...  This shows,
                       % multi tissue shift of NBFI is either rare,
                       % or zero.

kado = (sib == 2) + logmss;% 419, tis leaves 93 genes with wider
                           % shifts. I was right, genes do not diverge very
                       % much from their SOF behavior > they either
                       % switch much, or stay the same between the
                       % tissues. 

% >> now the group of genes which stay on the low
% exclusively. loglss - logmss - loghss
sib = loglss - logmss - loghss; 
sum(sib == 1) % 2592: 2/3 of the genes switch, Indeed, BUT, it is
              % shown that the switches are mostly between zero to
              % something, rather than between big and big. 

syms(myGenes(logical(sib)))

% >> group of genes which stay on the high exclusively. loghss -
% logmss
sib = loghss;
sum(loghss) % 632
sib = loghss; - loglss;
sum(sib == 1) % 48. 

book = syms(myGenes(sib == 1));
sofGenes = unique(jccSFInstances(:, 1));
sofGenesTID = myGenes(sofGenes);
a = sofGenesTID <= 13771;
b = sofGenesTID >= 13695;
sum((a + b) == 2) 

myList = syms(myGenes(jccSFInstances(selected, 1)));
unique(myList)
length(myList)

% 4. get the histograms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load count of enriched functions
load('~/resultsAndFigures/geneEnrichments/functionCount.mat')

load('~/resultsAndFigures/cGenes.mat')
 
cGenesInd = find(cGenes);

% confidently commonly expressed genes with shift of funciton in at
% least two tissues
cGenesSOF = firstList; % 324

% genes expressed in more than one tissue, with shift of function
% (loose selection of the other one) I have 12211 genes like this. 
oneTPlusSOF = thirdList; % 1599

shiftFlists.first = firstList';
shiftFlists.second = secondList; 
shiftFlists.third = thirdList;

save('~/resultsAndFigures/geneEnrichments/shiftFunctionGenes.mat', ...
     'shiftFlists')

load('~/resultsAndFigures/geneEnrichments/shiftFunctionGenes.mat')
% NOTE: those genes which are part of second but not third, are those
% which diverged from their pure partner, but did not diverge from
% their pure partner, but did not diverge from themselves much.

myGenes = shiftFlists.third;
saveFolder = ['~/resultsAndFigures/geneEnrichments/' ...
              'shiftOfFunctionEnrichments/morethanone/']
type = 'ATN'

for g = 1:length(myGenes)
    g
    ID = myGenes(g);
    functionalEnrichment_function(ID, type, saveFolder)
end

load('~/resultsAndFigures/commonAndPure/SOFPureNet.mat')

load('~/data/general/GOdata_GPL570_06.mat')

pureGenes = unique([SOFPureNet(:, 1), SOFPureNet(:, 2)]);

% >>>>>>> I realized I have 19 ribosome genes in the gene list, I
% wanted to look at the links for them. 
GOdata.geneSymbols{pureGenes(1:100)}
[a, b] = ismember('RPL13', GOdata.geneSymbols)

GOdata.geneSymbols{shiftFlists.second(1260:1278)}

% getting the rib genes which are part of pure!!! 
ribGenes = shiftFlists.second(1260:1278);

sGenes = ribGenes;
sGenes = 319;
[a1,b1] = ismember(SOFPureNet(:, 1), sGenes);
[a2,b2] = ismember(SOFPureNet(:, 2), sGenes);

links = SOFPureNet(((a1+a2)>0), :)

% so ADAR is differentiating from 40 of its brain links in
% other tissues. What are those genes doing in other tissues, wtf. 

% some ribosome links are marked as tissue specific. weird. 
l = 1;
[ GOdata.geneSymbols(SOFPureNet(l, 1))
 GOdata.geneSymbols(SOFPureNet(l, 2))]'
[(SOFPureNet(l, 1))
 (SOFPureNet(l, 2))]'

SOFPureNet(l, 3)
l = l + 1

% 5. get the null: links which are presenet in all the tissues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}
sumNets = zeros(18494, 18494);
nets = zeros(5, 18494, 18494);
for t = 1 : 5
    tissue = tissues{t};
    load( ['~/networks/tissues/' tissue '/' ...
           'binaryNet_FDR5e-5_0.8Expr_Ind0.10.mat'])
    nets(t, :, :) = binNet;
    sumNets = sumNets + binNet;
end

[a, b, c] = find(sumNets == 5);

sparseNets = [[a, b], zeros(length(a), 5)];
for i = 1:length(a)
    sparseNets(i, 3:7) = nets(:, a(i), b(i));
end

save('~/resultsAndFigures/tempResults/sparseNet_null.mat' , ...
     'sparseNets')
load('~/resultsAndFigures/tempResults/sparseNet_null.mat')

clear nets

% load(['~/resultsAndFigures/commonAndPure/' ...
%       'multiTissueExpGenesFunctions.mat'])

load(['~/resultsAndFigures/commonAndPure/' ...
      'multiTissueExpGenesFunctions_GO07_nonComp_expInOneT.mat'])

expMat = zeros(18494, 5);
for t = 1:5
    tissue = tissues{t};
    load(['~/data/general/tissueExpGenes/' tissue ...
          'ExpGenes0.8.mat'])
    expMat(:, t) = expGenesInd;
end

teCount = sum(expMat');
% for those genes, get this whole exp profile. 
% get the function of these genes in each of the tissues. 
tempMyGenes = find(teCount > 1);

% if the links are preseent in two tissues, how much do the genes
% overlap
%myGenes = unique([sparseNets(:,1), sparseNets(:,2)]);
%[aa, bb] = ismember(myGenes, tempMyGenes);
%mytsfs = tsfs(bb, :, :);
% clear tsfs
% tsfs = mytsfs;
% clear mytsfs; 

sz = length(a);
allol = zeros(sz, 5);
allmin = zeros(sz, 5);
shareFunction = zeros(sz, 5);
shareFunctionJK = zeros(sz, 5);
for i = 1:sz
    i
    g1 = a(i);
    g2 = b(i);
    tissueIDs = find(sparseNets(i, 3:7));
    g1Ind = find(myGenes == g1);
    g2Ind = find(myGenes == g2);
    g1fs = reshape(tsfs(g1Ind, :, :), [5, 3507]);
    g2fs = reshape(tsfs(g2Ind, :, :), [5, 3507]);
    
    % Identification of the change of function for true pure link
    % partners: pure link partners which share a lot of functions:
    % > 20
    
    diag1 = sum(g1fs');
    diag2 = sum(g2fs');
    
    overlap1 = g1fs * g1fs';
    overlap2 = g2fs * g2fs';
    overlaptwo = g1fs * g2fs';
    diag3 = diag(overlaptwo);
    shareFunction(i, tissueIDs) = diag3(tissueIDs);
    
    shareFunctionJK(i, tissueIDs) = shareFunction(i, tissueIDs) ./ (diag1(tissueIDs) + ...
        diag2(tissueIDs) - shareFunction(i, tissueIDs) + 1);

    % diag3 = diag(overlaptwo);
    
    % nd1 = atnds(:, g1);
    % nd2 = atnds(:, g2);
    
    % g1net = zeros(5, 18494);
    % g2net = zeros(5, 18494);
    % for t =1 :5
    %     fullNet = atn{t} + atn{t}';
    %     g1net(t, :) = fullNet(g1, :);
    %     g2net(t, :) = fullNet(g2, :);
    % end
    
    % netOverlap = g1net * g2net';
    % brainF = (g1fs(2, :) + g2fs(2,:)) ==2;    
    % liverF = (g1fs(3, :) + g2fs(3,:)) ==2;
    
    % similarity = liverF + brainF;
    % sum(similarity == 2);
    
    
    % 1. these links are pure, so the genes should not have similar
    % functions in other tissues. 
    
    % 2. The links, potentially, have other functions in other
    % tissues. 
    
    % 3. do they? Do they have any measureable difference in
    % functionality in other tissues?
    
    % 4. Both of them in different tissues?
    
    % If the two two genes have a lot in a tissue, but nothing in
    % common: any of them have more than 10 but none sharing.
    
    % let's see if they share < 1% 
    
    % compare the minimum with the shared. 
    % get the diag - 
    
    % checking what is the minimum fcount for each of the genes in
    % the other tissues.
    diags = [diag1(tissueIDs); diag2(tissueIDs)];
    minDs = min(diags);
    allmin(i, tissueIDs) = minDs;
    
    % checking what is the overlap of the genes in the other tissues
    allol(i, :) = diag(overlaptwo);
end

save('~/resultsAndFigures/tempResults/mainOverlapJK_null_nonComp_GO07.mat', ...
     'shareFunctionJK')

save('~/resultsAndFigures/tempResults/mainOverlap_null_nonComp_GO07.mat', ...
     'shareFunction')
save('~/resultsAndFigures/tempResults/allOverlap_null_nonComp_GO07.mat', ...
     'allol')
save('~/resultsAndFigures/tempResults/allmin_null_nonComp_GO07.mat', ...
     'allmin')

save('~/resultsAndFigures/tempResults/mainOverlap_null.mat', ...
     'shareFunction')
save('~/resultsAndFigures/tempResults/allOverlap_null.mat', ...
     'allol')
save('~/resultsAndFigures/tempResults/allmin_null.mat', ...
     'allmin')

% how many functions they share, in the tissue they have link in.
load('~/resultsAndFigures/tempResults/mainOverlap_null.mat')

% how many function they share (regardless if they share a tissue)
load('~/resultsAndFigures/tempResults/allOverlap_null.mat')

% what is their min function count in each tissue 

%load('~/resultsAndFigures/tempResults/allmin_null.mat')
clear
load('~/resultsAndFigures/tempResults/mainOverlapJK_null_nonComp_GO07.mat')

load('~/resultsAndFigures/tempResults/mainOverlap_null_nonComp_GO07.mat')
    
load('~/resultsAndFigures/tempResults/allOverlap_null_nonComp_GO07.mat')

load(['~/resultsAndFigures/tempResults/' ...
      'allmin_null_nonComp_GO07.mat'])

% getting the final plots for the null
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% >>>>> NUMBERS:

% >> How many links have > 5 FI in two tissues
passMin = allmin > 5;
sib = sum(passMin');
sum(sib > 1)
finalLinks = sib > 1;

% >>> GENERAL: get the base tissue and the other tissue. 
copyShareFunctionJCC = shareFunctionJK;
copyShareFunctionJCC(~passMin) = -1;
[maxJCC, maxInd] = max(copyShareFunctionJCC');

copyShareFunctionJCC = shareFunctionJK;
copyShareFunctionJCC(~passMin) = 1;
[minJCC, minInd] = min(copyShareFunctionJCC');

% >>>>>> PLOT3: Dist of minF for the base tissue for finalLinks: the tissue with
% max JCC 

minFs = zeros(size(maxJCC));
for i = 1:length(minFs)
    minFs(i) = allmin(i, maxInd(i));
end

h = figure
hist(minFs(finalLinks), 20)
title('Distribution of null base minFI for the finalLinks')
figFolder = '~/resultsAndFigures/functionalIdentityPureLinks/figures/'
fileName = ['nullBaseMinFICount_finalLinks_hist'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);

% >>>>>> PLOT1: The dist of JCC in the pure tissue. 
% we are looking at the final links. 

h = figure
hist(maxJCC(finalLinks), 20)
title('Null_ base _ Distribution of JCC for the finalHit')
figFolder = '~/resultsAndFigures/functionalIdentityPureLinks/figures/'
fileName = ['nullLinks_baseTissueJCC_finalHit__hist'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);
 
% >>>>>> PLOT2: Distribution of JCC in other links, for each
% linkMin
h = figure
hist(minJCC(finalLinks), 20)
title('Distribution of JCC in the other links, for eachlinkMin')
figFolder = '~/resultsAndFigures/functionalIdentityPureLinks/figures/'
fileName = ['nullLinks_otherTissueMinJCC_finalLinks_hist'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);
 
% >>>>>> PLOT4: Distribution of the FIcount for links from plot2
minFcounts = zeros(size(minInd));
for i = 1:length(minInd)
    minFcounts(i) = allmin(i, minInd(i));
end 
h = figure
hist(minFcounts(finalLinks), 20)
title('Distribution of FIcount in the other links, for eachlinkMin')
figFolder = '~/resultsAndFigures/functionalIdentityPureLinks/figures/'
fileName = ['nullLinks_otherTissueMinJCC_FIcount_finalLinks_hist'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf',[figFolder fileName '.pdf']);

% >>>>>> PLOT5: scatter of plots 3 and 4
h = figure
y = minFcounts(finalLinks);
x = minFs(finalLinks);
plot(x, y, '.')
corr(x', y') % they are correlated. The higher the first one, the
             % higher other one. How is it for the pure and semi
             % pure. This means that these are probably the group
             % of genes which do not switch function very much. 

% >>>>>> PLOT6: scatter of plots 1 and 2
w = maxJCC(finalLinks);
z = (minJCC(finalLinks));
plot(w, z, '.')
corr(w', z') % they are correlated, highly highly. The more similar
             % they are in base tissue, the more similar they are
             % in the other tissue. This shows that they do not
             % tend to change function, this is beautiful compared
             % to the pure and semi pure, since in pure and semi
             % pure they are around zero : random. the similarity
             % in the other tissue is RANDOM. 

% what it means: that pure links in general, do in fact say something about
% the functional similarity of the two genes involved. They are not
% just random.  
corr(x', z') % BUT in the null, they are correlated to the roof
             % (.73): the fact that the JCC is small in the other
             % tissue (4, z), is highly correlated with the small
             % FI in the alt tissue (2, x)

myMat = [x' y' w' z'];
nullXYZWmat = myMat;
save('~/resultsAndFigures/commonAndPure/nullXYZWmat.mat', 'nullXYZWmat')
sib = corr(myMat, 'type', 'Spearman')
addpath('~/codes/MATLAB/myCodes/general/')
h = figure
heatmap(sib, [],  [],  [], 'Colorbar', true, 'Colormap', ...
         'summer')
title('null')

h = figure
subplot(2, 2, 1)
hist(x, 20)
title('X')

subplot(2, 2, 2)
hist(y, 20)
title('Y')


subplot(2, 2, 3)
hist(w, 20)
title('W')

subplot(2, 2, 4)
hist(z, 20)
title('Z')

% end of plotting for the null
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% which of the links have the chance to have overlap in any tissue: 
passMins = allmin > 5;
sib = sum(passMins');

% which of them have the chance to diverge in another tissue (I can
% also get the overlap) >>>> cases of similarity vs cases of
% diverge 
sum(sib > 1)

% out of this, which of them show high similarity in one tissue: 
[maxJK, maxInds]= max(shareFunctionJK'); 
hist(maxJK(sib>1)) % genes with links presenet in more than one
                   % tissue, are much more likely to have bigger JC
                   % in one tissue ... But, how do they look in the
                   % other tissues?

% for those links, which have minimum of F for two tissues, I want
% to get the minimum overlap for the tissues they have function.
copyShareFunctionJK = shareFunctionJK;
copyShareFunctionJK(~passMins)= 1; % those without min F are now
                                    % 1, so they are max similar and
                                   % would not be selected

% getting the min for the overlap
[minJK, minInds]= min(copyShareFunctionJK'); 
hist(minJK(sib > 1)) % genes with links presenet in more than one

% obsolete from here on 
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% to find those who have overlap of function in one
% of the tissues they share link in: 
overlapMat = shareFunction > 5; 
sum(sum(overlapMat') > 0) / length(overlapMat) % these are .8021, while
                                           % for pure links it is
                                           % .771
% find diverge of function 
divergeMat = (shareFunction ./ (allmin + .03)) < .05;
c = 0;
hitlinks = zeros(1, length(overlapMat));
for i = 1 : length(overlapMat)
    % finding the tissues they have link in and share function in 
    overlapTs = find(overlapMat(i, :));
    for j = 1:length(overlapTs)
        myoverlap = divergeMat(i, :);
        myoverlap(overlapTs(j)) = []; % if in the other tissues,
                                      % they share a function ...
        mymin = passMins(i, :);
        mymin(overlapTs(j)) = [];
        if (max(mymin + myoverlap) == 2)
            c = c + 1;
            hitlinks(i) = hitlinks(i) + 1;
        end
    end
end

% hit links and c are different, as some of the links diverge more
% than once 
sum(hitlinks > 0) / length(overlapMat)
c / length(overlapMat) % only in 0.0226 of the cases where there is
                       % a link, we observe a shift of
                       % function. While for the pure links, we
                       % observed a shift of function for 0.25 of
                       % the links. 

% c should be divided by chances of overlap: 
book = sum(sparseNets(:, 3:end)');
instChances = 0;
sum(book == 5) * 5 * 4 /2 + ...
sum(book == 4) * 4 * 3 / 2 + ...
sum(book == 3) * 3 * 2 /2 + ...
sum(book == 2) 

c / (sum(sum(sparseNets(:, 3:end))) - length(sparseNets)) % 0.0422

sum(hitlinks) / ((sum(sum(sib))) - length(sib))
% for the commonlinks
load('~/resultsAndFigures/tempResults/sparseNet_null.mat')

sib = sparseNets(:, 3:7);

clinks = (sib == 5);

sib = passMins + passol;
sib = sib == 2;
halva = sum(sib') > 0;
halim = shareFunction > 10;
Inds = find((halva + halim') == 2);

% 6. show the diverge and play with it with the different thresholds. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6.1 load the pairs, with pure link in one tissue and shift of
% functionality 
% go up to section 2. to get the pure links 
clear
FDR = '0012'
FC = '3'

load(sprintf(['~/resultsAndFigures/TSlinks/finalTables/' ...
              'finalTable_CG13_FC%slog_FDR%s.mat'], FC, FDR))

% getting the ATND and TSND
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

% get all the genes in the pure net.
sumNet = zeros(18494, 18494);
for t = 1: 5
    tempPure = ((finalTable(t).cgNet + finalTable(t).noTSNet) == 2) ...
        .* t;
    sumNet = sumNet + tempPure;
end

book = sum(sumNet) + sum(sumNet');

[a, b, c] = find(sumNet);

load('~/resultsAndFigures/tempResults/mainOverlap.mat')
load('~/resultsAndFigures/tempResults/allOverlap.mat')
load('~/resultsAndFigures/tempResults/allmin.mat')

load('~/resultsAndFigures/tempResults/mainOverlapJK_GO07_nonComp.mat')
load('~/resultsAndFigures/tempResults/mainOverlap_GO07_nonComp.mat')
load('~/resultsAndFigures/tempResults/allOverlap_GO07_nonComp.mat')
load('~/resultsAndFigures/tempResults/allmin_GO07_nonComp.mat')

passMins = allmin > 5;
passol = (allol ./ (allmin + .05)) < .05;

sib = passMins + passol;
sib = sib == 2;
halva = sum(sib') > 0;
halim = shareFunction > 5; % <<<<<<===== THRESHOLD
% >>> ignore this part later
kad = sib(halim, :);
sum(sum(kad))  % we have 32128 incidents from the total of 80637*4
               % incidents : .0996
% <<< ignore this part later

% getting the links with two incidents: 
holu = sum(sib');
moreThanOne = holu > 1;
hit = find((moreThanOne + halim') == 2); 
hit = find((halva + halim') == 2); 

mynet = [a(hit), b(hit), c(hit)];

% selecting any tissue?
% mynetF1 = mynet(mynet(:, 3) == 2, :); 
% pureLinkMat = mynetF1;

% total count of pure links
pureLinkMat = mynet;
save('~/resultsAndFigures/TSlinks/pureLinkMat_V0.mat', ...
     'pureLinkMat')
load('~/resultsAndFigures/TSlinks/pureLinkMat_V0.mat')

save('~/resultsAndFigures/TSlinks/pureLinkMat_V1_GO07_noComp.mat', ...
     'pureLinkMat')

% for given two genes: 
% we want to have it's coexpression and it's
% functional profile in the networks: 

% >>>>>> load the gene function matrix ... 18k 8k
% in the network, in all the tissue networks. 
load('~/data/general/GOdata_GPL570_06.mat')
sumF = sum(GOdata.matF);
sumP = sum(GOdata.matP);
sumC = sum(GOdata.matC);
geneCounts = sumP + sumC;% + sumF;

% getting the 8611 functions in tsfs
inFunctions1 = geneCounts >= 5;
inFunctions2 = geneCounts <= 500;
inFunctions = (inFunctions1 + inFunctions2) == 2;
inFInds = find((inFunctions1 + inFunctions2) == 2);
fCount = sum(inFunctions)
inTerms = GOdata.GOTerms(inFInds);
inGOIDs = GOdata.GOID(inFInds);
geneCounts = geneCounts(inFInds);

fmat = GOdata.matP(:, inFInds)+ GOdata.matC(:, inFInds);
% >>>>>>> load the gene enriched function matrix ... 18k, 8 k
% get the genes expressed in more than one tissue (cause I have the
% enriched matrix for them only)

expMat = zeros(18494, 5);
for t = 1:5
    tissue = tissues{t};
    load(['~/data/general/tissueExpGenes/' tissue ...
          'ExpGenes0.8.mat'])
    expMat(:, t) = expGenesInd;
end
teCount = sum(expMat');
myGenes = find(teCount > 1);

% the enriched function in each tissue network
load(['~/resultsAndFigures/commonAndPure/' ...
      'multiTissueExpGenesFunctions.mat'])
fCount = size(tsfs, 3);

% %%% WITH THE NEW FILTER ON FUNCTIONS
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% I decided to add another filter: remove put filter
% for 10 and 200, finding the index of the new filter in the old
% one. 
inf1 = sum(fmat) >=10;
inf2 = sum(fmat) <=200;
inFnew = ((inf1 + inf2) == 2);
inFnewInds = find(inFnew); %***** THIS is important, pick the
                           %functions from this one 
inTermsNew = inTerms(inFnewInds);
inGOIDsNew = inGOIDs(inFnewInds);
geneCounts = geneCounts(inFnewInds);
newFmat = fmat(:, inFnewInds);

fmat = newFmat; % ***** fmat is fixed. 

% >>>>>>> load the gene enriched function matrix ... 18k, 8 k
% get the genes expressed in more than one tissue (cause I have the
% enriched matrix for them only)
expMat = zeros(18494, 5);
for t = 1:5
    tissue = tissues{t};
    load(['~/data/general/tissueExpGenes/' tissue ...
          'ExpGenes0.8.mat'])
    expMat(:, t) = expGenesInd;
end
teCount = sum(expMat');
myGenes = find(teCount > 1);

% the enriched function in each tissue network
load(['~/resultsAndFigures/commonAndPure/' ...
      'multiTissueExpGenesFunctions.mat'])
whos tsfs
tsfs = tsfs(:, :, inFnewInds); % ********* TSFS fixed. 
fCount = size(tsfs, 3);

inTerms = inTermsNew;
inGOIDs = inGOIDsNew;
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

% >>>> load the five ATN and TSN networks
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

% 6.2 extracting all the data for one pair
i = 122
g1 = pureLinkMat(i, 1)
g2 = pureLinkMat(i, 2)

i = 3
g1 = pureLinkMat(selectedLinks(i), 1)
g2 = pureLinkMat(selectedLinks(i), 2)

% getting the network neighbours
g1Net = zeros(5, 18494);
g2Net = zeros(5, 18494);
for t = 1:5
    net = atn{t};
    g1Net(t, :) = net(g1, :) + net(:, g1)';
    g2Net(t, :) = net(g2, :) + net(:, g2)';
end

sum(g1Net')
sum(g2Net')

% getting the enriched functions 
[a, b] = ismember(g1, myGenes)
g1Fs = reshape(tsfs(b, :, :), [5, fCount]);

[a, b] = ismember(g2, myGenes)
g2Fs = reshape(tsfs(b, :, :), [5, fCount]);

diag(g1Net * g2Net')'
diag(g1Fs * g2Fs')'

sum(g1Net')
sum(g2Net')

sum(g1Fs')
sum(g2Fs')
linkCode = i
i = i + 1
 
% 6.3 now, the network! 
% the link gens are connected to all their neighbours. But the
% neighbours connected if they share some enriched functions in
% either of the tissues: for each tissue we get a set of genes. 

% the links to be printed and count of enriched functions they
% share: all the neighbours: for each tissue a different network?
% ok. 
dataFolder = '~/data/general/'
load([dataFolder 'GPL570GemmaMapNEW.mat'])
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}

% building the sharedF matrix for the all the neighbor genes 
% for each pair of neighbors, see how many of thos Fs they share
% with each other
% writing the file 

clear tsfs
clear atn
clear tsnet
clear net 
clear passMins
clear pssol
clear shareFunction
clear tempPure
clear binNet
clear allmin
clear allol

for t = 1:5
    % GENE1 
    % for each pair of the neighbors, how many of their functions are
    % the enriched functions of the partner geens. 
    % get the enriched founctions of the partner genes:
    % get the functions
    templateFs = g1Fs(t, :);
    neighborFMatG1 = zeros(18494, 18494);
    neighborList = find(g1Net(t, :));
    len = length(neighborList)
    for i = 1:len
        i
        for j = (i+1):len 
            % getting the genes
            ni = neighborList(i);
            nj = neighborList(j);
            
            % getting their functions
            fi = fmat(ni, :);
            fj = fmat(nj, :);
            
            % getting their overlap of the shared for the three links
            % (partner with i, partner with j, i with j)
            neighborFMatG1(g1, ni) = max(fi * templateFs', .5);

            %            neighborFMatG1(ni, g1) = max(fi * templateFs', .5);

            neighborFMatG1(g1, nj) = max(fj * templateFs', .5);
            
            %            neighborFMatG1(nj, g1) = max(fj * templateFs', .5);

            neighborFMatG1(ni, nj) = max((((fj + fj) * templateFs') ...
                                          /2), .5);
            
            %            neighborFMatG1(nj, ni) = max(sum((fj + templateFs + fi) ...
            %                                == 3), .5);
        end
    end

    dg = diag(g1Fs * g2Fs');
    neighborFMatG1(g1, g2) = dg(t) ;
    
    % fix the matrix and merge the lower and upper indices
    tup = triu(neighborFMatG1, 1);
    tlo = tril(neighborFMatG1, -1);
    lowMat = tlo + tup';
    neighborFMatG1 = lowMat;
    clear lowMat tup tlo

    % GENE2
    % for each pair of the neighbors, how many of their functions are
    % the enriched functions of the partner geens. 
    % get the enriched founctions of the partner genes:
    templateFs = g2Fs(t, :);
    neighborFMatG2 = zeros(18494, 18494);
    neighborList = find(g2Net(t, :));
    len = length(neighborList)
    for i = 1:len
        i
        for j = (i+1):len
            % getting the genes
            ni = neighborList(i);
            nj = neighborList(j);
            
            % getting their functions
            fi = fmat(ni, :);
            fj = fmat(nj, :);
            
            % getting their overlap of the shared for the three links
            % (partner with i, partner with j, i with j)
                        
            % getting their overlap of the shared for the three links
            % (partner with i, partner with j, i with j)
            neighborFMatG2(g2, ni) = max(fi * templateFs', .5);

            %            neighborFMatG1(ni, g1) = max(fi * templateFs', .5);

            neighborFMatG2(g2, nj) = max(fj * templateFs', .5);
            
            %            neighborFMatG1(nj, g1) = max(fj * templateFs', .5);

            neighborFMatG2(ni, nj) = max((((fj + fj) * templateFs') ...
                                          /2), .5);
            
            %            neighborFMatG2(nj, ni) = max(sum((fj + templateFs + fi) ...
            %                                == 3), .5);

        end
    end
    
    % getting the g1 g2 value in G2
    neighborFMatG2(g1, g2) = dg(t) ;

    tup = triu(neighborFMatG2, 1);
    tlo = tril(neighborFMatG2, -1);
    lowMat = tlo + tup';
    neighborFMatG2 = lowMat;
    clear lowMat top tlo
    
    edgeThr = 3 
    book = (neighborFMatG2 > edgeThr) + (neighborFMatG1 ...
                                                      > edgeThr);
    sum(sum(book > 0))
    
    % making the link matrix    
    [a, b, c] = find(book);
    
    tissue = tissues{t};
    g1Sym = gpl570.uniqueSymbols{g1};
    g2Sym = gpl570.uniqueSymbols{g2};
    fileName = sprintf('%d_%s_%d_%d_%s_%s_edgeThr%.1f.txt', linkCode, tissue, g1, g2, g1Sym, ...
                       g2Sym, edgeThr)
    file = ['~/networks/cytoNetworks/functionalProfile/' fileName]
    fid = fopen(file, 'w')

    fprintf(fid, ['Gene1ID\t', 'Gene1Sym\t' 'Gene2ID\t', 'Gene2Sym\t', ...
                  'type\t', '#SharedEnrichedF1\t', '#SharedEnrichedF2\t', ...
                 '#sharedEnrichedMax\n'])
    for j = 1:length(a)
        sym1 = gpl570.uniqueSymbols{a(j)};
        sym2 = gpl570.uniqueSymbols{b(j)};
        type = 'f';
        if(((a(j) == g1) + (a(j) == g2) + (b(j) == g1) + (b(j) == g2)) > 0)
            type = 'c'
        end
        fprintf(fid, ['%d\t%s\t%d\t%s\t%s\t%d\t%d\t%d\n'], ...
                a(j), sym1, b(j), sym2, type,neighborFMatG1(a(j), ...
                                                          b(j)), ...
                neighborFMatG2(a(j), b(j)), max(neighborFMatG2(a(j), b(j)), neighborFMatG1(a(j), b(j))));
    end
    fclose(fid)
end

% 7. The examine the clusters from the clusterOne cytoscape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the problem with link 14 is that, they don't directly share a
% functional cluster of the genes in their tissue of pure
% links. Basically, their pure link contribute to small part of
% their functional identity. This coudl be learned from the
% difference between the count of genes they share in the networks,
% and the count of enriched functions they have : the genes they
% share do not contribute to a specific functional profile 

dataFolder = '~/data/general/'
load([dataFolder 'GPL570GemmaMapNEW.mat'])
n
file = '~/networks/localCyto/cluster2.txt'
t = 1
file = '~/networks/localCyto/blood_215_862.txt'

t = 2
file = '~/networks/localCyto/brain_215_862.txt'

t = 2

files = dir(['~/networks/localCyto/clusterOneResults_*' ...
             '_71_319_ABCF2_ADAR.csv'])

f = 4
file = ['~/networks/localCyto/' files(f).name]
fid = fopen(file)

% read one line 
tline = fgets(fid)

% parse the gene names: 
book = strsplit(tline, {', ', ' ' ,'\t', '"'});
[a, b] = ismember(book, gpl570.uniqueSymbols);
IDs = b(a);

% got the IDs, which functions they have? 
templateFs % has the enriched functions of g1. 
fmat % has the functions of the genes. 
tsfs % has the enriched functions of the genes. 

t = 4
templateFs2 = g2Fs(t, :);
templateFs1 = g1Fs(t, :);
both = (templateFs2 + templateFs1) == 2;

templateFs = both;
templateFs = templateFs2;
templateFs = templateFs1;
% now for my genes (IDs): 
sum(templateFs)
myfMat = fmat(IDs, logical(templateFs));
sum(myfMat)
interestingTerms = find(sum(myfMat) > (5))
myTerms = inTerms(logical(templateFs));
myTermIDs = find(templateFs);

finalt = myTerms(interestingTerms)
finalID = myTermIDs(interestingTerms)
finalt(:)

addpath('~/codes/MATLAB/myCodes/general')
h = figure; 
heatmap(myfMat, [],  [],  [], 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'summer')

% and we print the gene IDs and the gene Names 

% and for the links, we have their count: but also, which functions
% for each link? not now. I want them printable right away... 


% 7.1 For a given pair of genes, get the clusters for each
% tissue. The tissue orders are blood brain liver lung sm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
files = dir(['~/networks/localCyto/clusterOneResults_*' ...
             '_71_319_ABCF2_ADAR.csv'])

% get the list of all genes: 
genesInClusters = zeros(18494, 5);
for f = 1:length(files)
    file = ['~/networks/localCyto/' files(f).name]
    fid = fopen(file)
    % read one line 
    tline = fgets(fid)
    tline = fgets(fid)
    c = 1
    key = true
    while ((sum(tline) ~= -1) && key)
        book = strsplit(tline, {', ', ' ' ,'\t', '"'});
        [a, b] = ismember(book, gpl570.uniqueSymbols);
        IDs = b(a);
        if (sum(a) <= 5)
            key = false
        end
        genesInClusters(IDs , f) = c;
        tline = fgets(fid)
        c = c + 1;
    end
end 

genesIn = sum(genesInClusters') > 0;
geneIDs = find(genesIn);
finalMat = [geneIDs' genesInClusters(geneIDs, :)];

fileName = 'ADAR_ABCF2_71_319.txt'
file = ['~/networks/cytoNetworks/nodeAttribute/' fileName]
fid = fopen(file, 'w')

% printing the node attributes: 
fprintf(fid, ['Gene1ID\t', 'Gene1Sym\t' 'blood\t', 'brain\t', ...
              'liver\t', 'lung\t', 'muscle\n'])
for i = 1:length(finalMat)
    sym = gpl570.uniqueSymbols{finalMat(i, 1)};
    b = finalMat(i, :);
    fprintf(fid, ['%d\t%s\t%d\t%d\t%d\t%d\t%d\n'], b(1), sym, b(2), ...
            b(3), b(4), b(5), b(6))
end
fclose(fid)

% 8. the plots for the FI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
% >>>> loading and preparing data
load(['~/resultsAndFigures/commonAndPure/' ...
      'multiTissueExpGenesFunctions_GO07_nonComp_expInOneT.mat'])

% >>>> loading and preparing data
load(['~/resultsAndFigures/commonAndPure/' ...
      'multiTissueExpGenesFunctions.mat'])

% Tissues And final table
FC = '3'
FDR = '0012'
load(sprintf(['~/resultsAndFigures/TSlinks/finalTables/' ...
              'finalTable_CG13_FC%slog_FDR%s.mat'], FC, FDR))
tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'};

% Function matrix
gCount = 18494;
load('~/data/general/GOdata_GPL570_07_nonComp.mat')
% load('~/data/general/GOdata_GPL570_07.mat')
% load('~/data/general/GOdata_GPL570_06.mat')
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
inGOIDs = GOdata.GOID(inFInds);
geneCounts = geneCounts(inFInds);
%save('~/data/clusterings/inTerms8611.mat', 'inTerms')

% function specificity:
%sumF = sum(GOdata.matF(:, inFunctions));
sumP = sum(GOdata.matP(:, inFunctions));
%sumC = sum(GOdata.matC(:, inFunctions));

fMat = GOdata.matP(:, inFunctions);

%specF = sumF / gCount;
specP = sumP / gCount;
%specC = sumC / gCount;

totalSpec =  specP;% + specC;% + specF;

% get the genes expressed in more than one tissue. 
expMat = zeros(18494, 5);
for t = 1:5
    tissue = tissues{t};
    load(['~/data/general/tissueExpGenes/' tissue ...
          'ExpGenes0.8.mat'])
    expMat(:, t) = expGenesInd;
end

teCount = sum(expMat');
% for those genes, get this whole exp profile. 
% get the function of these genes in each of the tissues. 
myGenes = find(teCount > 0);

% >>>> load the five ATN and TSN networks
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

% 8.1 for the genes which are in the network for each tissue, get the
% dist of FI, 
% >>>>>>>>>>>>>>> TESTS for the dist of FI, fcoun and ND dists,
% their correlatin and stuff

for t = 1:5
    % get the genes that are xpressed in blood.
    tgenes = find(atnds(t, :) > 0);
    [a, b] = ismember(tgenes, myGenes);
    myMat = reshape((tsfs(b(a), 1, :)), [sum(a), size(tsfs, 3)]);
    FIcount = sum(myMat'); % count of FI for each gene in blood
    hist(log10(FIcount+1), 40)
    kado = fMat(tgenes, :);
    tfCounts = sum(kado');
    
    fiStr(t).tissue = tissues{t};
    fiStr(t).inNetGenes = tgenes;
    fiStr(t).FIcount = FIcount;
    fiStr(t).tfCounts = tfCounts;
    fiStr(t).nds = atnds(t, tgenes);
end

save(['~/resultsAndFigures/functionalIdentityForGenes/' ...
      'tissuesFIgeneralInf.mat'], 'fiStr')

load(['~/resultsAndFigures/functionalIdentityForGenes/' ...
      'tissuesFIgeneralInf.mat'])

h = figure
hist(log10(tfCounts+1), 40)

tnds = atnds(t, tgenes);
h = figure
log10(hist(tnds+1), 40)

mytsnds = tsnds(t, tgenes);
h = figure
hist(log10(mytsnds+1), 40)

plot(tnds(a), FIcount, '.', 'color',c2)
corr(tnds', FIcount')
h = figure
plot(tnds, tfCounts, '.')
corr(tnds', tfCounts')
plot(tfCounts, FIcount, '.')
% get the corr of FI and ND

i = 30

myGenes(hitGenes(i))

hitGenes = find(shiftFunctionGenes);
book = reshape(tsfs(hitGenes(i), :, :), [5, 3507]);

book * book'
i = i + 1

acadsb 17
acat1 18

% 8.2 PLOT
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% for each tissue, get the FIdist, NDdist, fcountDist. 

% How much do the TS links contribute to the FI.

% How many genes show shift of functionality. 

% Get some role of the TS links. 

% Get the shift of functionality for instances: the overlap, the
% switch: ok. 
shitfFunctionGenes;
sfInstances; % this one has the valued shift. 
myGenes; % this one has the real gene IDs for the secondary IDs in
         % the sfInstances
tsfs; % this one has the FI
fMat; % this one has the GOL: GO label

n = 100
h = figure
plot(sfInstances(1:n, 4), '*')
xlim([0, n + 1])
ylim([0, 100])
hold on
plot(sfInstances(1:n, 5), '*')
hold on
plot(sfInstances(1:n, 6), '*', 'color','r')

% now sort them, put them to lines

% get their max 
% sort them by their max

sib = max(sfInstances(:, 4:5)');
[a, b] = sort(sib, 'descend');

sortedSFInst = sfInstances(b, :);
% draw a line from the overlap to min : black 
% draw a line from the overlap to max : blue

n = length(sortedSFInst)

s = 1
e = n
h = figure
hold all
for i = s:e
    y1 = sortedSFInst(i, 6);
    y2 = min(sortedSFInst(i, 4:5));
    y3 = max(sortedSFInst(i, 4:5));
    line([i, i], [y1, y2], 'LineWidth', 3,'Color', 'r')
    line([i, i], [y2, y3], 'LineWidth', 3 ,'Color', 'k')
    line([i, i], [0, y1], 'LineWidth', 3 ,'Color', 'b')
end
xlim([0, n+1])
s = s + 100
xlim([0, e+1])

sortedSFInst([8, 10, 15], :)
myGenes(sortedSFInst([8 , 10, 15], 1))

% which of them have FI assigned to them, in any of the tissues? 

% for any of the instances, check to see that. 
hit = zeros(length(sortedSFInst), 1);
for i = 1: length(sortedSFInst)
    g = sortedSFInst(i, 1);
    t1 = sortedSFInst(i, 2);
    t2 = sortedSFInst(i, 3);
    FI1 = reshape(tsfs(g, t1, :), [1, size(tsfs, 3)]);
    FI2 = reshape(tsfs(g, t2, :), [1, size(tsfs, 3)]);
    fs = fMat(myGenes(g), :);
    sortedSFInst(i, :);
    sum(fs);
    if sum((fs + FI1) == 2) > 0 
        hit(i) = hit(i) + 1;
    end
    if  sum((fs + FI2) == 2) > 0
        hit(i) = hit(i) + 1;
    end
end

sortedSFInst(i, :)
GOdata.geneSymbols(myGenes(g))
i

myNet = zeros(5, 18494);
geneID = myGenes(g);
for t = 1:5
    net = atn{t};
    myNet(t, :) = net(geneID, :) + net(:, geneID)';
end

myNet * myNet'

FI1TermIDs = find((fs + FI1) == 2)
[inTerms(FI1TermIDs)]
inTerms(FI1TermIDs(1))

FI2TermIDs = find((fs + FI2) == 2)
inTerms(FI2TermIDs)
inTerms(FI2TermIDs(3))

hitInds = find(hit == 2);
%hind = 7
hitInds(hind)
i = hitInds(hind)
hind = hind + 1

%>>> For a given gene, give me the count of TSN and ATN links in
% any of the instances: 

theDs = zeros(length(sortedSFInst), 4);
for i = 1 : length(sortedSFInst)
    id = myGenes(sortedSFInst(i, 1))
    t1 = sortedSFInst(i, 2);
    t2 = sortedSFInst(i, 3); 
    atnd1 = atnds(t1, id);
    tsnd1 = tsnds(t1, id);
    atnd2 = atnds(t2, id);
    tsnd2 = tsnds(t2, id);
    theDs(i, 1) = atnd1; % - tsnd1;
    theDs(i, 2) = tsnd1;
    theDs(i, 3) = atnd2; %- tsnd2;
    theDs(i, 4) = tsnd2;
end

% modifications on the theDs for plotting:
logDs = log2(theDs + 1);

% 1.
[a, b] = sort([logDs(:, 1) + logDs(:, 3)]);
soloDs = logDs(b, :);

% 2. 
temp = max(logDs(:, [2, 4])');
[a, b] = sort(temp);
mmsoloDs = logDs(b, :);

% selecting one (also change the plot loop)
Ds = mmsoloDs;
Ds = soloDs;
sp = .5;

c1 = [200, 200, 200]/255
c2 = [80, 80, 80]/255
h = figure
set(h, 'Position', [100, 100, 450, 750])
hold all
for i = 1:length(Ds)
    x01 = sp;
    x02 = -sp
    
    % >> way one
    % x1 = (Ds(i, 1));
    % x2 = (Ds(i, 2));
    % x3 = Ds(i, 3);
    % x4 = Ds(i, 4);
    
    % way two % works with mmsoloDs
    [a, b] = max([Ds(i, [2, 4])]);
    [c, d] = min([Ds(i, [2, 4])]);
    x1 = (Ds(i, b*2-1));
    x2 = a;
    x3 = Ds(i, d*2-1);
    x4 = c;
    
    line([x01, x1+sp], [i, i], 'LineWidth', 1,'Color', c1)
    line([x01, x2+sp], [i, i], 'LineWidth', 1 ,'Color', c2)
    line([x02, -x3-sp], [i, i], 'LineWidth', 1 ,'Color', c1)
    line([x02, -x4-sp], [i, i], 'LineWidth', 1 ,'Color', c2)
end 
ylim([-90, 2800])

set(gca, 'XTick', [-15.5, -10.5, -5.5, -.5, .5, 5.5, 10.5, 15.5 ...
                   ])
set(gca, 'XTickLabel', {'15', '10', '5', '0', '0', '5', '10', ...
                    '15'})
xlim([-12, 12])
set(h, 'PaperPositionMode', 'auto')

figFolder = '~/resultsAndFigures/functionalIdentityForGenes/figures/'
fileName = ['SOF_Incidents_GO07_nonComp_5functionThr_maxTsndSorted'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf', [figFolder fileName '.pdf']);


% 8.3 PLOT: for the extreme cases of SOF, 
% give me the expression level of the genes involved in each of the
% tissues. Give me the colors as well, you can do it in a heatmap
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% > preparation: expression variance for each tissue: 

load('~/data/general/tissueExpGenes/allDSExprInfo.mat')
myMat = zeros(18494, 53);
for i = 1: 53
    myMat(:, i) = allDSExprInfo(i).exprLevel;
end

myMat = quantilenorm(myMat);

avgExp = zeros(18494, 5);
varExp = zeros(18494, 5);

avgExp(:, 1) = mean(myMat(:, 13:21)'); % blood
avgExp(:, 2) = mean(myMat(:, 1:12)'); % brain
avgExp(:, 3) = mean(myMat(:, 29:38)'); % liver
avgExp(:, 4) = mean(myMat(:, 39:53)'); % lung
avgExp(:, 5) = mean(myMat(:, 22:28)'); % muscle

varExp(:, 1) = var(myMat(:, 13:21)'); % blood
varExp(:, 2) = var(myMat(:, 1:12)'); % brain
varExp(:, 3) = var(myMat(:, 29:38)'); % liver
varExp(:, 4) = var(myMat(:, 39:53)'); % lung
varExp(:, 5) = var(myMat(:, 22:28)'); % muscle

myGenes; % the list of 14301 genes expressed in more than one
         % tissue. 
sfInstances; % cases of filtered SFO
tsfs; % the enriched Fs for each of the genes in each of the tissues
fMat; % the GOFs for each of the genes
atnds; %
tsnds; %
atn; % the whole networks, to identify which genes are involved in
     % the identity 

% how would you sort the SHIFT in the networks 

SOFutil.myGenes = myGenes;
SOFutil.sfInstances = sfInstances;
SOFutil.tsfs = tsfs;
SOFutil.fMat = fMat;
SOFutil.atnds = atnds;
SOFutil.tsnds = tsnds;

save('~/resultsAndFigures/functionalIdentityForGenes/sofUtil.mat', ...
     'SOFutil')

load('~/resultsAndFigures/functionalIdentityForGenes/sofUtil.mat')

myGenes = SOFutil.myGenes;
sfInstances = SOFutil.sfInstances;
tsfs = SOFutil.tsfs;
fMat = SOFutil.fMat;
atnds = SOFutil.atnds;
tsnds = SOFutil.tsnds;

thr = 6
t1pass = sfInstances(:, 4) >= thr;
t2pass = sfInstances(:, 5) >= thr;

selected = ((t1pass + t2pass) == 2);
selectedInst = sfInstances(selected, :);

load(['~/resultsAndFigures/functionalIdentityForGenes/' ...
      'plotData_SOF.mat'])
selectedInst = plotData.jccSFInstances;

% for each instance, get the gene partners which have the enriched
% functions, in the given tissue.
selectedInst = bigJCCInstances;
selectedInst = sfInstances;
lineVs = zeros(length(selectedInst), 4);
for i = 1: length(selectedInst)
    mgid = selectedInst(i, 1);
    gid = myGenes(mgid);
    t1 = selectedInst(i, 2);
    t2 = selectedInst(i, 3);
    
    fi1 = reshape(tsfs(mgid, t1, :), [1, 3507]);
    fi2 = reshape(tsfs(mgid, t2, :), [1, 3507]);
    
    net1 = atn{t1};
    nghs1 = net1(gid, :) + net1(:, gid)';
    
    net2 = atn{t2};
    nghs2 = net2(gid, :) + net2(:, gid)';
    
    % getting the genes involved in FI1
    n1fMat = fMat(logical(nghs1), :);
    nghs1IDs = find(nghs1);
    temp = (n1fMat * fi1');
    funcNghs1 = nghs1IDs(logical(temp)); % list of gene IDs ,
                                         % involving in the FI1
    
    % getting the genes involved in FI2
    n2fMat = fMat(logical(nghs2), :);
    nghs2IDs = find(nghs2);
    temp = (n2fMat * fi2');
    funcNghs2 = nghs2IDs(logical(temp)); % list of gene IDs
                                         % involving in the FI2
    
    % now their average and their variance for tissue exp: load
    % that file ...
    gExpfi1t1 = avgExp(funcNghs1, t1);
    gExpfi1t2 = avgExp(funcNghs1, t2);
    % [a, b] = sort(gExpfi1t2);
    % sgef1t2 = a;
    % [a, b] = sort(gExpfi1t1);
    % sgef1t1 = a;
    lineVs(i, 1) = length(gExpfi1t1);
    lineVs(i, 2) = sum(gExpfi1t2 < 4.9393);
    
    gExpfi2t1 = avgExp(funcNghs2, t1);
    gExpfi2t2 = avgExp(funcNghs2, t2);   
    % [a, b] = sort(gExpfi2t1);
    % sgef2t1 = a;
    % [a, b] = sort(gExpfi2t2);
    % sgef2t2 = a;
    lineVs(i, 3) = length(gExpfi2t2);
    lineVs(i, 4) = sum(gExpfi2t1 < 4.9393);
    
    % h = figure; 
    % subplot(1, 2, 1)
    % sib = [sgef1t1, sgef1t2]
    % %sib = [gExpfi2t1 gExpfi2t2]
    % heatmap(sib, [],  [],  [], 'TickAngle', 45, 'ShowAllTicks', ...
    %         true, 'Colorbar', true, 'Colormap', 'summer')
    % subplot(1, 2, 2)
    % sib = [sgef2t2, sgef2t1]
    % heatmap(sib, [],  [],  [], 'TickAngle', 45, 'ShowAllTicks', ...
    %         true, 'Colorbar', true, 'Colormap', 'summer')
end

save('~/resultsAndFigures/functionalIdentityForGenes/allLineVs.mat', ...
     'lineVs')

addpath('~/codes/MATLAB/myCodes/general/')

sib = [sgef1t1, sgef1t2]
sib = [sgef2t2, sgef2t1]
%sib = [gExpfi2t1 gExpfi2t2]
h = figure; 
heatmap(sib, [],  [],  [], 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'summer')

% if you mark it as NOT EXPRESSED, then you can do the line thing
% like you did for the other thing, using two colors 

h = figure
hold all
sp = .5
set(h, 'Position', [100, 100, 450, 750])
c1 = [200, 200, 200]/255
c2 = [80, 80, 80]/255
for i = 1:length(lineVs)
    x01 = sp;
    x02 = -sp
    
    x1 = log2(lineVs(i, 1) + 1);
    x2 = log2(lineVs(i, 2) + 1);;
    x3 = log2(lineVs(i, 3) + 1);
    x4 =  log2(lineVs(i, 4) + 1);
    
    line([x01, x1+sp], [i, i], 'LineWidth', 2,'Color', c1)
    line([x01, x2+sp], [i, i], 'LineWidth', 2 ,'Color', c2)
    line([x02, -x3-sp], [i, i], 'LineWidth', 2 ,'Color', c1)
    line([x02, -x4-sp], [i, i], 'LineWidth', 2 ,'Color', c2)
end 
title(['the expressed vs non expressed genes in the opposite ' ...
       'tissues'])
ylim([0, 101])

% getting the NDs for the distances
sfInstances = bigJCCInstances;
theDs = zeros(length(sfInstances), 4);
for i = 1 : length(sfInstances)
    id = myGenes(sfInstances(i, 1))
    t1 = sfInstances(i, 2);
    t2 = sfInstances(i, 3); 
    atnd1 = atnds(t1, id);
    tsnd1 = tsnds(t1, id);
    atnd2 = atnds(t2, id);
    tsnd2 = tsnds(t2, id);
    theDs(i, 1) = atnd1; % - tsnd1;
    theDs(i, 2) = tsnd1;
    theDs(i, 3) = atnd2; %- tsnd2;
    theDs(i, 4) = tsnd2;
end

selectedInst;
selectedDs = theDs(selected, :);

selectedAll = [selectedInst selectedDs];

sp = .5
Ds = selectedDs;
h = figure
set(h, 'Position', [100, 100, 450, 750])
hold all
for i = 1:length(Ds)
    x01 = sp;
    x02 = -sp
    
    % >> way one
    x1 = log2(Ds(i, 1));
    x2 = log2(Ds(i, 2));
    x3 = log2(Ds(i, 3));
    x4 = log2(Ds(i, 4));
    
    % way two % works with mmsoloDs
    % [a, b] = max([selectedDs(i, [2, 4])]);
    % [c, d] = min([selectedDs(i, [2, 4])]);
    % x1 = (Ds(i, b*2-1));
    % x2 = a;
    % x3 = Ds(i, d*2-1);
    % x4 = c;
    
    line([x01, x1+sp], [i, i], 'LineWidth', 2,'Color', c1)
    line([x01, x2+sp], [i, i], 'LineWidth', 2 ,'Color', c2)
    line([x02, -x3-sp], [i, i], 'LineWidth', 2 ,'Color', c1)
    line([x02, -x4-sp], [i, i], 'LineWidth', 2 ,'Color', c2)
end 
title('tsn and tns of the top 100 SOF')
ylim([0, 101])

% 8.35 studying the selected instances: 
% >>>>>
whos myGenes % list of genes 
whos selectedDs % node degree af atn1, tsn1, at2, tsn2
whos lineVs % exp of t1, exp of g2 in t1, exp of t2, exp of g1 in t2
whos selectedInst % the instances
gpl570

G = 3
ids = myGenes(selectedInst(:,1));
syms = gpl570.uniqueSymbols(ids);
syms(G)
whos syms
selectedInst(G, :)

g1 = ids(G)
sum(fMat(g1, :))

disID1 = selectedInst(G,1)
t = 2
book = inTerms(logical(reshape(tsfs(disID1, t, :), [1, 3507])));
book{:}

g1 = ids(G)
% get the foverlap for fMat
fs = fMat(g1, :);
OS = zeros(5, 1);
overlaps = zeros(3507, 5);
for i = 1:5
    temp = reshape(tsfs(disID1, i, :), [1, 3507]);
    OS(i) = temp * fs';
    overlaps(:, i) = (temp + fs) == 2;
end

overlaps' * overlaps

% getting the network neighbours
g1Net = zeros(5, 18494);
for t = 1:5
    net = atn{t};
    g1Net(t, :) = net(g1, :) + net(:, g1)';
end

sum(g1Net')
g1Net * g1Net'
% getting the enriched functions 
[a, b] = ismember(g1, myGenes)
g1Fs = reshape(tsfs(b, :, :), [5, fCount]);

g1Fs * g1Fs'
sum(g1Net')

sum(g1Fs')

% together = [sfInstances theDs];
% selectedto = together(selected, :);

% 8.4 PLOT: get the JACCARD dist and the portion DIST for the genes expressed
% in more than one tissue. 
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% what we need: 
myGenes; % the 14301 genes expressed in all the tissues
fMat; % the GO mat 
tsfs; % the tsfs mat for each of the tissues. 

allFIMat = reshape(sum(tsfs, 2), [14301, 3507]);
allFIMat(allFIMat > 0) = 1;

myfMat = fMat(myGenes, :);

fifOverlapCount = diag(allFIMat * myfMat'); % overlap count of f
                                            % and FI

fCounts = sum(myfMat'); % count of f for the genes
sum(fCounts>0)
allFICounts = sum(allFIMat'); % count of FI for the genes
sum(allFICounts >0)

% get the count of genes with both f and FI
sum(((fCounts>0) + (allFICounts>0)) == 2) 
myGenes2 = find((((fCounts>0) + (allFICounts>0)) == 2)); % index of
                                                         % genes in 14301

% get the jaccard: 
myGenes2fsOverlap = fifOverlapCount(myGenes2);
myGenes2fsSum = fCounts(myGenes2) + allFICounts(myGenes2);
JCC = myGenes2fsOverlap' ./ (myGenes2fsSum - myGenes2fsOverlap');

[a, b] = ksdensity(JCC, 'support', [-.01, 1.01])
plot(b, a)

% plotting the Jaccard distance of fi and f
h = figure
mybars = hist(JCC,  [0:.05:1]);
bar(log2(mybars))
set(gca, 'XTick', [0:2:21], 'XTickLabel', [0:.1:1])
xlim([0, 21.1])
figFolder = '~/resultsAndFigures/functionalIdentityForGenes/figures/'
fileName = ['fif_JCC'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf', [figFolder fileName '.pdf']);


% the dist of all FI
h = figure
allFICounts;
mybars = hist((allFICounts), [0:10:200])
bar(log2(mybars))
set(gca, 'XTick', [1:2:21], 'XTickLabel', [0:20:200])
xlim([0 22])
figFolder = '~/resultsAndFigures/functionalIdentityForGenes/figures/'
fileName = ['allFICounts'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf', [figFolder fileName '.pdf']);

% the dist of all FC 
h = figure
fCounts;
mybars = hist((fCounts), [0:10:200])
bar(log2(mybars))
set(gca, 'XTick', [1:2:21], 'XTickLabel', [0:20:200])
xlim([0 22])
figFolder = '~/resultsAndFigures/functionalIdentityForGenes/figures/'
fileName = ['fCounts'];
print(h, '-depsc', [figFolder fileName '.eps']);
print(h, '-dpdf', [figFolder fileName '.pdf']);

% 8.5 just getting everything for local plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% based on sfInstances and myGenes (no sorting and no nothing)
sfInstances;
myGenes;

% get the Ds 
% getting the NDs for the distances
theDs = zeros(length(sfInstances), 4);
for i = 1 : length(sfInstances)
    id = myGenes(sfInstances(i, 1))
    t1 = sfInstances(i, 2);
    t2 = sfInstances(i, 3); 
    atnd1 = atnds(t1, id);
    tsnd1 = tsnds(t1, id);
    atnd2 = atnds(t2, id);
    tsnd2 = tsnds(t2, id);
    theDs(i, 1) = atnd1; % - tsnd1;
    theDs(i, 2) = tsnd1;
    theDs(i, 3) = atnd2; %- tsnd2;
    theDs(i, 4) = tsnd2;
end


% get lineVs (expression of the genes which are part of the
% functionality)
selectedInst = sfInstances;
lineVs = zeros(length(selectedInst), 4);
for i = 1: length(selectedInst)
    mgid = selectedInst(i, 1);
    gid = myGenes(mgid);
    t1 = selectedInst(i, 2);
    t2 = selectedInst(i, 3);
    
    fi1 = reshape(tsfs(mgid, t1, :), [1, 3507]);
    fi2 = reshape(tsfs(mgid, t2, :), [1, 3507]);
    
    net1 = atn{t1};
    nghs1 = net1(gid, :) + net1(:, gid)';
    
    net2 = atn{t2};
    nghs2 = net2(gid, :) + net2(:, gid)';
    
    % getting the genes involved in FI1
    n1fMat = fMat(logical(nghs1), :);
    nghs1IDs = find(nghs1);
    temp = (n1fMat * fi1');
    funcNghs1 = nghs1IDs(logical(temp)); % list of gene IDs ,
                                         % involving in the FI1
    
    % getting the genes involved in FI2
    n2fMat = fMat(logical(nghs2), :);
    nghs2IDs = find(nghs2);
    temp = (n2fMat * fi2');
    funcNghs2 = nghs2IDs(logical(temp)); % list of gene IDs
                                         % involving in the FI2
    
    % now their average and their variance for tissue exp: load
    % that file ...
    gExpfi1t1 = avgExp(funcNghs1, t1);
    gExpfi1t2 = avgExp(funcNghs1, t2);
    lineVs(i, 1) = length(gExpfi1t1);
    lineVs(i, 2) = sum(gExpfi1t2 < 4.9393);
    
    gExpfi2t1 = avgExp(funcNghs2, t1);
    gExpfi2t2 = avgExp(funcNghs2, t2);   
    lineVs(i, 3) = length(gExpfi2t2);
    lineVs(i, 4) = sum(gExpfi2t1 < 4.9393);
end

% We want to show what role TS links play in the SOF. We have the
% TSN ATN thing, sorted nicely: (got this part from above)
logDs = log2(theDs + 1);

% 1.
[a, b] = sort([logDs(:, 1) + logDs(:, 3)]);
soloDs = logDs(b, :);

% 2. 
temp = max(logDs(:, [2, 4])');
[a, b] = sort(temp);
mmsoloDs = logDs(b, :);

% sort the sfInstances and lineVs
mmsoloSFI = sfInstances(b, :);
maxf = zeros(length(sfInstances), 1);
minf = zeros(length(sfInstances), 1);

mmsoloLineVs = lineVs(b, :);
fixedLineVs = zeros(size(lineVs));

extSOF.myGenes = myGenes;
extSOF.mmsoloSFI = mmsoloSFI;
extSOF.mmsoloLineVs = mmsoloLineVs;
extSOF.mmsoloDs = mmsoloDs;
save('~/resultsAndFigures/functionalIdentityForGenes/lineVs(exp)_Ds_sorted.mat', ...
     'extSOF')

load(['~/resultsAndFigures/functionalIdentityForGenes/' ...
      'lineVs(exp)_Ds_sorted.mat'])
myGenes = extSOF.myGenes;
mmsoloSFI = extSOF.mmsoloSFI;
mmsoloLineVs = extSOF.mmsoloLineVs;
mmsoloDs = extSOF.mmsoloDs;

i = 4
k = length(mmsoloSFI) - i;
%k = 748
mmsoloSFI(k, :)
mmsoloLineVs(k, :)
mmsoloDs(k, :)
geneID = (myGenes(mmsoloSFI(k, 1)))
geneSym = GOdata.geneSymbols(geneID)
 i = i + 1;

% how many of ATP stuff are in here 

% selecting one (also change the plot loop)
Ds = mmsoloDs;
Ds = soloDs;
sp = .5;

c1 = [200, 200, 200]/255
c2 = [80, 80, 80]/255
h = figure
set(h, 'Position', [100, 100, 450, 750])
hold all
for i = 1:length(Ds)
    x01 = sp;
    x02 = -sp
    
    % >> way one
    % x1 = (Ds(i, 1));
    % x2 = (Ds(i, 2));
    % x3 = Ds(i, 3);
    % x4 = Ds(i, 4);
    
    % way two % works with mmsoloDs
    [a, b] = max([Ds(i, [2, 4])]);
    [c, d] = min([Ds(i, [2, 4])]);
    x1 = (Ds(i, b*2-1));
    x2 = a;
    x3 = Ds(i, d*2-1);
    x4 = c;
    
    if (b == 1)
        maxf(i) = mmsoloSFI(i, 4);
        minf(i) = mmsoloSFI(i, 5);
        fixedLineVs(i, :) = mmsoloLineVs(i, :);
    else
        maxf(i) = mmsoloSFI(i, 5);
        minf(i) = mmsoloSFI(i, 4);
        fixedLineVs(i, 3:4) = mmsoloLineVs(i, 1:2);
        fixedLineVs(i, 1:2) = mmsoloLineVs(i, 3:4);
    end
    
    % line([x01, x1+sp], [i, i], 'LineWidth', 1,'Color', c1)
    % line([x01, x2+sp], [i, i], 'LineWidth', 1 ,'Color', c2)
    % line([x02, -x3-sp], [i, i], 'LineWidth', 1 ,'Color', c1)
    % line([x02, -x4-sp], [i, i], 'LineWidth', 1 ,'Color', c2)
end 
h = figure
plot(log2(maxf), '.')

h = figure
plot(log2(minf), '.')

thefCounts.minf = minf;
thefCounts.maxf = maxf;
save('~/resultsAndFigures/functionalIdentityForGenes/sofUtil_sortedMinMaxFs.mat', ...
     'thefCounts')

save('~/resultsAndFigures/functionalIdentityForGenes/sofUtil_sortedLineVs.mat', ...
     'fixedLineVs')

load('~/resultsAndFigures/functionalIdentityForGenes/sofUtil_sortedLineVs.mat')


c1 = [165, 189, 219]/256
c2 = [12, 44, 132]/256 
h = figure
sp = .5
set(h, 'Position', [100, 100, 450, 750])
hold all
for i = 1:length(fixedLineVs)
    x01 = sp;
    x02 = -sp
    
    % >> way one
    x1 = log2(fixedLineVs(i, 1));
    x2 = log2(fixedLineVs(i, 2));
    x3 = log2(fixedLineVs(i, 3));
    x4 = log2(fixedLineVs(i, 4));
    
    line([x01, x1+sp], [i, i], 'LineWidth', 1,'Color', c1)
    line([x01, x2+sp], [i, i], 'LineWidth', 1 ,'Color', c2)
    line([x02, -x3-sp], [i, i], 'LineWidth', 1 ,'Color', c1)
    line([x02, -x4-sp], [i, i], 'LineWidth', 1 ,'Color', c2)
end 

ylim([-90, 2800])
xlim([-10.5, 10.5])
set(gca, 'XTick', [ -10.5, -5.5, -.5, .5, 5.5, 10.5 ...
                   ])
set(gca, 'XTickLabel', {'10', '5', '0', '0', '5', '10' ...
                    })

figFolder = '~/resultsAndFigures/functionalIdentityForGenes/figures/'
set(h, 'PaperPositionMode', 'auto')
file = sprintf('%sSOFInstances_geneExpThings', figFolder);
print(h, '-depsc', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])

% examining the SOF cases for the overlap with go terms for the top
% 200 based on ND. 
% 8.6 just getting everything for local plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the last 200 extreme cases 
myGenes;
mmsoloSFI;
mmsoloLineVs;
mmsoloDs;

stepup = 2

c = 0
for stepup = 64:220
    j = 1
    stepup = sup(j);
    j = j + 1
    G = length(mmsoloSFI) - stepup;
    ids = myGenes(mmsoloSFI(:,1));
    syms = gpl570.uniqueSymbols(ids);
    syms(G)
    %whos syms;
    mmsoloSFI(G, :)

    g1 = ids(G)
    sum(fMat(g1, :))

    disID1 = mmsoloSFI(G,1);
    t = 2;
    book = inTerms(logical(reshape(tsfs(disID1, t, :), [1, 3507])));
    book{:}
    book = inTerms(logical(fMat(g1, :)));
    book{:}

    g1 = ids(G)
    % get the foverlap for fMat
    fs = fMat(g1, :);
    OS = zeros(5, 1);
    overlaps = zeros(3507, 5);
    for i = 1:5
        temp = reshape(tsfs(disID1, i, :), [1, 3507]);
        OS(i) = temp * fs';
        overlaps(:, i) = (temp + fs) == 2;
    end

    observe = overlaps' * overlaps
    syms(G)
    sum(fMat(g1, :)) 
    if (sum(sum(observe))>0)
        c = c + 1;
        sup(c) = stepup;
    end
end

% getting the network neighbours
g1Net = zeros(5, 18494);
for t = 1:5
    net = atn{t};
    g1Net(t, :) = net(g1, :) + net(:, g1)';
end

sum(g1Net')
g1Net * g1Net'

% getting the enriched functions 
[a, b] = ismember(g1, myGenes)
g1Fs = reshape(tsfs(b, :, :), [5, fCount]);

g1Fs * g1Fs'
sum(g1Net')
sum(g1Fs')

% >>>> bonus: getting the correlation 
myGenes;
mmsoloSFI;
mmsoloLineVs;
mmsoloDs;

% correlation of max for that. or this. 
ar = zeros(1, 3);
% ar(1) = corr(mmsoloSFI(:, 4), mmsoloDs(:,2)) % f count and tsnd
% ar(2) = corr(mmsoloSFI(:, 5), mmsoloDs(:,4)) % f count and tsnd

ar(1) = corr([sfInstances(:, 4); sfInstances(:, 5)], [theDs(:,2); theDs(:,4)]) % f count and tsnd

% ar(3) = corr(mmsoloDs(:, 2), mmsoloLineVs(:,2)) % tsnd and non exp
% ar(4) = corr(mmsoloDs(:, 4), mmsoloLineVs(:,4)) % tsnd and non exp

ar(2) = corr([theDs(:, 2); theDs(:, 4)], [lineVs(:,2); lineVs(:,4)]) % tsnd and non exp

% ar(5) = corr(mmsoloSFI(:, 4), mmsoloLineVs(:,2)) % f count non exp
% ar(6) = corr(mmsoloSFI(:, 5), mmsoloLineVs(:,4)) % f count non
% exp 

ar(3) = corr([sfInstances(:, 4); sfInstances(:, 5)],[lineVs(:, 2); ...
                    lineVs(:, 4)]); % f count non

save(['~/resultsAndFigures/functionalIdentityForGenes/' ...
      'FIcount_tsnd_nonExp_corrs_JCC0.1_3inst.mat'], 'ar')

save(['~/resultsAndFigures/functionalIdentityForGenes/' ...
      'FIcount_tsnd_nonExp_corrs_JCC0.01_3inst.mat'], 'ar')

ar1 = ar;
load(['~/resultsAndFigures/functionalIdentityForGenes/' ...
      'FIcount_tsnd_nonExp_corrs_JCC0.01_3inst.mat'])

% 9. examining the moderate JCC instances for the cases
% I can use almost all the above code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fMat;
myGenes;
bigJCCInstances;

ids = myGenes(bigJCCInstances(:,1));
syms = gpl570.uniqueSymbols(ids);

hits = zeros(length(bigJCCInstances), 1);
sumFs = zeros(length(bigJCCInstances), 1);
diagOnly = zeros(length(bigJCCInstances), 1);
for i = 1:length(bigJCCInstances)
    g = ids(i);
    sum(fMat(g, :))
    disID1 = bigJCCInstances(i,1);
    t = 2;
    book = inTerms(logical(reshape(tsfs(disID1, t, :), [1, 3507])));
    book{:};
    book = inTerms(logical(fMat(g1, :)));
    book{:};

    % get the foverlap for fMat
    fs = fMat(g, :);
    OS = zeros(5, 1);
    overlaps = zeros(3507, 5);
    for j = 1:5
        temp = reshape(tsfs(disID1, j, :), [1, 3507]);
        OS(j) = temp * fs';
        overlaps(:, j) = (temp + fs) == 2;
    end

    observe = overlaps' * overlaps;
    ll = logical(diag(observe));
    sumFs(i) = sum(sum(observe));
    hits(i) = sum(ll);
    diagOnly(i) = sum(diag(observe)) - sumFs(i) + sum(diag(observe));
end

finalHits = (diagOnly > 3) + (hits>1);

hist(hits)

kado = find(finalHits == 2);
i = 1
i = i + 1
bigJCCInstances(kado(i), 1:6)
sumFs(kado(i))

disID1 = bigJCCInstances(kado(i), 1)
g = myGenes(disID1);
fs = fMat(g, :);
OS = zeros(5, 1);
overlaps = zeros(3507, 5);
for j = 1:5
    temp = reshape(tsfs(disID1, j, :), [1, 3507]);
    OS(j) = temp * fs';
    overlaps(:, j) = (temp + fs) == 2;
end
observe = overlaps' * overlaps
gpl570.uniqueSymbols(g)

% 10. getting the tsn for all vs things : also, this is for JCC <
% 0.01, I should do it for JCC < .1 really and then it is done. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myGenes;
theDs;
sfInstances = bigJCCInstances;
atnds;
tsnds;

vgtsnds = tsnds(:, myGenes);
alltsnds = vgtsnds(:);

allSelected = ([theDs(:, 2); theDs(:, 4)]);

h = figure
b1 = hist(maxSelected, [0:50:1500]);%/length(maxSelected);

b2 = hist(maxtsnds, [0:50:1500]);%/length(maxtsnds); 

plot((log2(b1+1)/sum(log2(b1+1))))
hold on
plot(log2(b2+1)/sum(log2(b2+1)), 'r')

cdfplot(allSelected)
hold on
cdfplot(alltsnds)

set(gca, 'XTick', [0:200:1500], 'XTickLabel', [0:200:1500])

vgtsnds = tsnds(:, myGenes);
vgatnds = atnds(:, myGenes);
vgp = vgtsnds./(vgatnds+1);
allvp = vgp(:);
b1 = hist(allvp, [0:.1:1])/length(allvp)

dstsn = ([theDs(:, 2); theDs(:, 4)]);
dsatn = ([theDs(:, 1); theDs(:, 3)]);
selecteddp = dstsn./(dsatn);
b2 = hist(selecteddp, [0:.1:1])/length(selecteddp);

h = figure
bar([b1; b2]')
set(gca, 'XTick', [1:11], 'XTickLabel', [0:.1:1])

% this is ok. but I want to find the selected cases and remove them
% from the overlal. 

vgtsnds = tsnds(:, myGenes);
vgatnds = atnds(:, myGenes);

% marking the selected cases as -1 in vgtsnds and vgatnds
for i = 1:length(sfInstances)
    ind = sfInstances(i, 1);
    t1 = sfInstances(i, 2);
    t2 = sfInstances(i, 3);
    vgtsnds(t1, ind) = -1;
    vgtsnds(t2, ind) = -1;
    vgatnds(t1, ind) = -1;
    vgatnds(t2, ind) = -1;
end

sibtsnds = vgtsnds(:);
sibatnds = vgatnds(:);

book = sibtsnds > -1;
sibtsnds = sibtsnds(book);
sibatnds = sibatnds(book);

vgp = sibtsnds ./ (sibatnds+1);
b1 = hist(vgp, [0:.1:1])/length(vgp);
b2 = hist(selecteddp, [0:.1:1])/length(selecteddp);

h = figure
mybars = log2(([b1; b2]') * 10000);
bar(mybars)

sum(vgp > .70) / length(vgp)
sum(selecteddp > .70)/length(selecteddp)

hold on
plot(mybars(:,1))
hold on
plot(mybars(:,2))
set(gca, 'XTick', [1:11], 'XTickLabel', [0:.1:1])
legend('in all cases', 'in extreme SOF cases')
xlabel('portion of tsnd to atnd')
ylabel('dist scaled to 10000, log2 transferred')

figFolder = '~/resultsAndFigures/functionalIdentityForGenes/figures/'
file = sprintf('%sSOFInstances_tsatnds_JCC.1', figFolder);
print(h, '-depsc', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])

% 11. that tiny bar plot of count of genes expressed in multiple
% tissues 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expMat;
book = sum(expMat');
mybar = [sum(book ==2), sum(book == 3), su(book ==4), sum(book == ...
                                                  5)]
h = figure
mybar = [mybar; [0 0 0 0]]
bar(mybar, 'stacked')

figFolder = '~/resultsAndFigures/functionalIdentityForGenes/figures/'
file = sprintf('%scountOfGenesExp_multiTissue_bar', figFolder);
print(h, '-depsc', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])

% 12. The violin plot of the JCC and FI# for the pure, semi-pure and
% null
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% >>>>>>>>>>>>>>>> Load the three groups: 

% >>> GROUP 1. null 

clear
load('~/resultsAndFigures/tempResults/mainOverlapJK_null_nonComp_GO07.mat')

load('~/resultsAndFigures/tempResults/mainOverlap_null_nonComp_GO07.mat')
    
load('~/resultsAndFigures/tempResults/allOverlap_null_nonComp_GO07.mat')

load(['~/resultsAndFigures/tempResults/' ...
      'allmin_null_nonComp_GO07.mat'])

% >> How many links have > 5 FI in two tissues
passMin = allmin > 5;
sib = sum(passMin');
sum(sib > 1)
finalLinks = sib > 1;

% >>> GENERAL: get the base tissue and the other tissue. 
copyShareFunctionJCC = shareFunctionJK;
copyShareFunctionJCC(~passMin) = -1;
[maxJCC, maxInd] = max(copyShareFunctionJCC');

copyShareFunctionJCC = shareFunctionJK;
copyShareFunctionJCC(~passMin) = 1;
[minJCC, minInd] = min(copyShareFunctionJCC');

minFs = zeros(size(maxJCC));
for i = 1:length(minFs)
    minFs(i) = allmin(i, maxInd(i));
end

minFcounts = zeros(size(minInd));
for i = 1:length(minInd)
    minFcounts(i) = allmin(i, minInd(i));
end 

x = minFcounts(finalLinks); % plot 4
y = minFs(finalLinks); % plot 3
plot(x, y, '.')
corr(x', y') % they are correlated. The higher the first one, the
             % higher other one. How is it for the pure and semi
             % pure. This means that these are probably the group
             % of genes which do not switch function very much. 

w = maxJCC(finalLinks); % plot 1
z = (minJCC(finalLinks)); % plot 2
plot(w, z, '.')
corr(w', z') % they are correlated, highly highly. The more similar
             % they are in base tissue, 

null1 = w;
null2 = z;
null3 = y;
null4 = x;

% 2. pure
load('~/resultsAndFigures/tempResults/pureData_allThings.mat')

% >>> getting count of genes in the tissues 
% f1s = zeros(size(a));
% f2s = zeros(size(b));
% for i = 1:length(a)
%     ind = find(myGenes == a(i));
%     t = c(i);
%     f1s(i) = sum(tsfs(ind, t, :));
    
%     ind = find(myGenes == b(i));
%     t = c(i);
%     f2s(i) = sum(tsfs(ind, t, :));
% end 
% hist(f1s)
% book = min([f1s'; f2s']);
% h = figure
% hist(book)
% myb = hist(book);
% bar(log2(myb))

pureData
% >> How many links have > 5 FI in the pure tissue. 
book = pureData.links(:, [3,4]);
minPureFI = min(book');
hist(minPureFI)
sum(minPureFI > 5)
pureMinFIpass = minPureFI > 5;

% >> How many links have > 5 FI in one of the other tissues. 
book = pureData.allmin;
minPass = book > 5;
sib = sum(minPass');
sum(sib > 0)
pureMinOtherPass = sib > 0;

% >> How many links have both. 
kado = pureMinOtherPass + pureMinFIpass;
sum(kado == 2)
finalLinks = kado == 2;

copyJCothers = pureData.JCothers;
copyJCothers(~minPass) = 1;
[minJCothers, minInd] = min(copyJCothers');
minFcounts = zeros(size(minInd));
for i = 1:length(minInd)
    minFcounts(i) = pureData.allmin(i, minInd(i));
end 

% not correlated. It is 0.03
x = minFcounts(finalLinks);
y = minPureFI(finalLinks);
plot(x, y, '.')
corr(x', y') % the correlation is about zeros. that means that
             % there is no pattern, genes just tend to have more
             % functions annotated to them in the tissue that they
             % have the pure link in. same is for the null? 

% >>>>>> PLOT6: scatter of plots 1 and 2
z = pureData.sfJC(finalLinks);
w = (minJCothers(finalLinks));
plot(z, w, '.')
corr(z, w') % They are FAR less Correlated in this case too... it

pure1 = z;
pure2 = w;
pure3 = y;
pure4 = x;

% 3. semi pure
load('~/resultsAndFigures/tempResults/semiPureData_allThings.mat')
semiPureData

minPureFC = zeros(165844, 1); 
for i = 1:length(minPureFC)
    pti = semiPureData.links(i, 3);
    minPureFC(i) = min(semiPureData.links(i, (3 + pti)), semiPureData.links(i, (8 + pti)));
end

% getting the final plots for the semi pure
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% >>>>> NUMBERS:
semiPureData

% >> How many links have > 5 FI in the pure tissue. 
sum(minPureFC >  5)
pureMinFIpass = minPureFC > 5;

% >> How many links have > 5 FI in one of the other tissues. 
book = semiPureData.allmin;
minPass = book > 5;
sib = sum(minPass');
sum(sib > 0)
pureMinOtherPass = sib > 0;

% >> How many links have both. 
kado = pureMinOtherPass + pureMinFIpass';
sum(kado == 2)
finalLinks = kado == 2;

copyJCothers = semiPureData.JCothers;
copyJCothers(~minPass) = 1;
[minJCothers, minInd] = min(copyJCothers');

minFcounts = zeros(size(minInd));
for i = 1:length(minInd)
    minFcounts(i) = semiPureData.allmin(i, minInd(i));
end 

% not correlated. It is zero
x = minFcounts(finalLinks);
y = minPureFC(finalLinks);
plot(x, y, '.')
corr(x', y) % they are NOT correlated (just like pure). The higher the first one, the
             % higher other one. How is it for the pure and semi
             % pure ...

% >>>>>> PLOT6: scatter of plots 1 and 2
w = semiPureData.sfJC(finalLinks);
z = (minJCothers(finalLinks));
plot(x, y, '.')
corr(x, y') % They are FAR less Correlated in this case too... pure

semiPure1 = w;
semiPure2 = z;
semiPure3 = y;
semiPure4 = x;

% I want two violin plots, comparison of pairs for each case : 6
% violin in each: I need a three column violin files : the JCC and
% the FICounts 

hist(semiPure1)
hist(pure1)
hist(null1)

pl = length(pure1);
spl = length(semiPure1);
nl = length(null1);
fcLen = (pl + spl + nl) * 2;
fcVio = zeros(fcLen, 3); % 3 columns: values, base/other tissue,
                         % pure/semiPure/null links

fcVio(1:pl, 1) = pure1;
fcVio(1:pl, 2) = 1; % "base" tissue 
fcVio(1:pl, 3) = 1; % "pure" links

s = pl + 1
e = 2 * pl
fcVio(s:e, 1) = pure2;
fcVio(s:e, 2) = 2; % "other" tissue
fcVio(s:e, 3) = 1; % "pure" links

s = e + 1
e = s + spl - 1
fcVio(s:e, 1) = semiPure1;
fcVio(s:e, 2) = 1; % "base" tissue
fcVio(s:e, 3) = 2; % "semiPure" links

s = e + 1
e = s + spl - 1
fcVio(s:e, 1) = semiPure2;
fcVio(s:e, 2) = 2; % "other" tissue
fcVio(s:e, 3) = 2; % "semiPure" links

s = e + 1
e = s + nl - 1
fcVio(s:e, 1) = null1;
fcVio(s:e, 2) = 1; % "base" tissue
fcVio(s:e, 3) = 3; % "null" links

s = e + 1;
e = s + nl - 1;
fcVio(s:e, 1) = null2;
fcVio(s:e, 2) = 2; % "other" tissue
fcVio(s:e, 3) = 3; % "null" links

% save the file 

FIcVio = zeros(fcLen, 3); % 3 columns: values, base/other tissue,
                         % pure/semiPure/null links

FIcVio(1:pl, 1) = pure3;
FIcVio(1:pl, 2) = 1; % "base" tissue 
FIcVio(1:pl, 3) = 1; % "pure" links

s = pl + 1
e = 2 * pl
FIcVio(s:e, 1) = pure4;
FIcVio(s:e, 2) = 2; % "other" tissue
FIcVio(s:e, 3) = 1; % "pure" links

s = e + 1
e = s + spl - 1
FIcVio(s:e, 1) = semiPure3;
FIcVio(s:e, 2) = 1; % "base" tissue
FIcVio(s:e, 3) = 2; % "semiPure" links

s = e + 1
e = s + spl - 1
FIcVio(s:e, 1) = semiPure4;
FIcVio(s:e, 2) = 2; % "other" tissue
FIcVio(s:e, 3) = 2; % "semiPure" links

s = e + 1
e = s + nl - 1
FIcVio(s:e, 1) = null3;
FIcVio(s:e, 2) = 1; % "base" tissue
FIcVio(s:e, 3) = 3; % "null" links

s = e + 1;
e = s + nl - 1;
FIcVio(s:e, 1) = null4;
FIcVio(s:e, 2) = 2; % "other" tissue
FIcVio(s:e, 3) = 3; % "null" links

h = figure;
boxplot(fcVio(:, 1), [fcVio(:, 3) fcVio(:,2)])
figFolder = '~/resultsAndFigures/functionalIdentityPureLinks/figures/'
file = sprintf('%sFCs_boxplots', figFolder);
print(h, '-depsc', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])

h = figure;
boxplot(FIcVio(:, 1), [FIcVio(:, 3) FIcVio(:,2)])
figFolder = '~/resultsAndFigures/functionalIdentityPureLinks/figures/'
file = sprintf('%sFICs_boxplots', figFolder);
print(h, '-depsc', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])


fileName = 'pythonViolin_JCC.csv'
file = ['~/resultsAndFigures/functionalIdentityPureLinks/' fileName]
fid = fopen(file, 'w')

% printing the Jaccard Sim
fprintf(fid, ['data,', 'tissue,' 'link\n'])

for i = 1:fcLen
    fprintf(fid, ['%d,%d,%d\n'], fcVio(i, 1), fcVio(i, 2), fcVio(i, ...
                                                      3));
end
fclose(fid)

% printing the FICount 
fileName = 'pythonViolin_FIC.csv'
file = ['~/resultsAndFigures/functionalIdentityPureLinks/' fileName]
fid = fopen(file, 'w')

fprintf(fid, ['data,', 'tissue,' 'link\n'])

for i = 1:fcLen
    fprintf(fid, ['%d,%d,%d\n'], FIcVio(i, 1), FIcVio(i, 2), FIcVio(i, ...
                                                      3));
end
fclose(fid)

% 13. getting the heatmap of XYWZ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
load('~/resultsAndFigures/commonAndPure/pureXYWZmat.mat')
load('~/resultsAndFigures/commonAndPure/nullXYWZmat.mat')
load('~/resultsAndFigures/commonAndPure/semiXYWZmat.mat')

addpath('~/codes/MATLAB/myCodes/general')

h = figure; 
set(h, 'Position', [100, 100, 950, 300])
subplot(1, 3, 1)
sib = corr(nullXYWZmat, 'type', 'Spearman')
sib(2,2) = 0.011;
sib(1,1) = 1
HeatMap(sib)
heatmap(sib, [],  [],  [], 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'gray')
title('null')

subplot(1, 3, 2)
sib = corr(pureXYWZmat, 'type', 'Spearman')
sib(1,1) = 1
heatmap(sib, [],  [],  [], 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'gray')
title('pure')

subplot(1, 3, 3)
sib = corr(semiXYWZmat, 'type', 'Spearman')
sib(1,1) = 1
heatmap(sib, [],  [],  [], 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'gray')
title('semi')

set(h, 'PaperPositionMode', 'auto')
figFolder = '~/resultsAndFigures/functionalIdentityPureLinks/figures/'
file = sprintf('%sXYWZCorrs_heatmap', figFolder);
print(h, '-depsc', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])

% 14. Individual examples in the last two result parts (SOF and pure)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
load('~/resultsAndFigures/functionalIdentityForGenes/allLineVs.mat')
load('~/data/general/GOdata_GPL570_07_nonComp.mat')

load(['~/resultsAndFigures/functionalIdentityForGenes/' ...
      'plotData_SOF.mat'])

load(['~/resultsAndFigures/functionalIdentityForGenes/' ...
      'sofUtil_sortedLineVs.mat'])

load(['~/resultsAndFigures/functionalIdentityForGenes/' ...
      'lineVs(exp)_Ds_sorted.mat'])

load('~/resultsAndFigures/functionalIdentityForGenes/sofUtil.mat')
% load random file list 
file = '~/data/ShamsGeneLists/highconfSFARIgenes.txt'
file = '~/data/ShamsGeneLists/sanders.txt'
file = '~/data/ShamsGeneLists/iosssifov.txt'
fid = fopen(file)
myList = textscan(fid, ['%s', repmat('%d', 1, 4), '%s', '%d'], ...
                       'headerlines', 0, 'Delimiter', '\t');
myList = myList{1}

%load('~/resultsAndFigures/functionalIdentityForGenes/sofUtil_lineVs.mat')

geneSyms = GOdata.geneSymbols;

% >>>>>> Example 1: the SOF extreme cases. 

% Genes with low JCC and >50 both sides
jccSFInstances = plotData.jccSFInstances;
a = jccSFInstances(:, 7) < .1;
sum(a)
b = jccSFInstances(:, [4,5]) >  10;
sum(b)

sib = a' + sum(b');
sum(sib == 3)

selected = jccSFInstances((sib==3), :);
trueID = plotData.myGenes(selected(:, 1));

geneSyms(trueID)
mySyms = (unique(geneSyms(trueID)))

[a, b] = ismember(myList, mySyms);
sum(a)

ADSL = 392;
ADAR = 319;
kado1 = find(plotData.myGenes == 319)
kado2 = find(jccSFInstances(:, 1) == kado1)

jccSFInstances(kado2, :)

% terms in the tissues 

% Genes with low JCC and one side big one side small

% Genes with low JCC and both sides small..., but still
% interesting. 

% >>>>>> Example 2: the dark blue cases. 
load(['~/resultsAndFigures/functionalIdentityForGenes/' ...
      'lineVs(exp)_Ds_sorted.mat']) 
fixedLineVs;
mmsoloLineVs;

% >>>>>> Example 3: the top volcano. 

% >>>>>> Example 4: the pure things ...

% >>>>>> The TF check 

file = '~/data/TFTargetList/Cell_10.txt'
fid = fopen(file)
myList = textscan(fid, ['%s', repmat('%d', 1, 4), '%s', '%d'], ...
                       'headerlines', 0, 'Delimiter', '\t');
list1 = cell(17473, 1);
list2 = cell(17473, 1);
for i = 1:length(myList{1})
    sib = strsplit(myList{1}{i}, ' ');
    list1{i} = sib{1};
    list2{i} = sib{2};
end

length(unique(list1))
length(unique(list2))

[a, b] = ismember(mySyms, list1);
sum(a)
a1 = a;
[a, b] = ismember(mySyms, list2);
sum(a)
a2 = a;

toTestList = mySyms(a);
i = 7
toTestList(i)
[c, d] = ismember(toTestList(i), geneSyms)
geneID = d

% Now I want to see if the clusters these TSFs belong to, are
% sort of consistent. 
load('~/data/clusterings/clusterAccessories.mat')
load('~/data/clusterings/totalMat_FDR0.050_noMF.mat')
load('~/data/clusterings/allClusters')

myCInds = allClusters(d, :);
myTNames = cAccess.typeNames(logical(myCInds));
myTisNames = cAccess.tissueNames(logical(myCInds));

myCs = allClusters(:, logical(myCInds));

mcCount = size(myCs, 2)
cjccsim = zeros(mcCount, mcCount);
for i = 1:mcCount
    for j = (i+1):mcCount
        book = myCs(:, i) + myCs(:, j);
        cjccsim(i,j) = sum(book == 2) / sum(book >0);
    end
end

addpath('~/codes/MATLAB/myCodes/general')
heatmap(cjccsim, [],  [],  [], 'TickAngle', 45, 'ShowAllTicks', ...
        true, 'Colorbar', true, 'Colormap', 'gray')

[k1, k2] = find(cjccsim > .1);

hit  = [k1, k2];
size(hit)

% printing out the name of clusters which are close to each other,
% then finding out where the SOF happened 
pairs = cell(length(hit), 1);

for i = 1: length(hit)
    s = sprintf('%s\t%s\t%s\t%s', myTNames{hit(i, 1)}, ...
                myTNames{hit(i, 2)},  myTisNames{hit(i, 1)},myTisNames{hit(i, 2)});
    pairs{i} = s;
end

pairs

%>>>>>>> things to check: 

% for big SOF: the vlines, blue thing. 

% the TF partners? 
[x1, y1] = ismember(list1, 'HLF');
[x2, y2] = ismember(list2, 'HLF');

tlist1 = list2(x1);
tlist2 = list1(x2);

[x, y] = ismember(tlist1, geneSyms);
targets = y(x);
[x, y] = ismember(tlist2, geneSyms);
targets = [targets' y(x)']
geneSyms(targets)

% what are its targets enriched in? Are they its partners? Its TS
% partners? 

neighbours = zeros(5, 18494);

for i = 1:5
    tempnet = atn{i};
    neighbours(i, :) = tempnet(6923, :)' + tempnet(:, 6923);
    %    neighbours(i, targets) = neighbours(i, targets) + 2;
end

% >>>> about ATFs

[c, d] = ismember('ATF1', geneSyms)

myNet = zeros(10 ,10)
for i = 1:5
    tempNet = atn{i};
    myNet = myNet + tempNet(1136: 1145, 1136:1145);
    find(myNet)
end
 
% +++ some random plot and numbers, this is for poster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:5
    %    sum(sum(finalTable(i).noTSNet))
    %    sum(sum(finalTable(i).wholeNet))
     holu = finalTable(i).noTSNet + finalTable(i).cgNet;
    sum(sum(holu == 2))
end

h = figure
book = hist(plotData.jccSFInstances(:, end), 20) / length(plotData.jccSFInstances);
bar(book)

[xi, fi] = ksdensity(plotData.jccSFInstances(:,end), 'support', ...
                    'positive')

plot(fi, xi)

% for a given gene, get the neighbours, get the NBFI overlap in the
% tissues. 
[a, g1] = ismember('ATF2', geneSyms)

g1Net = zeros(5, 18494);
g1TS = zeros(5, 18494);
for t = 1:5
    net = atn{t};
    g1Net(t, :) = net(g1, :) + net(:, g1)';
    net = finalTable(t).wholeNet;
    g1TS(t, :) =  net(g1, :) + net(:, g1)';
end

sum(g1TS')

sum(g1Net')
book = g1Net * g1Net';
((sum(g1Net')/400) .* (50*50)).^.5

% JCC
JCCnet = zeros(5,5);
for i = 1:5
    for j = 1:5
        JCCnet(i,j) = book(i, j) /(book(i,i) + book(j, j) - book(i,j));
    end
end

% getting the enriched functions 
[a, b] = ismember(g1, SOFutil.myGenes)
g1Fs = reshape(SOFutil.tsfs(b, :, :), [5, fCount]);

book = g1Fs * g1Fs'
sum(g1Net')
sum(g1Fs')
((sum(g1Fs')/400) .* (50*50)).^.5

% JCC
JCCnet = zeros(5,5);
for i = 1:5
    for j = 1:5
        JCCnet(i,j) = book(i, j) /(book(i,i) + book(j, j) - book(i,j));
    end
end
JCCnet

jccSFInstances = plotData.jccSFInstances;

mySFin = find(jccSFInstances(:, 1) == b) ;

jccSFInstances(mySFin, :)

[a, b, c] = find(holu);

inds = a == 554;
geneSyms(b(inds))

b(inds)

[a, b] = ismember('ADAR', geneSyms)
g1 = b

[gar1, gar2, gar] = find(holu);
int2 = find(gar2 == g1);
int1 = find(gar1 == g1);

[a, b] = ismember('ABCF2', geneSyms)
g2 = b

g1Net = zeros(5, 18494);
g1TS = zeros(5, 18494);
for t = 1:5
    net = atn{t};
    g1Net(t, :) = net(g1, :) + net(:, g1)';
    g2Net(t, :) = net(g2, :) + net(:, g2)';
    net = finalTable(t).wholeNet;
    g1TS(t, :) =  net(g1, :) + net(:, g1)';
    g2TS(t, :) =  net(g2, :) + net(:, g2)';
end

[a, b] = ismember(g1, SOFutil.myGenes)
g1Fs = reshape(SOFutil.tsfs(b, :, :), [5, fCount]);

[a, b] = ismember(g2, SOFutil.myGenes)
g2Fs = reshape(SOFutil.tsfs(b, :, :), [5, fCount]);

book = g1Fs * g1Fs';
diag(book)
((diag(book)/) * (60*60)).^ .5
book = g2Fs * g2Fs';
diag(book)

book = g1Net * g1Net'
book = full(g2Net * g2Net')

book = g2Fs * g1Fs';
diag(book)


% 15. links that represent functional similarity in ATN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('~/resultsAndFigures/functionalIdentityForGenes/sofUtil.mat')
fMat = SOFutil.fMat;
tsfs = SOFutil.tsfs;
myGenes = SOFutil.myGenes;

% Tissues And final table
FC = '3'
FDR = '0012'
load(sprintf(['~/resultsAndFigures/TSlinks/finalTables/' ...
              'finalTable_CG13_FC%slog_FDR%s.mat'], FC, FDR))

atnSums = zeros(18494, 18494);
tsSum = zeros(18494, 18494);
for t = 1:5
    atnSums = atnSums + atn{t};
    tsSum = tsSum + (finalTable(t).wholeNet .* t);
end

sum(sum(atnSums > 0))
sum(sum(tsSum >0))

sib = fMat * fMat';
sd = diag(sib);
isd = diag(sd, 0);
sib = sib - isd;
sum(sum(sib > 0))
[a, b] = find(sib == 151);

binAtnSum = atnSums > 0;
binTsSum = tsSum > 0;
binAtnSum = binAtnSum - binTsSum;

fSimAtn = sib .* binAtnSum;
sum(sum(fSimAtn > 0))

sum(sum(fSimAtn > 0))/sum(sum(binAtnSum)) % about 3% of the ATN
                                          % links have some f similarity

[a, b, c] = find(fSimAtn);
h = figure
book = hist(c, [5, 10, 15, 20])/5;
bar(book)

fSimTsn = sib .* binTsSum;
sum(sum(fSimTsn > 0))

sum(sum(fSimTsn > 0))/sum(sum(binTsSum))
[a, b, c] = find(fSimTsn);
h = figure
book = hist(c, [5, 10, 15, 20]);
bar(book)

% 16. links that represent functional similarity in TSN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
% Tissues And final table
FC = '3'
FDR = '0012'
load(sprintf(['~/resultsAndFigures/TSlinks/finalTables/' ...
              'finalTable_CG13_FC%slog_FDR%s.mat'], FC, FDR))

% >>>> load the five ATN and TSN networks
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

% pairwise functional similarity of the genes
sib = fMat * fMat';

for t = 1:5
    net = atn{t} + atn{t}';
    kado = sib .* net;
    [a, b, c] = find(kado);
    length(a)
    h = figure;
    hist(c, 30)
    
    tsnet = finalTable(t).wholeNet + finalTable(t).wholeNet';
    kado = sib .* tsnet;
    [a, b, c] = find(kado);
    length(a)
    h = figure;
    hist(c, 30)

    %    hist(kado(kado > 0))

    expTemp = zeros(18494, 18494);
    expTemp(cGenes, cGenes) = 1;
    
    pureNet = expTemp .* finalTable(t).noTSNet;
    pureNet = (finalTable(t).cgNet + finalTable(t).noTSNet) == 2;
    kado = sib .* pureNet;

    [a, b, c] = find(kado);
    length(a)
    h = figure;
    hist(c, 30)
end

% from the pure links which represent functional similarity, I want
% to know if this similarity is represented through any other ATN
% or TSN or Pure link (any other ATN link?)
% so, basically, is it positive in any other tissue, that exact
% similarity? 
% 

pt = 2

NBFIfs1 = zeros(length(a), 5);
NBFIfs2 =  zeros(length(a), 5);

tsfsTotal1 = zeros(length(a), 5);
tsfsTotal2 = zeros(length(a), 5);

fsTotal1 = zeros(length(a), 5);
fsTotal2 = zeros(length(a), 5);

totalCommonFs = zeros(length(a), 1);
otherTissueFs1 = zeros(length(a), 5);
otherTissueFs2 = zeros(length(a), 5);

for i = 1:length(a)
    neighbourFs1 = zeros(5, 3507);
    neighbourFs2 = zeros(5, 3507);
    % for each link: 
    % 1. get the function(s) 
    g1 = a(i);
    g2 = b(i);
    f1s = fMat(g1, :);
    f2s = fMat(g2, :);
    
    fs = (f1s + f2s) == 2;
    totalCommonFs(i) = sum(fs);
    sum(fs)
    
    % 2. do they come up in the tissue or any other partner they
    % have in any other tissue: for each of them, in the other
    % tissues, do any of these terms come up? 
    
    % Checking NBFI: 
    % g1
    tg1 = find(SOFutil.myGenes == g1);
    tst1 = reshape(SOFutil.tsfs(tg1, :, :), 5, 3507);
    tsfsTotal1(i, :) = sum(tst1');
    NBFIfs1(i, :) = tst1 * fs' ;
    
    % g2
    tg2 = find(SOFutil.myGenes == g2);
    tst2 = reshape(SOFutil.tsfs(tg2, :, :), 5, 3507);
    tsfsTotal2(i, :) = sum(tst2');
    NBFIfs2(i, :) =     tst2 * fs' ;
    
    % partner:
    % get all the functions the partners have: does any other
    % partner have the TS function?
    g1nfs = zeros(5, 3507);
    g2nfs = zeros(5, 3507);
    
    g1ns = zeros(5, 18494);
    g2ns = zeros(5, 18494);
    for t = 1:5
        book = atn{t} + atn{t}';
        g1ns(t, :) = book(g1, :);
        g2ns(t, :) = book(g2, :);
        
        neighbourFs1(t, :) = sum(fMat(logical(g1ns(t, :)), :));
        neighbourFs2(t, :) = sum(fMat(logical(g2ns(t, :)), :));
    end
    
    neighbourFs1(pt, :) = neighbourFs1(pt, :) - fs;
    neighbourFs2(pt, :) = neighbourFs2(pt, :) - fs;
    
    % these sum1 and sum2 will show me that for gene 1 and two, how
    % many of the fs functions are recognized in each tissue. 
    sum1 = ((neighbourFs1' > 0)+0)' * fs';
    sum2 = ((neighbourFs2' > 0)+0)' * fs';
    
    % I want to see, which links, have no other representation of
    % these functions in the other tissues. first. 
    
    otherTissueFs1(i, :) = sum1;
    otherTissueFs2(i, :) = sum2;
    
    % totFunction presented 
    sum1 = ((neighbourFs1' > 0)+0)' * f1s';
    sum2 = ((neighbourFs2' > 0)+0)' * f2s';
    
    fsTotal1(i, :) = sum1;
    fsTotal1(i, :) = sum2;
end

% the result matrix: 
% G1: did we find "howmany" of the "pure shared" ( sum(fs) )
% functions in the NBFI of the same tissue, NBFI of other tissues, 
% the neighbours of same tissue, the neighbours of other tissues? -
% do the same for G2. 8 columns result, and these are for the pure
% links for which the both gene partners DO share some similar
% functions. 

resultMat1 = zeros(length(a), 4);
resultMat1(:, 1) = NBFIfs1(:, 1);
resultMat1(:, 2) = sum(NBFIfs1(:, 2:end)');
resultMat1(:, 3) = otherTissueFs1(:, 1);
resultMat1(:, 4) = sum(otherTissueFs1(:, 2:end)');
find(sum(resultMat1') == 0)

resultMat2 = zeros(length(a), 4);
resultMat2(:, 1) = NBFIfs2(:, 1);
resultMat2(:, 2) = sum(NBFIfs2(:, 2:end)');
resultMat2(:, 3) = otherTissueFs2(:, 1);
resultMat2(:, 4) = sum(otherTissueFs2(:, 2:end)');
find(sum(resultMat2') == 0)

% the count of function, from the common funtions, that came up
% somehow in the other tissues
book1 = sum(otherTissueFs1(:, 2:end)');
book2 = sum(otherTissueFs2(:, 2:end)');
sum(book1 == 0)
sum(book2 == 0)

find(book1 == 0)
find(book2 == 0)
k = 405
g1 = a(k)
g2 = b(k)
f1s = fMat(g1, :);
f2s = fMat(g2, :);
fs = (f1s + f2s) == 2;
sum(f1s)
sum(f2s)
sum(fs)
[a(k) b(k)] 

% getting the function
myTerms = inTerms(fs)
f1ts = inTerms(logical(fMat(a(k), :) - fs)); 
f2ts = inTerms(logical(fMat(b(k), :) - fs)); 
% how much they generally share in other tissues. 

g1 = g4
% getting their GTEx rep.
addpath('~/codes/MATLAB/myCodes/general')
finalMat = zeros(1, 17);
for t = 1:17
    fName = GTExNet(t).tissue;
    load(['~/data/GTEx/matFormat_GPL570/' fName '.mat'])
    dsGT = dataSet.mat;
    size(dsGT)
    smallDS = dsGT([g1, g2], :);
    sib = corr(smallDS');
    
    myCorrs = GTExNet(t).corrQs;
    
    corrQMat = 500 .* ones(2, 2);
    for i = 1:(2)
        for j = (i + 1):(2)
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
finalMat
% what are the functions that a pair of gens share, and why only in
% blood? 

g1fs = fMat(14361, :);
g2fs = fMat(17935, :);

% 17. pure links that represenet functional similarity 

sib = finalTable(1).noTSNet - pureNet;

[a, b, c] = find(sib);

% now the network of their common genes in brain ... 

%18. functional implication of different groups (networks) of
%links: I find how well functions are represented in the networks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
cGenes = sum(expMat') == 5;

load('~/data/newPPIN_18494Net.mat');
load('~/resultsAndFigures/functionalIdentityForGenes/sofUtil.mat')
fMat = SOFutil.fMat;
whos fMat
sib = fMat * fMat';
netFunctionOverlap_func(newPPI, sib)
sib = newPPI;

netFunctionOverlap_func(book, newPPI)

atnRs = zeros(1, 6);
for t = 1:5
    %  temp = netFunctionOverlap_func(atn{t}, sib);
    temp = netFunctionOverlap_func_halfNet(atn{t}, sib);
    atnRs(t) = temp(1,1);
end
temp = netFunctionOverlap_func(atnSums, sib);
atnRs(6) = temp;

tsnRs = zeros(1, 6);
for t = 1:5
    %  temp = netFunctionOverlap_func(finalTable(t).wholeNet, sib);
    temp = netFunctionOverlap_func_halfNet(finalTable(t).wholeNet, sib);
    tsnRs(t) = temp(1,1);
end
temp = netFunctionOverlap_func(tsSum, sib);
tsnRs(6) = temp;

GTExRs = zeros(1, 9);
dsIDs = [5, 6, 7, 8, 14, 15, 16, 17];
dsIDs = [ 6,14, 15, 16, 17];
sumGTEx = zeros(18494, 18494);
for t = 1:5
    %  temp = netFunctionOverlap_func(GTExNet(dsIDs(t)).net, sib);
    temp = netFunctionOverlap_func_halfNet(GTExNet(dsIDs(t)).net, sib);
  GTExRs(t) = temp(1,1);
  sumGTEx = sumGTEx + GTExNet(dsIDs(t)).net;
end
temp = netFunctionOverlap_func(sumGTEx, sib);
GTExRs(9) = temp;

expTemp = zeros(18494, 18494);
expTemp(cGenes, cGenes) = 1;

pureRs = zeros(1, 6);
pureSum = zeros(18494, 18494);
for t  = 1:5
    pureNet = expTemp .* finalTable(t).noTSNet;
    pureNet = (finalTable(t).cgNet + finalTable(t).noTSNet) == 2;
    %    temp = netFunctionOverlap_func(pureNet, sib);
        temp = netFunctionOverlap_func_halfNet(pureNet, sib);
    pureRs(t) = temp(1,1);
    pureSum = pureSum + pureNet;
end
temp = netFunctionOverlap_func(pureSum, sib);
pureRs(6) = temp;

funcInR.pure = pureRs;
funcInR.ts = tsnRs;
funcInR.atn = atnRs;
funcInR.GTEx = GTExRs;

PPINR.pure = pureRs;
PPINR.ts = tsnRs;
PPINR.atn = atnRs;
PPINR.GTEx = GTExRs;

save('~/resultsAndFigures/PPINIndicationOfTheLinks.mat', 'PPINR')
save('~/resultsAndFigures/FunctionalIndicationOfTheLinks.mat', ...
     'funcInR')

% how does the ratio grow for the atnSum 
rs = zeros(1, 5); 
for i = 1:5
    book = sumAtn == i;
    rs(i) = netFunctionOverlap_func_halfNet(book, sib);
end
save('~/resultsAndFigures/PPINIndicationOfTheLinks_sumAtn.mat', ...
     'rs')

rs = zeros(1, 5); 
for i = 1:5
    book = sumGTEx == i;
    rs(i) = netFunctionOverlap_func_halfNet(book, sib);
end
 save('~/resultsAndFigures/PPINIndicationOfTheLinks_sumGTEx.mat', ...
     'rs')
   
% 19. Which functions are represented in WHICH networks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

% get fMat, fCount
% Tissues And final table

tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}

expMat = zeros(18494, 5);
for t = 1:5
    tissue = tissues{t};
    load(['~/data/general/tissueExpGenes/' tissue ...
          'ExpGenes0.8.mat'])
    expMat(:, t) = expGenesInd;
end


load('~/resultsAndFigures/firstProject/functionalIdentityForGenes/sofUtil.mat')
fMat = SOFutil.fMat;
clear SOFutil;
fCount = size(fMat, 2);

FC = '3'
FDR = '0012'
load(sprintf(['~/resultsAndFigures/firstProject/TSlinks/finalTables/' ...
              'finalTable_CG13_FC%slog_FDR%s.mat'], FC, FDR))

% >>>> load the five ATN and TSN networks
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

tsnLC = sum(tsnds')./2;

fCount = size(fMat, 2);
gc = zeros(5, fCount);
bgExp = zeros(5, fCount);
obCount = zeros(5, fCount);
LRs = zeros(5, fCount);
hgps = zeros(5, fCount);
binps = zeros(5, fCount);
netDensities = [.007, .020, .008, .012, .006];
netDensities = [.0004, .004, .0011, .0005,.0004];

% each of this is one round of sampling...

for t = 1:5
    %myNet = finalTable(t).wholeNet;
    myNet = atn{t};

    for f = 1:fCount
        myTerm = fMat(:, f);
        sum(myTerm)

        % how many genes are marked as expressed in this term?
        eGenesCount = sum((myTerm + expMat(:, t)) == 2)
        eGenes = ((myTerm + expMat(:, t)) == 2);
        gc(t, f) = eGenesCount;
        % for ppin >>>>>>>>>>>
        % eGenesCount = sum((myTerm + expArr') == 2);
        % eGenes = ((myTerm + expArr') == 2);
        % <<<<<<<<<<<<<<<<<<<<

        % 81 genes with the term, 50 genes expressed. 
        % if this FSim was fully represented in this network, it would have
        % all the 50*49/2 links 

        % now let's see how many links are there, and how is it bigger than
        % the background expectations. 
        % the background expectation for the brain network is 0.02 

        if(eGenesCount >= 20)

            bgExp(t, f) = (eGenesCount * (eGenesCount-1)/2) * netDensities(t);
            litNet = myNet(eGenes, eGenes);
            obCount(t, f) = sum(sum(litNet));
            LRs(t, f) = (sum(sum(litNet))/((eGenesCount * (eGenesCount-1)/2) ...
                                           * netDensities(t) + 1));
            
            N = (eGenesCount * (eGenesCount-1)/2);
            k = N * netDensities(t);
            x = sum(sum(litNet));
            popSize = size(myNet, 1) * (size(myNet,1)-1)/2;
            hgps(t, f) = 1 - hygecdf(x, popSize, sum(sum(myNet)),N);
            %            binps(t, f) = 1 - binocdf(x, N, netDensities(t));
        else 
            bgExp(t, f) = -1;
            obCount(t, f) = -1;
            LRs(t, f) = -1;
            %            binps(t, f) = -1;
            hgps(t, f) = -1;
        end
        % so this function is not adding to that expectation, but which one
        % is? 
    end
end

% selected Fs for each of the tissues: 
sigTerms = zeros(5, fCount);
for t = 1:5
    selected = hgps(t, :) > -1;
    sLRs = LRs(t, selected);
    shgps = hgps(t, selected);
    sterms = find(selected);

    % correcting pvalues
    shgps = shgps * sum(selected);
    logps = log10(shgps);
    logps(logps < -10) = -10;
    logps = logps .* -1;
    %    plot(logRs, logps, '.')
    temp1 = (logps >= 2);
    sigTerms(t, sterms(temp1)) = 1;
end

% which one is missed why. 
t = 4
% the tsRes and other ..Res could be loaded, look down 
orSigTerms = tsRes.sigTerms(t, :);
tsSigTerms = sigTerms(t, :);
pureSigTerms = sigTerms(t, :);

mat = tsRes.sigTerms; 
t = 4;
sumMat = sum(mat);
kado = orSigTerms .* sumMat;
forSigTerms = kado == 1;
sum(forSigTerms)
% this gives me the terms whith are exclusive among the TSNs. If I
% want the terms exclusive between the TANs, I should also look at TANS
load('~/resultsAndFigures/functionsRepresented/taRes.mat')
sib = sum(taRes.sigTerms) - taRes.sigTerms(t, :);
kado = forSigTerms .* ~(logical(sib));
forSigTerms = kado;

kado = tsSigTerms .* sumMat;
ftsSigTerms = kado == 1;
sum(ftsSigTerms)

kado = pureSigTerms .* sumMat;
fpureSigTerms = kado == 1;
sum(fpureSigTerms)

book = ftsSigTerms .* fpureSigTerms;
sum(book)  
common = inTerms(book == 1);

% which functions are pureloss specific ...
purelossS = forSigTerms - fpureSigTerms;
fspu = inTerms(purelossS == 1);
fspu{:}

% which functions are TS loss specific ... 
tslossS = forSigTerms - ftsSigTerms;
fsts = inTerms(tslossS == 1);
fsts{:}

% >>>>>>>>>>>>>>>>>> Here I am getting the info for the pure and TS
% removals. 
% there are in total, 211 terms. 159 are still presenet with ts
% removal, 151 are still presenet with pure removal. pure removal
% causes less loss of links compared to the ts removal (7k vs 15k
% of links), but my guess is that the genes are focused. The
% question is, how different are these functions from each other? 
% let's see how many links each function lost, and what portion for
% each of the removals? 
lcorg = zeros(1, 211); % link counts...
lctsMod = zeros(1, 211); % link counts...
lcpureMod = zeros(1, 211); % link count ...

% for these 211 functions:
orgFInds = find(forSigTerms);
tsNet = finalTable(t).wholeNet;
tsMod = modNet;
pureMod = modNet;
subNets = zeros(211, 211);
geneOverlaps = fMat' * fMat;
miniGeneOverlaps = geneOverlaps(orgFInds, orgFInds);
parents = zeros(211, 211);
for f = 1: length(orgFInds)
    
    inGenes = logical(fMat(:, orgFInds(f)));
    % how many links it had in the ts network : whole net
    lcorg(f) = sum(sum(tsNet(inGenes, inGenes)));
    
    % how many links it have in the tsMod network
    lctsMod(f) = sum(sum(tsMod(inGenes, inGenes)));
    
    % how many links it have in the pureMod network
    lcpureMod(f) = sum(sum(pureMod(inGenes, inGenes)));
    
    % Also, whos is sub of who in the functions? get all the
    % parents of nes' terms , that is: fMat * fMat == diag()
    
    hisParents = find(miniGeneOverlaps(f, :) == miniGeneOverlaps(f, ...
                                                      f));
    parents(f, hisParents) = hisParents;
        
end

t
lungPureMod.bgExp = bgExp(t, :);
lungPureMod.obCount = obCount(t, :);
lungPureMod.LRs = LRs(t, :);
lungPureMod.hgps = hgps(t, :);
lungPureMod.binps = binps(t, :);
lungPureMod.sigTerms = sigTerms(t, :);
lungPureMod.terms = inTerms;
lungPureMod.termIDs = inGOIDs;
save('~/resultsAndFigures/functionsRepresented/lung_pureMod_Res.mat', ...
     'lungPureMod')
 
save('~/resultsAndFigures/functionsRepresented/brain_pureMod_Res.mat', ...
     'brainPureMod')

t
liverTsMod.bgExp = bgExp(t, :);
liverTsMod.obCount = obCount(t, :);
liverTsMod.LRs = LRs(t, :);
liverTsMod.hgps = hgps(t, :);
liverTsMod.binps = binps(t, :);
liverTsMod.sigTerms = sigTerms(t, :);
liverTsMod.terms = inTerms;
liverTsMod.termIDs = inGOIDs;
save('~/resultsAndFigures/functionsRepresented/liver_tsMod_Res.mat', ...
     'liverTsMod')

save('~/resultsAndFigures/functionsRepresented/brain_tsMod_Res.mat', 'brainTsMod')
 
ppiRes.bgExp = bgExp(1, :);
ppiRes.obCount = obCount(1, :);
ppiRes.LRs = LRs(1, :);
ppiRes.hgps = hgps(1, :);
ppiRes.binps = binps(1, :);
ppiRes.sigTerms = sigTerms(1, :);
ppiRes.terms = inTerms;
ppiRes.termIDs = inGOIDs;
save('~/resultsAndFigures/functionsRepresented/ppiRes.mat', 'ppiRes')
load('~/resultsAndFigures/functionsRepresented/ppiRes.mat')

% and who is sub port of who. 
% just saving it for the TS and TA networks.
tsRes.gc = gc;
tsRes.bgExp = bgExp;
tsRes.obCount = obCount;
tsRes.LRs = LRs;
tsRes.hgps = hgps;
%tsRes.binps = binps;
tsRes.sigTerms = sigTerms;
tsRes.terms = tsRes.terms; %inTerms;
tsRes.termIDs = tsRes.termIDs; %inGOIDs;

save('~/resultsAndFigures/functionsRepresented/tsRes.mat', ...
     'tsRes')
load('~/resultsAndFigures/firstProject/functionsRepresented/tsRes.mat')
% the second one has the genecount
save('~/resultsAndFigures/firstProject/functionsRepresented/tsRes02_HG.mat', ...
     'tsRes')

load('~/resultsAndFigures/firstProject/functionsRepresented/tsRes02.mat')

taRes.bgExp = bgExp;
taRes.obCount = obCount;
taRes.LRs = LRs;
taRes.hgps = hgps;
taRes.binps = binps;
taRes.sigTerms = sigTerms;
taRes.terms = tsRes.terms;
taRes.termIDs = tsRes.termIDs;

save('~/resultsAndFigures/firstProject/functionsRepresented/taRes_HG.mat', 'taRes')
load('~/resultsAndFigures/firstProject/functionsRepresented/taRes.mat')

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% I have a few things to investigate: 
% 1. average expresssion level of the genes for the TSly
% represented functions in other tissues, and in the tissue of
% study. OK. Get me the Exp file. 
load('~/data/general/linkExprInfo/dataSetExpInf.mat') % load
                                                      % dataSetExpInf
totalExp = zeros(18494, 53);
for d = 1:53
    totalExp(:, d) = dataSetExpInf(d).meanExp;
end
qnTotalExp = quantilenorm(totalExp);

t = 4;
mat = tsRes.sigTerms; 
mat = sigTerms;
% finding the functions specific to the tissue
sumMat = sum(mat);
kado = mat(t, :) .* sumMat;
sum(kado == 1)

tissueSFs = find(kado == 1);
whos tissueSFs

% finding the genes with have any of the functions
inGenes = logical(sum(fMat(:, tissueSFs)')); 
inGenesInds = find(inGenes);
sum(inGenes)

% getting the expression levels for those genes involved in the functions
avgExp = qnTotalExp(inGenes, :);

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% for genes who have a function:

% @@@ Second round: mod net
litNet = modNet(inGenes, inGenes);

litNet = finalTable(t).wholeNet(inGenes, inGenes);

fullNet = litNet + litNet';
nds = sum(fullNet);
%heatmap(fullNet)
ins = (nds > 0);
% recording the genes:
inGenesInds2 = inGenesInds(ins);

%colormap(gray)
smallnds = nds(ins);
%[a, b] = sort(smallnds, 'descend');
smallexp = avgExp(ins, :);
%smallexpSort = smallexp(b, :);
% h = figure;
% colormap(gray)
% heatmap(smallexpSort(1:50, :))
% heatmap(avgExp(ins, :)) 
% h = figure
% colormap(gray)
% heatmap(litNet(ins, ins)) 
heatmap(smallexp)
colormap(gray)

% lung
thisT = [10:24]
otherT = [1:9, 25:53];

% liver
thisT = [32:41]
otherT = [1:31, 42:53];

% brain
thisT = [42:53]
otherT = [1:41];

hs = zeros(length(inGenesInds2), 1);
ps = zeros(length(inGenesInds2), 1);
for i = 1: length(inGenesInds2)
    [h, p] = ttest2(smallexp(i, thisT), smallexp(i, otherT));
    ps(i) = p;
    hs(i) = h;
end

hist(ps)
ps = ps * length(inGenesInds2);
h = figure
hist(ps)
m1 = median(smallexp(:, otherT)');
m2 = median(smallexp(:, thisT)');
fch = m2./m1;

[a, b] = sort(ps);
thr = max(find(a <= .05))
sortsmallexp = smallexp(b, :);
sortsmallnds = smallnds(b);
sortfch = fch(b);
inGenesInds3 = inGenesInds2(b);

passfch = sortfch >= 2;
book = sortsmallnds .* passfch;
sum(book(1:thr) > 0)
sum(book(1:thr))

inGenesInds4 = inGenesInds3(passfch(1:thr));

h = figure
colormap(gray)
heatmap(qnTotalExp(inGenesInds4, :));
heatmap(sortsmallexp(1500:end, :)) 

inGenesInds4 = inGenesInds3(1158:end);

% if we remove these links, how many functions do we lose? 

% for a given netork, remove the links from inGenesInds4. 
t
modNet = finalTable(t).wholeNet;
sum(sum(modNet))
for i = 1:length(inGenesInds4)
    g = inGenesInds4(i);
    modNet(g, :) = 0;
    modNet(:, g) = 0;
end
sum(sum(modNet))

% for the given network, remove pure links
t
modNet = finalTable(t).wholeNet;
modNet = modNet - pureNet;
sum(sum(modNet))

lit2 = modNet(inGenes, inGenes);
sum(sum(lit2))

% now let's remove pure genes and their links from the whole
% net... 
% >> get the pure net from this

expTemp = zeros(18494, 18494);
expTemp(cGenes, cGenes) = 1;

pureNet = (finalTable(t).cgNet + finalTable(t).noTSNet) == 2;

litPure = pureNet(inGenes, inGenes);

% <<<

% getting the count of TS links between genes marked as expressed
% only in one tissue

sumExp = sum(expMat');
sum(sumExp == 1)
sum(sumExp >= 1)

sumTSNet = zeros(18494, 18494);
for t = 1:5
    sumTSNet = sumTSNet + finalTable(t).wholeNet;
end

% getting the count of tissues the paird genes are expressed in: 
tcountMat = zeros(18494, 18494);
expGenes = find(sumExp >= 1);
for i = 1 : length(expGenes)
    for j = i+1:length(expGenes)
        iec = sumExp(expGenes(i));
        jec = sumExp(expGenes(j));
        tcountMat(expGenes(i), expGenes(j)) = min(iec, jec);
    end
end


sum(sum(tcountMat == 1))
book = sumTSNet .* tcountMat;

for i = 1:5
    sum(sum(book == i)) / sum(sum(tcountMat == i))
end

% comparing the two sets of functions
% >>>>>>>>>>>>>>>>>>>>>>>>>>>

whos orSigTerms
whos tsSigTerms
whos pureSigTerms

brainLists.orSigTerms = orSigTerms;
brainLists.tsSigTerms = tsSigTerms;
brainLists.pureSigTerms = pureSigTerms;

save('~/resultsAndFigures/functionsRepresented/brainSigTerms.mat', ...
     'brainLists')

load(['~/resultsAndFigures/functionsRepresented/' ...
      'brain_pureMod_Res.mat'])

load(['~/resultsAndFigures/functionsRepresented/' ...
      'brain_tsMod_Res.mat'])

load('~/resultsAndFigures/functionsRepresented/tsRes02.mat')

% I want the TS functions which are not enriched in any other TAN
% either
t = 4; 
taMat = taRes.sigTerms;
tsMat = tsRes.sigTerms;

otherSumTAMat = sum(taMat) - taMat(t, :);
otherSumTSMat = sum(tsMat) - tsMat(t, :);

thisTerms = tsMat(t, :);

sib = (thisTerms .* -10) + otherSumTSMat + otherSumTAMat;
sum(sib == -10)

tissueTSfilter = sib == -10;

% brainTsMod.bgExp = bgExp(2, :);
% brainTsMod.obCount = obCount(2, :);
% brainTsMod.LRs = LRs(2, :);
% brainTsMod.hgps = hgps(2, :);
% brainTsMod.binps = binps(2, :);
% brainTsMod.sigTerms = sigTerms(2, :);
% brainTsMod.terms = inTerms;
% brainTsMod.termIDs = inGOIDs;

% OK, now that they are loaded, make me two lists. For each of the
% pure and and ts, give me the list of functions and how many edges
% were missing from the original one. 

% I want a matrix: 
% for each of the original 211 functions: 
% 0. the ID
% 1. the name, 
% 2. count of genes,
% 3. count of observed links in the TS 
% 4. count of observed links in the pureMod 
% 5. did we lose it in the pureMod
% 6. count of observed links in the TSgeneMod
% 7. did we lose it in the TSgeneMod

fnames = tsRes.terms(tissueTSfilter);
fIDs = tsRes.termIDs(tissueTSfilter);
gc = tsRes.gc(t, tissueTSfilter);
obsL = tsRes.obCount(t, tissueTSfilter);
pureObsL = lungPureMod.obCount(tissueTSfilter);
pureSig = lungPureMod.sigTerms(tissueTSfilter);
tsObsL = lungTsMod.obCount(tissueTSfilter);
tsSig = lungTsMod.sigTerms(tissueTSfilter);
pureObratio = pureObsL ./obsL;
tsObratio = tsObsL ./obsL;
joined =  [pureObratio' tsObratio'];
colormap(gray)
heatmap(joined) 

% the 211 are those with some TAN overlap. these are those with no
% TAN overlap
thelung91.fnames = fnames;
thelung91.fIDs = fIDs;
thelung91.gc = gc;
thelung91.obsL = obsL;
thelung91.pureObsl = pureObsL;
thelung91.pureObratio = pureObratio;
thelung91.pureSig = pureSig;
thelung91.tsObsL = tsObsL;
thelung91.tsObratio = tsObratio;
thelung91.tsSig = tsSig;
save('~/resultsAndFigures/functionsRepresented/thelung91.mat', ...
     'thelung91')

% the 211 are those with some TAN overlap. these are those with no
% TAN overlap
theliver43.fnames = fnames;
theliver43.fIDs = fIDs;
theliver43.gc = gc;
theliver43.obsL = obsL;
theliver43.pureObsl = pureObsL;
theliver43.pureObratio = pureObratio;
theliver43.pureSig = pureSig;
theliver43.tsObsL = tsObsL;
theliver43.tsObratio = tsObratio;
theliver43.tsSig = tsSig;
save('~/resultsAndFigures/functionsRepresented/theliver43.mat', ...
     'theliver43')

% the 211 are those with some TAN overlap. these are those with no
% TAN overlap
thebrain146.fnames = fnames;
thebrain146.fIDs = fIDs;
thebrain146.gc = gc;
thebrain146.obsL = obsL;
thebrain146.pureObsl = pureObsL;
thebrain146.pureObratio = pureObratio;
thebrain146.pureSig = pureSig;
thebrain146.tsObsL = tsObsL;
thebrain146.tsObratio = tsObratio;
thebrain146.tsSig = tsSig;
save('~/resultsAndFigures/functionsRepresented/thebrain146.mat', ...
     'thebrain146')

load('~/resultsAndFigures/functionsRepresented/thelung91.mat')
file = ['~/resultsAndFigures/functionsRepresented/thelung91_table.txt']
fid = fopen(file, 'w')
fprintf(fid, 'functionID\tfunctionName\tgeneCount\tlinkCount(lc)\tlc_minusPure\tlc_minusTsGenes\tpureRatio\ttsRatio\tPuresig\ttsSig\n')
for i = 1:91
    fprintf(fid, ['%d\t%s\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%d\t%d\n'], ...
            thelung91.fIDs(i),...
            thelung91.fnames{i},...
            thelung91.gc(i),...
            thelung91.obsL(i),...
            thelung91.pureObsl(i),...
            thelung91.tsObsL(i),...
            thelung91.pureObratio(i),...
            thelung91.tsObratio(i),...
            thelung91.pureSig(i),...
            thelung91.tsSig(i));
end
fclose(fid)

% basically, tell me which genes these are. remove them. what
% happens for the ts functions? 

% OK. I learned that the expression thing is true. Now I have two
% things: 1. how many links are between geens which are not
% expressed in other tissues : I got this. 
% 2. what about the TS functions which are shared? what is
% happening there : these are basically ts expressions as well, for
% the most parts. 
% 3. if we remove the links in 1, how many functions will we lose? 

[a, b, c] = find(fullNet);

% SCRAP
%%%%%%%%%%%%%%

% this bleongs to the top part, but I want to see how many
% functions have same group of genes
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
kado = fMat' * fMat;
families = zeros(size(kado));

famCount = zeros(1, length(kado));

observed = zeros(1, length(kado));
for i = 1:length(kado)
    if (observed(i) == 0)
        groups = kado(i, :);
        diff = groups - kado(i, i);
        sames1 = (diff == 0);
        D = diag(kado);
        sames2 = ((D - kado(i, i)) == 0);
        sames = (sames1 + sames2') == 2;
        families(sames, i) = 1;
        families(i, sames) = 1;
        observed(sames) = i;
    end
end

% groups of similar functions 
% get one from each. 

nonRinTermsInds = unique(observed);
nonRinTerms = zeros(1, length(kado));
nonRinTerms(nonRinTermsInds) =  1;
nonRinTerms = logical(nonRinTerms);
% not much overlap ... 

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

book = sib .* atn{5};
[a, b, c] = find(book);
boxplot(c);

h = figure
book = sib .* finalTable(5).wholeNet;
[a, b, c] = find(book);
boxplot(c);
ylim([-2, 70])

% which terms are better represented in the network then? 
% for each term, I have some genes. They could all be
% linked together, or none. 

% 20. The network numbers: how the TS link portion decrease 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sumNet = zeros(18494, 18494);
for t = 1:5
    %   tnet = atn{t};
    tnet = finalTable(t).wholeNet;
    sumNet = tnet + sumNet;
end
sum(sum(sumNet>0))

for t = 1:5
    expFilter = zeros(18494, 18494);
    inGenes = expMat(:,t);
    %tNet = sumNet - atn{t};
    tNet = sumNet - finalTable(t).wholeNet;
    expFilter(logical(inGenes), logical(inGenes)) = 1;

    sib = tNet .* expFilter;
    1 - sum(sum(sib > 0))/ sum(sum(tNet > 0))
end

sum(sum(sumNet > 0))

expFilter = zeros(18494, 18494);
inGenes = expMat(:,t);
expFilter(logical(inGenes), logical(inGenes)) = 1;

sib = sumNet .* expFilter;

sum(sum(sib > 0))/ sum(sum(sumNet > 0))

% How many of the brain TSNs, have a link which is not in at least
% one of the other tissues, How many of links have genes in not 2,
% in not 3, not 4, not 5 ... 

% Given a blood TAN, if one or two of genes not expressed in 
t = 1
atnNet = atn{t};

% for the links in atn, how many of them has both genes in one
% tissue, in two tissues, in three tissues... 

% get the exp matrix: 
expNet = zeros(18494, 18494);
for t = 1:5
    inGenes = logical(expMat(:,t));
    expNet(inGenes, inGenes) = expNet(inGenes, inGenes) + 1;
end

t = 1;
indAtn = atn{t};
indTsn = finalTable(t).wholeNet;
inGenes = logical(expMat(:,t));
expFilter(inGenes, inGenes) = 1;

sib = indAtn .* expNet;

% what percent of sib == 1, are TSN links ...
k = 1
kado = sib == k;
port = kado .* indTsn;
sum(sum(sib == k))
sum(sum(indTsn))
sum(sum(port))
sum(sum(port))/sum(sum(sib == k))

% what percent of sib == 2, are TSN links ... 
k = 2
kado = sib == k;
port = kado .* indTsn;
sum(sum(sib == k))
sum(sum(indTsn))
sum(sum(port))
sum(sum(port))/sum(sum(sib == k))

% what percent of sib == 3, are TSN links ...
k = 3
kado = sib == k;
port = kado .* indTsn;
sum(sum(sib == k))
sum(sum(indTsn))
sum(sum(port))
sum(sum(port))/sum(sum(sib == k))

% what percent of sib == 4 are TSN
k = 4
kado = sib == k;
port = kado .* indTsn;
sum(sum(sib == k))
sum(sum(indTsn))
sum(sum(port))
sum(sum(port))/sum(sum(sib == k))

% what percent of sib == 5 are TSN
k = 5
kado = sib == k;
port = kado .* indTsn;
sum(sum(sib == k))
sum(sum(indTsn))
sum(sum(port))
sum(sum(port))/sum(sum(sib == k))

% >>>>>> PORTIONS 1

portions = zeros(5,5);
for t = 1:5
    indTsn = finalTable(t).wholeNet;
    %    indTsn = atn{t};

    book = sum(expMat') .* expMat(:, t)';
    sib = indTsn .* expNet;

    for k = 1:5
        portions(t, k) = sum(sum(sib == k)) / (sum(book == k) * (sum(book==k)-1)/ 2);
    end
end

c1 = [228, 26, 28]/256
c2 = [254, 217, 166]/256
c3 = [166, 86, 40]/256
c4 = [55, 126, 184]/256
c5 = [247, 129, 191]/256
h = figure; 
hold all
plot(log10(portions(1, :)), 'color', c1)
plot(log10(portions(2, :)), 'color', c2)
plot(log10(portions(3, :)), 'color', c3)
plot(log10(portions(4, :)), 'color', c4)
plot(log10(portions(5, :)), 'color', c5)

set(gca, 'XTick', [1:5], 'XTickLabel', [1:5])
xlim([1 5.25])

legend(tissues)

figFolder = '~/resultsAndFigures/TSlinks/figures/'
file = sprintf('%sdropBasedOnExp', figFolder);
print(h, '-depsc', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])


% brain TSN links have %this of links with genes expressed in all
% tissues. 

% however, this is disproportional cause brain ATN has %this
% percent of links with genes expressed in all tissues. 

% >>>>>> PORTIONS 2
portions = zeros(5,5);
for t = 1:5
    indTsn = finalTable(t).wholeNet;
    %    indTsn = atn{t};

    % get all the links that have one or two of the InGenes
    thisInGenes = expMat(:, t);
    sumInGenes = sum(expMat');
    for k = 1:5
        
        % Getting the count of repeats for genes expressed in the
        % tissue == t
        thisSumInGenes = thisInGenes .* sumInGenes';
        
        % Getting the genes with count of repeats == k from the
        % thisSumInGenes
        kGenes = thisSumInGenes == k;
        
        % Getting the benes with count of repeats > k from the
        % thisSumInGenes
        kGreaterGenes = thisSumInGenes > k;
        
        gc1 = sum(kGenes);
        
        % How many links expected for the kGenes
        ls1Expected = gc1 * (gc1 -1) /2;
        
        % how many links expected between those and others
        ls2Expected = gc1 * sum(kGreaterGenes);
        
        % count of links between kGenes only 
        ls1Count = sum(sum(indTsn(kGenes, kGenes)));
        
        % count of links between kGenes and others
        ls2Count = sum(sum(indTsn(kGenes, kGreaterGenes))) + ...
            sum(sum(indTsn(kGreaterGenes, kGenes)));
        
        %portions(t, k) = ls1Count / ls1Expected;
        portions(t, k) = (ls1Count + ls2Count)/(ls1Expected + ls2Expected);
    end
end

c1 = [228, 26, 28]/256
c2 = [254, 217, 166]/256
c3 = [166, 86, 40]/256
c4 = [55, 126, 184]/256
c5 = [247, 129, 191]/256
h = figure; 
hold all
plot(log10(portions(1, :)), 'color', c1)
plot(log10(portions(2, :)), 'color', c2)
plot(log10(portions(3, :)), 'color', c3)
plot(log10(portions(4, :)), 'color', c4)
plot(log10(portions(5, :)), 'color', c5)

set(gca, 'XTick', [1:5], 'XTickLabel', [1:5])
xlim([1 5.25])

legend(tissues)

figFolder = '~/resultsAndFigures/TSlinks/figures/'
file = sprintf('%sdropBasedOnExp', figFolder);
print(h, '-depsc', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])

% 21. The node degree distribution for genes which have pure link
% vs those who don't, between the TS links
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 22. The functional enrichment for genes which have pure links,
% versus those who do not
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t = 1:5
    tsNet = finalTable(t).wholeNet;
    
    tsNds = sum(tsNet + tsNet');
    expTemp = zeros(18494, 18494);
    expTemp(cGenes, cGenes) = 1;

    pureNet = expTemp .* finalTable(t).noTSNet;
    pureNet = (finalTable(t).cgNet + finalTable(t).noTSNet) == 2;
    pureNds = sum(pureNet + pureNet');
end

inGenes = tsNds > 0;
subNDs = tsNds > 0;

nonPure = tsNds - pureNds;

a = full(tsNds(inGenes)');
b = full(pureNds(inGenes)');
corr(a, b)

% >>> hubs of genes with no pure are enriched in 
inGenes = ((tsNds > 0) +  ~(pureNds > 0)) == 2;
sum(inGenes) 
hist(tsNds(inGenes))
subNDs = tsNds(inGenes);
thr = quantile(subNDs, .90)
sum(subNDs > thr)
genesInds = find(inGenes);
finalListTS = genesInds(subNDs > thr);
hist(tsNds(finalListTS))

brainHubs.TS = finalListTS;
bloodHubs.TS = finalListTS;
liverHubs.TS = finalListTS;
lungHubs.TS = finalListTS;
muscleHubs.TS = finalListTS;
% they are enriched in...

% >>> hubs of genes with pure enriched with
inGenes = pureNds > 0;
sum(inGenes)
genesInds = find(inGenes);
subNDs = pureNds(inGenes);
thr = quantile(subNDs, .50)
hist(subNDs)
sum(subNDs > thr)
finalListPure = genesInds(subNDs > thr);

bloodHubs.pure = finalListPure;
brainHubs.pure = finalListPure;
liverHubs.pure = finalListPure;
lungHubs.pure = finalListPure;
muscleHubs.pure = finalListPure;


save('~/resultsAndFigures/geneEnrichments/pureHubsTSHubs/muscleHubs.mat', ...
     'muscleHubs')

save('~/resultsAndFigures/geneEnrichments/pureHubsTSHubs/lungHubs.mat', ...
     'lungHubs')

save('~/resultsAndFigures/geneEnrichments/pureHubsTSHubs/liverHubs.mat', ...
     'liverHubs')

save('~/resultsAndFigures/geneEnrichments/pureHubsTSHubs/bloodHubs.mat', ...
     'bloodHubs')

save('~/resultsAndFigures/geneEnrichments/pureHubsTSHubs/brainHubs.mat', ...
     'brainHubs')

geneSet = brainHubs.pure;
geneSet = bloodHubs.pure;
geneSet = liverHubs.TS;
geneSet = lungHubs.TS;
geneSet = muscleHubs.pure;
 
load('~/resultsAndFigures/functionalIdentityForGenes/sofUtil.mat')

k = 2
g = geneSet(k)
tempID = find(SOFutil.myGenes == g)
geneFMat = reshape(SOFutil.tsfs(tempID, :, :), 5, 3507);
sum(geneFMat')
geneFMat * geneFMat'
GOdata.geneSymbols(g)
k = k + 1

% 23. TANs and TSNs node degree distributions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

finalTable
atn

for t = 1:5
    fullNet = finalTable(t).wholeNet + finalTable(t).wholeNet';
    nds = sum(fullNet);
    [fi, xi] = ksdensity(nds(nds>0), 'Support', 'positive');
end
h = figure
plot(xi, fi)

% 24. Getting the plot for the STD of the partners of pure links
% and their means. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the pure links

pureSum = zeros(18494, 18494);
for t  = 1:5
    pureNet = (finalTable(t).cgNet + finalTable(t).noTSNet) == 2;
    pureSum = pureSum + pureNet;
end
% get the genes with link in the pureNet 
fullPureSum = pureSum + pureSum';
book = sum(pureSum + pureSum');

% get the quantile normalized expressions 
load('~/data/general/linkExprInfo/dataSetExpInf.mat') % load
                                                      % dataSetExpInf
totalExp = zeros(18494, 53);
for d = 1:53
    totalExp(:, d) = dataSetExpInf(d).meanExp;
end
qnTotalExp = quantilenorm(totalExp);

sib = std(qnTotalExp');

pureGenes = book > 0;
dummy = sib(logical(pureGenes));
h = figure
hist(dummy)
title(['Distribution of std for the 6658 genes which are partners ' ...
       'of pure links'])
figFolder = '~/resultsAndFigures/commonAndPure/figures/'
file = sprintf('%spureGenesDists', figFolder);
print(h, '-depsc', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])

find(((dummy<.3) + (dummy >0)) == 2)

find(fullPureSum(2757, :))

% 25. Getting the TS genes. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the quantile expression matrix for datasets. 
load('~/data/general/linkExprInfo/dataSetExpInf.mat') % load
                                                      % dataSetExpInf
totalExp = zeros(18494, 53);
for d = 1:53
    totalExp(:, d) = dataSetExpInf(d).meanExp;
end
qnTotalExp = quantilenorm(totalExp);


% get the genes expressed in more than one tissue. 
expMat = zeros(18494, 5);
for t = 1:5
    tissue = tissues{t};
    load(['~/data/general/tissueExpGenes/' tissue ...
          'ExpGenes0.8.mat'])
    expMat(:, t) = expGenesInd;
end

% blood 
thisT = [1:9]
otherT = [10:53];
t = 1

% lung
thisT = [10:24]
otherT = [1:9, 25:53];
t = 4

% muscle
thisT = [25:31]
otherT = [1:24, 32:53];
t = 5

% liver
thisT = [32:41]
otherT = [1:31, 42:53];
t = 3

% brain
thisT = [42:53]
otherT = [1:41];
t = 2

hs = zeros(18494, 1);
ps = zeros(18494, 1);
for i = 1:18494
    [h, p] = ttest2(qnTotalExp(i, thisT), qnTotalExp(i, otherT));
    ps(i) = p;
    hs(i) = h;
end

inGenes = logical(expMat(:, t));
hist(ps)
ps = ps * sum(inGenes);
ps = ps(inGenes);
h = figure
hist(ps)
m1 = median(qnTotalExp(inGenes, otherT)');
m2 = median(qnTotalExp(inGenes, thisT)');
fch = m2./m1;

inGeneList = find(inGenes);

[a, b] = sort(ps);
thr = max(find(a <= .1))
sortfch = fch(b);

passfch = sortfch >= 2;
sum(passfch(1:thr) > 0)

sortGeneList = inGeneList(b);
tissueSpecificGenes = sortGeneList(passfch(1:thr) > 0);

% blood specific genes

% 26. Overlap of the neighbours for partners of pure links versus
% semi pure and others. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% 27. DE genes between the tissues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
load('~/data/general/linkExprInfo/dataSetExpInf.mat') % load
                                                      % dataSetExpInf
totalExp = zeros(18494, 53);
for d = 1:53
    totalExp(:, d) = dataSetExpInf(d).meanExp;
end
qnTotalExp = quantilenorm(totalExp);

avgExp = zeros(18494, 5);
varExp = zeros(18494, 5);

avgExp(:, 1) = mean(qnTotalExp(:, 1:10)'); % blood
avgExp(:, 2) = mean(qnTotalExp(:, 11:25)'); % lung
avgExp(:, 3) = mean(qnTotalExp(:, 26:31)'); % muscle
avgExp(:, 4) = mean(qnTotalExp(:, 32:41)'); % liver
avgExp(:, 5) = mean(qnTotalExp(:, 42:end)'); % brain

varExp(:, 1) = var(qnTotalExp(:, 1:10)'); % blood
varExp(:, 2) = var(qnTotalExp(:, 11:25)'); % lung
varExp(:, 3) = var(qnTotalExp(:, 26:31)'); % muscle
varExp(:, 4) = var(qnTotalExp(:, 32:41)'); % liver
varExp(:, 5) = var(qnTotalExp(:, 42:end)'); % brain

expMat = zeros(18494, 5);
for t = 1:5
    tissue = tissues{t};
    load(['~/data/general/tissueExpGenes/' tissue ...
          'ExpGenes0.8.mat'])
    expMat(:, t) = expGenesInd;
end

groups = [1*ones(1,9) 2*ones(1, 15) 3*ones(1, 7) 4*ones(1,10), ...
          5*ones(1,12)]

ps = zeros(1, 18494);
for i = 1:18494
    ps(i) = anova1(qnTotalExp(i, :), groups, 'off');
end

inGenes = sum(expMat') > 0;

hist(ps(inGenes))

maxInds = zeros(18494, 2);
maxVs = zeros(18494, 2);
for i = 1:18494
    book = avgExp(i, :);
    [a, b] = sort(book, 'descend'); 
    maxInds(i, :) = b(1:2);
    maxVs(i, :) = a(1:2);
end

fcs = maxVs(:, 1) ./ maxVs(:, 2);

expGenesInds = find(inGenes);
selectedPs = ps(inGenes);

correctedPs = selectedPs .*length(selectedPs);
selectedFCs = fcs(inGenes);

passPs = correctedPs < .1;
passFCs = selectedFCs >= 1.5;

book = (passPs + passFCs') == 2;

finalGeneList = expGenesInds(book);
finalGeneListTissue = maxInds(book , 1);

TSGenes.GeneList = finalGeneList;
TSGenes.tissue = finalGeneListTissue;
TSGenes.rawPs = ps;
TSGenes.fc = fcs; 

save('~/resultsAndFigures/TSgenes.mat', 'TSGenes')

% 28. the functions enriched in tissues 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
% 1. get the functions which are exclusively enriched in each of
% the TS networks: 

load('~/resultsAndFigures/firstProject/functionsRepresented/tsRes02_HG.mat')
load('~/resultsAndFigures/firstProject/functionsRepresented/taRes_HG.mat')

tsExpSigTerms = zeros(5, 3507);
for t = 1:5
    tsst = tsRes.sigTerms;
    tast = taRes.sigTerms;
    tsst(t, :) = [];
    tast(t, :) = [];
    book = tsRes.sigTerms(t, :);
    sib = sum(tsst) + sum(tast);
    final = (book - sib) > 0;
    tsExpSigTerms(t, :) = final;
end

% 2. remove the TS genes, as the genes marked as expressed only in
% one tissue and test which functions are still enriched. 

tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'}
expMat = zeros(18494, 5);
for t = 1:5
    tissue = tissues{t};
    load(['~/data/general/tissueExpGenes/' tissue ...
          'ExpGenes0.8.mat'])
    expMat(:, t) = expGenesInd;
end

TSExpMarkedGenes = zeros(18494, 5);
for t = 1:5
    book = expMat;
    book(:, t) = [];
    sib = sum(book');
    final = (expMat(:, t) - sib') > 0;
    TSExpMarkedGenes(:, t) = final;
end

save('~/resultsAndFigures/expInOneTissueGenes.mat', 'TSExpMarkedGenes')

% load the TSNs 
FC = '3'
FDR = '0012'
load(sprintf(['~/resultsAndFigures/firstProject/TSlinks/finalTables/' ...
              'finalTable_CG13_FC%slog_FDR%s.mat'], FC, FDR))

% remove the ExpInduced links
for t = 1:5
    tsNet = finalTable(t).wholeNet;
    tsNet(:, logical(TSExpMarkedGenes(:, t))) = 0;
    tsNet(logical(TSExpMarkedGenes(:, t)), :) = 0;
    tsNetNonExpInd{t} = tsNet;
end

load('~/resultsAndFigures/firstTSgenes.mat')
for t = 1:5
    tsNet = finalTable(t).wholeNet;
    mytgenes = TSGenes.GeneList(TSGenes.tissue == t);
    tsNet(:, mytgenes) = 0;
    tsNet(mytgenes, :) = 0;
    tsNetNotTSgenes{t} = tsNet;
end

% 3. remove the pure link
for t = 1:5
    tsNet = finalTable(t).wholeNet;
    pureNet = (finalTable(t).cgNet + finalTable(t).noTSNet) == 2;
    tsNetNonPure{t} = tsNet - pureNet;
end

% 4. get the sigterms for non pure and non exp Induced
load('~/resultsAndFigures/firstProject/functionalIdentityForGenes/sofUtil.mat')
fMat = SOFutil.fMat;
clear SOFutil

fCount = size(fMat, 2);
gc = zeros(5, fCount);
bgExp = zeros(5, fCount);
obCount = zeros(5, fCount);
LRs = zeros(5, fCount);
hgps = zeros(5, fCount);
binps = zeros(5, fCount);
netDensities = [.0004, .004, .0011, .0005,.0004]; % densities for
                                                  % the TSNs
for t = 1:5
    % >> for non pure:
    %      myNet = tsNetNonPure{t};
    
    % >> for removal of genes marked as expressed only in one tissue
        myNet = tsNetNonExpInd{t};
    
    % >> for remval of the genes marked as TS 
    %     myNet = tsNetNotTSgenes{t};
     sum(sum(myNet))
    for f = 1:fCount
        myTerm = fMat(:, f);
        sum(myTerm)

        % how many genes are marked as expressed in this term?
        eGenesCount = sum((myTerm + expMat(:, t)) == 2)
        eGenes = ((myTerm + expMat(:, t)) == 2);
        gc(t, f) = eGenesCount;
        % for ppin >>>>>>>>>>>
        % eGenesCount = sum((myTerm + expArr') == 2);
        % eGenes = ((myTerm + expArr') == 2);
        % <<<<<<<<<<<<<<<<<<<<

        % 81 genes with the term, 50 genes expressed. 
        % if this FSim was fully represented in this network, it would have
        % all the 50*49/2 links         % 81 genes with the term, 50 genes expressed. 

        % now let's see how many links are there, and how is it bigger than
        % the background expectations. 
        % the background expectation for the brain network is 0.02 

        if(eGenesCount >= 20)

            bgExp(t, f) = (eGenesCount * (eGenesCount-1)/2) * netDensities(t);
            litNet = myNet(eGenes, eGenes);
            obCount(t, f) = sum(sum(litNet));
            LRs(t, f) = (sum(sum(litNet))/((eGenesCount * (eGenesCount-1)/2) ...
                                        * netDensities(t) + 1));
                        
            N = (eGenesCount * (eGenesCount-1)/2);
            k = N * netDensities(t);
            x = sum(sum(litNet));
            popSize = size(myNet, 1) * (size(myNet,1)-1)/2;
            hgps(t, f) = 1 - hygecdf(x, popSize, sum(sum(myNet)),N);
            %            binps(t, f) = 1 - binocdf(x, N, netDensities(t));
        else 
            bgExp(t, f) = -1;
            obCount(t, f) = -1;
            LRs(t, f) = -1;
            %            binps(t, f) = -1;
            hgps(t, f) = -1;
        end
    end
end

% selected Fs for each of the tissues: 
sigTerms = zeros(5, fCount);
for t = 1:5
    selected = hgps(t, :) > -1;
    sLRs = LRs(t, selected);
    shgps = hgps(t, selected);
    sterms = find(selected);

    % correcting pvalues
    shgps = shgps * sum(selected);
    logps = log10(shgps);
    logps(logps < -10) = -10;
    logps = logps .* -1;
    %    plot(logRs, logps, '.')

    temp1 = (logps >= 2);

    sigTerms(t, sterms(temp1)) = 1;
end

tsNonTSGenesRes.bgExp = bgExp;
tsNonTSGenesRes.obCount = obCount;
tsNonTSGenesRes.LRs = LRs;
tsNonTSGenesRes.hgps = hgps;
tsNonTSGenesRes.binps = binps;
tsNonTSGenesRes.sigTerms = sigTerms;
tsNonTSGenesRes.terms = tsRes.terms;
tsNonTSGenesRes.termIDs = tsRes.termIDs;
tsNonTSGenesRes.gc = gc;

save('~/resultsAndFigures/firstProject/functionsRepresented/tsNonTSGenesRes_HG.mat', ...
     'tsNonTSGenesRes')

load('~/resultsAndFigures/functionsRepresented/tsNonTSGenesRes.mat')

tsNonExpIndRes.bgExp = bgExp;
tsNonExpIndRes.obCount = obCount;
tsNonExpIndRes.LRs = LRs;
tsNonExpIndRes.hgps = hgps;
tsNonExpIndRes.binps = binps;
tsNonExpIndRes.sigTerms = sigTerms;
tsNonExpIndRes.terms = tsRes.terms;
tsNonExpIndRes.termIDs = tsRes.termIDs;
tsNonExpIndRes.gc = gc;

save('~/resultsAndFigures/firstProject/functionsRepresented/tsNonExpIndRes_HG.mat', ...
     'tsNonExpIndRes')
load('~/resultsAndFigures/functionsRepresented/tsNonExpIndRes.mat')

tsNonPureRes.bgExp = bgExp;
tsNonPureRes.obCount = obCount;
tsNonPureRes.LRs = LRs;
tsNonPureRes.hgps = hgps;
tsNonPureRes.binps = binps;
tsNonPureRes.sigTerms = sigTerms;
tsNonPureRes.terms = tsRes.terms;
tsNonPureRes.termIDs = tsRes.termIDs;
tsNonPureRes.gc = gc;

save('~/resultsAndFigures/firstProject/functionsRepresented/tsNonPureRes_HG.mat', ...
     'tsNonPureRes')

load('~/resultsAndFigures/firstProject/functionsRepresented/tsNonPureRes_HG.mat')
load('~/resultsAndFigures/firstProject/functionsRepresented/tsNonExpIndRes_HG.mat')
load('~/resultsAndFigures/firstProject/functionsRepresented/taRes_HG.mat')
load('~/resultsAndFigures/firstProject/functionsRepresented/tsRes02_HG.mat')

tsExpSigTerms = zeros(5, 3507);
for t = 1:5
    tsst = tsRes.sigTerms;
    tast = taRes.sigTerms;
    tsst(t, :) = [];
    tast(t, :) = [];
    book = tsRes.sigTerms(t, :);
    sib = sum(tsst) + sum(tast);
    final = (book - sib) > 0;
    tsExpSigTerms(t, :) = final;
end

% we are working with <<< tsExpSigTerms >>>> (for the exclusive
% terms), and tsNonPureRes (for the attributes) and the
% tsNonExpIndRes (for the attributes)

% load(['~/resultsAndFigures/functionsRepresented/' ...
%       'sampleTSGenesRemainSigTerms.mat'])
% tsRandomSigTerms = sampleRemainSigTerms;
% clear sampleRemainSigTerms;
clear
load(['~/resultsAndFigures/functionsRepresented/' ...
      'sampleEIlink_countOfPure_RemainSigTerms.mat'])
eiPureCountRandomSigTerms = sampleRemainSigTerms;
clear sampleRemainSigTerms;

load(['~/resultsAndFigures/functionsRepresented/' ...
      'sampleEIlinkRemainSigTerms.mat'])
eiRandomSigTerms = sampleRemainSigTerms;
clear sampleRemainSigTerms;

load(['~/resultsAndFigures/functionsRepresented/' ...
      'samplePureRemainSigTerms.mat'])
pureRandomSigTerms = sampleRemainSigTerms;
clear randomSigTerms

tsRandomRemainCount = zeros(1, 3507);

% load(['~/resultsAndFigures/functionsRepresented/' ...
%       'timesTermsMissed_sampleLinkSelection_TSGenesCount_.mat'])
% tsIndFloss = indFloss;
% clear indFloss;
load(['~/resultsAndFigures/functionsRepresented/' ...
      'timesTermsMissed_sampleLinkSelectionFromEILinks_pureLinkCount.mat'])
eiPureCountIndFloss = indFloss;
clear indFloss;
     
load(['~/resultsAndFigures/functionsRepresented/' ...
      'timesTermsMissed_sampleLinkSelection_EIlinkCount.mat'])
eiIndFloss = indFloss;
clear indFloss;

load(['~/resultsAndFigures/functionsRepresented/' ...
      'timesTermsMissed_sampleLinkSelection_pureCounts.mat'])
pureIndFloss = indFloss;

tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'};
totalContributingEICount = zeros(1, 5);
totalContributingPureCount = zeros(1, 5);
for t = 2:4
    file = sprintf(['~/resultsAndFigures/functionsRepresented/' ...
    '%sTSNExclusiveTerms_table_withLog10Pvalues_withCounts_withEIPureCountRemoval.txt'], tissues{t})
    fid = fopen(file, 'w')
    fprintf(fid, ['functionID\tfunctionName\tgeneCount\' ...
                  'tlinkCount(lc)\tlc_minusPurLinks\' ...
                  'tlc_minusExpressionInducedLinks\t1-ratioOfPureLinks' ...
                  '\t1-RatioOfExpressionInducedLinks\' ...
                  'tSigAfterPureRemoval\' ...
                  'tSigAfterExpressionInducedRemoval\t'...
                  'log10Pvalue_plus1e-20\t'...
                  'log10PvalueAfterPureRemoval_plus1e-20\t'...
                  'log10PvalueAfterEIremovalplus1e-20\t'...
                  'significantCountIn100_randomPureLinkRemoval\t'...
                  'significantCountIn100_randomEILinkRemoval_pureCount\t'...
                  'significantCountIn100_randomEILinkRemoval\n'])
    termCount = sum(tsExpSigTerms(t, :));
    terms = find(tsExpSigTerms(t, :));
    for i = 1:termCount
        % for the first few I can use tsRes or tsNonPureRes or
        % tsNonExpIndRes. they have the same values
        thisTerm = terms(i);
        nonPureObRatio = tsNonPureRes.obCount(t, thisTerm) / ...
            tsRes.obCount(t, thisTerm);
        nonTsObRatio = tsNonExpIndRes.obCount(t, thisTerm) / ...
            tsRes.obCount(t, thisTerm);
        
        totalContributingEICount(t) = totalContributingEICount(t) + ...
            (tsRes.obCount(t, thisTerm) - tsNonExpIndRes.obCount(t, thisTerm));
        totalContributingPureCount(t) = totalContributingPureCount(t) + ...
                    (tsRes.obCount(t, thisTerm) - tsNonPureRes.obCount(t, thisTerm));
                
        fprintf(fid, ['%d\t%s\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%d\t'...
                      '%d\t%.2f\t%.2f\t%.2f\t%d\t%d\t%d\n'], ...
                tsRes.termIDs(thisTerm),...  
                tsRes.terms{thisTerm},...
                tsRes.gc(t, thisTerm),...
                tsRes.obCount(t, thisTerm),...
                tsNonPureRes.obCount(t, thisTerm),...
                tsNonExpIndRes.obCount(t, thisTerm),...
                nonPureObRatio,...
                nonTsObRatio,...
                tsNonPureRes.sigTerms(t, thisTerm),...
                tsNonExpIndRes.sigTerms(t, thisTerm),...
                log10(tsRes.binps(t, thisTerm) + 1e-20),...
                log10(tsNonPureRes.binps(t, thisTerm) + 1e-20),...
                log10(tsNonExpIndRes.binps(t, thisTerm) + 1e-20),...
                pureIndFloss(t, thisTerm),...
                eiPureCountIndFloss(t, thisTerm),...
                eiIndFloss(t, thisTerm));
    end
    fclose(fid)
end

% 29. print the TSN and TAN functions in tables 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
load('~/resultsAndFigures/functionsRepresented/taRes.mat')
load('~/resultsAndFigures/functionsRepresented/tsRes.mat')

sib = taRes.sigTerms;
sib = sib + tsRes.sigTerms;

kado = sum(sib) > 0;
inds = find(kado);

file = sprintf(['~/resultsAndFigures/functionsRepresented/' ...
                'enrichedTerms_TAN_TSN_table.csv'])
fid = fopen(file, 'w')
fprintf(fid, ['functionID\tfunctionTerm\tblood_TAN\tblood_TSN\t' ...
              'brain_TAN\tbrainTSN\tliver_TAN\tliver_TSN\tlung_TAN\t' ...
              'lung_TSN\tmuscle_TAN\tmuscle_TSN\n'])

for i = 1:length(inds)

    actInd = inds(i);
    thisTerm = tsRes.terms(actInd);
    thisTermID = tsRes.termIDs(actInd);
    
    fprintf(fid, ['%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n'], ...
            taRes.termIDs(actInd),...  
            taRes.terms{actInd},...
            taRes.sigTerms(1, actInd),...
            tsRes.sigTerms(1, actInd),...
            taRes.sigTerms(2, actInd),...
            tsRes.sigTerms(2, actInd),...
            taRes.sigTerms(3, actInd),...
            tsRes.sigTerms(3, actInd),...
            taRes.sigTerms(4, actInd),...
            tsRes.sigTerms(4, actInd),...
            taRes.sigTerms(5, actInd),...
            tsRes.sigTerms(5, actInd))
end

% 30. random selection of genes instead of TSN and pure. I need to
% do 30 for each.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'};

expMat = zeros(18494, 5);
for t = 1:5
    tissue = tissues{t};
    load(['~/data/general/tissueExpGenes/' tissue ...
          'ExpGenes0.8.mat'])
    expMat(:, t) = expGenesInd;
end

% get fMat, fCount
% Tissues And final table
load('~/resultsAndFigures/functionalIdentityForGenes/sofUtil.mat')
fMat = SOFutil.fMat;
clear SOFutil;
fCount = size(fMat, 2);

FC = '3'
FDR = '0012'
load(sprintf(['~/resultsAndFigures/TSlinks/finalTables/' ...
              'finalTable_CG13_FC%slog_FDR%s.mat'], FC, FDR))

% >>>> load the five ATN and TSN networks

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

tsnLC = sum(tsnds')./2;

% getting pure nets
for t = 1:5
    pureNet{t} = (finalTable(t).cgNet + finalTable(t).noTSNet) == ...
        2;
    pureLC(t) = sum(sum(pureNet{t}));
end

% >>>>> pick one
% LCs = tsnLC;
LCs = pureLC;

fCount = size(fMat, 2);
gc = zeros(5, fCount);
bgExp = zeros(5, fCount);
obCount = zeros(5, fCount);
LRs = zeros(5, fCount);
hgps = zeros(5, fCount);
binps = zeros(5, fCount);
netDensities = [.007, .020, .008, .012, .006];
netDensities = [.0004, .004, .0011, .0005,.0004];

% for ppin
% netDensities = sum(sum(newPPI)) / (15039 * 15039 /2);
% expArr = (sum(myNet) + sum(myNet')) >0;

for sit = 1:30
    gc = zeros(5, fCount);
    bgExp = zeros(5, fCount);
    obCount = zeros(5, fCount);
    LRs = zeros(5, fCount);
    hgps = zeros(5, fCount);
    binps = zeros(5, fCount);

    sit
    for t = 1:5
        
        % >>>>> for the TSN
        % myNet = finalTable(t).wholeNet;
        % baseNet = atn{t};
        
        % >>>>> for the pure
        myNet = pureNet{t};
        baseNet = finalTable(t).wholeNet;
        
        % sampling count of links for TSN
        
        [l1, l2, gar] = find(myNet);
        sampleLs = datasample([1:length(l1)], LCs(t));
        myNet = sparse(l1(sampleLs), l2(sampleLs), ones(length(sampleLs), 1), ...
                       18494, 18494);
        
        % >>>>> for the pure 
        

        for f = 1:fCount
            myTerm = fMat(:, f);
            sum(myTerm)

            % how many genes are marked as expressed in this term?
            eGenesCount = sum((myTerm + expMat(:, t)) == 2)
            eGenes = ((myTerm + expMat(:, t)) == 2);
            gc(t, f) = eGenesCount;
            % for ppin >>>>>>>>>>>
            % eGenesCount = sum((myTerm + expArr') == 2);
            % eGenes = ((myTerm + expArr') == 2);
            % <<<<<<<<<<<<<<<<<<<<

            % 81 genes with the term, 50 genes expressed. 
            % if this FSim was fully represented in this network, it would have
            % all the 50*49/2 links 

            % now let's see how many links are there, and how is it bigger than
            % the background expectations. 
            % the background expectation for the brain network is 0.02 

            if(eGenesCount >= 20)

                bgExp(t, f) = (eGenesCount * (eGenesCount-1)/2) * netDensities(t);
                litNet = myNet(eGenes, eGenes);
                obCount(t, f) = sum(sum(litNet));
                LRs(t, f) = (sum(sum(litNet))/((eGenesCount * (eGenesCount-1)/2) ...
                                               * netDensities(t) + 1));
                
                %            tGenes = sum(expMat(:, tissueID));
                N = (eGenesCount * (eGenesCount-1)/2);
                %            M = tGenes * (tGenes - 1) / 2;
                %            k = M * netDensity;
                x = sum(sum(litNet));
                %        hgps(f) = 1 - hygecdf(x, M, k,N);
                binps(t, f) = 1- binocdf(x, N, netDensities(t));
            else 
                bgExp(t, f) = -1;
                obCount(t, f) = -1;
                LRs(t, f) = -1;
                binps(t, f) = -1;
            end
            % so this function is not adding to that expectation, but which one
            % is? 
        end
    end

    % selected Fs for each of the tissues: 
    sigTerms = zeros(5, fCount);
    for t = 1:5
        selected = binps(t, :) > -1;
        sLRs = LRs(t, selected);
        sbinps = binps(t, selected);
        sterms = find(selected);

        % correcting pvalues
        sbinps = sbinps * sum(selected);
        logps = log10(sbinps);
        logps(logps < -10) = -10;
        logps = logps .* -1;
        %    plot(logRs, logps, '.')

        temp1 = (logps >= 2);

        sigTerms(t, sterms(temp1)) = 1;
    end
    randomSigTerms{sit} = sigTerms;
end

save('~/resultsAndFigures/functionsRepresented/randomTSNsigTerms.mat', ...
     'randomSigTerms')

% 31. This is in response to reviewers comment "the change in the
% number of enriched terms due to the loss of TSNs does not appear
% to be that high. If one sampled the initial set of edges but to a
% number that matches the reduction cased by removing pure links,
% would the count of enriched functions be similarly reduced?" 
% Here I am basically getting the enriched functions upon the
% removal of similar count of links as the pure links. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'};

expMat = zeros(18494, 5);
for t = 1:5
    tissue = tissues{t};
    load(['~/data/general/tissueExpGenes/' tissue ...
          'ExpGenes0.8.mat'])
    expMat(:, t) = expGenesInd;
end

% get fMat, fCount
% >>>>> pick for GO or PHENO
% GO Tissues And final table
load('~/resultsAndFigures/firstProject/functionalIdentityForGenes/sofUtil.mat')
fMat = SOFutil.fMat;
clear SOFutil;
fCount = size(fMat, 2);

% PHENO 
load('~/resultsAndFigures/firstProject/diseaseRepresented/dMat_pheno.mat')
fMat = phenoDs.mat > 0;
% <<<<<<<<<<

FC = '3'
FDR = '0012'
load(sprintf(['~/resultsAndFigures/firstProject/TSlinks/finalTables/' ...
              'finalTable_CG13_FC%slog_FDR%s.mat'], FC, FDR))

% >>>> load the five ATN and TSN networks

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

tsnLC = sum(tsnds')./2;

% getting pure nets
for t = 1:5
    pureNet{t} = (finalTable(t).cgNet + finalTable(t).noTSNet) == ...
        2;
    pureLC(t) = sum(sum(pureNet{t}));
end

% getting the expression induced link count (based on the marked expressed
% gene)
load('~/resultsAndFigures/firstProject/expInOneTissueGenes.mat')
for t = 1:5
    tsnet = finalTable(t).wholeNet;
    l1 = sum(sum(tsnet));
    tsnet(logical(TSExpMarkedGenes(:, t)), :) = 0;
    tsnet(:, logical(TSExpMarkedGenes(:, t))) = 0;
    l2 = sum(sum(tsnet));
    eiLC(t) = l1 - l2;
    eiNet{t} = finalTable(t).wholeNet - tsnet;
end

baseNetLCs = tsnLC;
fCount = size(fMat, 2);
netDensities = [.0004, .004, .0011, .0005,.0004];
for sit = 1:300
    gc = zeros(5, fCount);
    bgExp = zeros(5, fCount);
    obCount = zeros(5, fCount);
    LRs = zeros(5, fCount);
    hgps = zeros(5, fCount);
    binps = zeros(5, fCount);

    sit
    for t = 2:4
        
        % pick one
        %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        % >> for the TSN
        %myNet = finalTable(t).wholeNet;
        % baseNet = atn{t};
        
        % >> for the pure
        baseNet = finalTable(t).wholeNet;
        % <<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        % pick one
        %>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
        % >>>>>> sampling count of links for TSN (either pure or eiLC)
        % [l1, l2, gar] = find(baseNet);
        % %LC = pureLC(t);
        % LC = eiLC(t);
        % sampleLs = datasample([1:length(l1)], LC, 'Replace', ...
        %                       false);
        % l1(sampleLs) = [];
        % l2(sampleLs) = [];
        % myNet = sparse(l1, l2, ones(length(l1), 1), ...
        %                18494, 18494);
        
        % >>>>>> sampling EI links to the count of pure links (this
        %should be done only for brain, liver and lung)
        tempNet = baseNet + eiNet{t};
        [l1, l2, l3] = find(tempNet);
        
        % get the ei links inds
        sinds = find(l3 == 2);
        LC = pureLC(t);
        sampleLs = datasample(sinds, LC, 'Replace', false);
        l1(sampleLs) = [];
        l2(sampleLs) = [];
        myNet = sparse(l1, l2, ones(length(l1), 1), ...
                       18494, 18494);
        % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        for f = 1:fCount
            myTerm = fMat(:, f);
            sum(myTerm);

            % how many genes are marked as expressed in this term?
            eGenesCount = sum((myTerm + expMat(:, t)) == 2);
            eGenes = ((myTerm + expMat(:, t)) == 2);
            gc(t, f) = eGenesCount;
            % for ppin >>>>>>>>>>>
            % eGenesCount = sum((myTerm + expArr') == 2);
            % eGenes = ((myTerm + expArr') == 2);
            % <<<<<<<<<<<<<<<<<<<<

            % 81 genes with the term, 50 genes expressed. 
            % if this FSim was fully represented in this network, it would have
            % all the 50*49/2 links 

            % now let's see how many links are there, and how is it bigger than
            % the background expectations. 
            % the background expectation for the brain network is 0.02 

            if(eGenesCount >= 5)

                bgExp(t, f) = (eGenesCount * (eGenesCount-1)/2) * netDensities(t);
                litNet = myNet(eGenes, eGenes);
                obCount(t, f) = sum(sum(litNet));
                LRs(t, f) = (sum(sum(litNet))/((eGenesCount * (eGenesCount-1)/2) ...
                                               * netDensities(t) + 1));

                
                N = (eGenesCount * (eGenesCount-1)/2);
                k = N * netDensities(t);
                x = sum(sum(litNet));
                popSize = size(myNet, 1) * (size(myNet,1)-1)/2;
                hgps(t, f) = 1 - hygecdf(x, popSize, sum(sum(myNet)),N);
                %            binps(t, f) = 1 - binocdf(x, N, netDensities(t));
            else 
                bgExp(t, f) = -1;
                obCount(t, f) = -1;
                LRs(t, f) = -1;
                %            binps(t, f) = -1;
                hgps(t, f) = -1;
            end
        end
    end

    % selected Fs for each of the tissues: 
    sigTerms = zeros(5, fCount);
    for t = 1:5
        selected = hgps(t, :) > -1;
        sLRs = LRs(t, selected);
        shgps = hgps(t, selected);
        sterms = find(selected);

        % correcting pvalues
        shgps = shgps * sum(selected);
        logps = log10(shgps);
        logps(logps < -10) = -10;
        logps = logps .* -1;
        %    plot(logRs, logps, '.')

        temp1 = (logps >= 2);

        sigTerms(t, sterms(temp1)) = 1;
    end
    sampleRemainSigTerms{sit} = sigTerms;
end

save('~/resultsAndFigures/functionsRepresented/samplePureRemainSigTerms_Pheno.mat', ...
     'sampleRemainSigTerms')

save('~/resultsAndFigures/functionsRepresented/sampleEIlinkRemainSigTerms_pheno.mat', ...
     'sampleRemainSigTerms')

save('~/resultsAndFigures/firstProject/functionsRepresented/sampleEIlinkRemainSigTerms300_HG.mat', ...
     'sampleRemainSigTerms')

save('~/resultsAndFigures/functionsRepresented/sampleEIlink_countOfPure_RemainSigTerms_pheno.mat', ...
     'sampleRemainSigTerms')

save('~/resultsAndFigures/functionsRepresented/samplePureRemainSigTerms_Pheno.mat', ...
     'sampleRemainSigTerms')

save('~/resultsAndFigures/firstProject/functionsRepresented/samplePureRemainSigTerms300_HG.mat', ...
     'sampleRemainSigTerms')

load('~/resultsAndFigures/functionsRepresented/samplePureRemainSigTerms300.mat')

save('~/resultsAndFigures/firstProject/functionsRepresented/sampleEIlink_countOfPure_RemainSigTerms300_HG.mat', ...
     'sampleRemainSigTerms')

load(['~/resultsAndFigures/functionsRepresented/' ...
      'sampleEIlink_countOfPure_RemainSigTerms.mat'])

save('~/resultsAndFigures/functionsRepresented/sampleEIlinkRemainSigTerms.mat', ...
     'sampleRemainSigTerms')

load(['~/resultsAndFigures/functionsRepresented/' ...
      'sampleEIlinkRemainSigTerms.mat'])

save('~/resultsAndFigures/functionsRepresented/samplePureRemainSigTerms.mat', ...
     'sampleRemainSigTerms')

load(['~/resultsAndFigures/functionsRepresented/' ...
      'samplePureRemainSigTerms.mat'])

% I need to load the functions for original TSN and the minus-pure links
load('~/resultsAndFigures/functionsRepresented/tsNonPureRes.mat')
load('~/resultsAndFigures/functionsRepresented/sampleEIlink_countOfPure_RemainSigTerms.mat')
load('~/resultsAndFigures/firstProject/functionsRepresented/taRes_HG.mat')
load('~/resultsAndFigures/firstProject/functionsRepresented/tsRes02_HG.mat')

% to get the ts exclusive terms
tsExpSigTerms = zeros(5, 3507);
for t = 1:5
    tsst = tsRes.sigTerms;
    tast = taRes.sigTerms;
    tsst(t, :) = [];
    tast(t, :) = [];
    book = tsRes.sigTerms(t, :);
    sib = sum(tsst) + sum(tast);
    final = (book - sib) > 0;
    tsExpSigTerms(t, :) = final;
end

floss = zeros(300, 5);
% count of functions which are not significant anymore after
% removal of the random sets
for i = 1:300
    book = sampleRemainSigTerms{i};
    kado = (tsExpSigTerms - book) >0;
    floss(i, :) = sum((kado>0)');
end

save('~/resultsAndFigures/firstProject/functionsRepresented/affectedTerms_sampleLinkSelectionFromEILinks_pureCounts300_HG.mat', ...
     'floss')

load(['~/resultsAndFigures/functionsRepresented/' ...
      'affectedTerms_sampleLinkSelectionFromEILinks_pureCounts.mat'])
pcEIFloss = floss;
clear floss

save('~/resultsAndFigures/functionsRepresented/affectedTerms_sampleLinkSelection_pureCounts.mat', ...
     'floss') 

load(['~/resultsAndFigures/functionsRepresented/' ...
      'affectedTerms_sampleLinkSelection_pureCounts.mat'])
pcFloss = floss(101:200, :);

h = figure
hist(pcFloss(:, 2))
figFolder = '~/resultsAndFigures/functionsRepresented/figures/'
file = sprintf('%sbrainFloss_pureCount', figFolder);
print(h, '-depsc', [file '.eps'])
print(h, '-dpdf', [file '.pdf'])



indFloss = zeros(5, 3507);
for i = 1:100
    book = sampleRemainSigTerms{i};
    kado = (tsExpSigTerms .* book);
    indFloss = indFloss + kado;
end

save('~/resultsAndFigures/functionsRepresented/timesTermsMissed_sampleLinkSelectionFromEILinks_pureLinkCount.mat', ...
     'indFloss')

save('~/resultsAndFigures/functionsRepresented/timesTermsMissed_sampleLinkSelection_EIlinkCount.mat', ...
     'indFloss')

save('~/resultsAndFigures/functionsRepresented/timesTermsMissed_sampleLinkSelection_pureCounts.mat', ...
     'indFloss')

load('~/resultsAndFigures/functionsRepresented/timesTermsMissed_sampleLinkSelection_pureCounts.mat')

% that plot
book = tsNonPureRes.sigTerms;
kado = (tsExpSigTerms - book) >0;
pureLoss = sum((kado>0)');

for t = 2:4
    tissue = tissues{t};
    myfloss = floss(:, t);
    [f, x] = ksdensity(floss(:, t));
    h = figure
    plot(x, f)
    hold on
    line([pureLoss(t), pureLoss(t)], [0 ,.05])
    title(tissue)
    xlabel('#drop-out functions')
    figFolder = '~/resultsAndFigures/TSlinks/figures/'
    file = sprintf('%s%sRandomFunctionLoss', figFolder, tissue);
    print(h, '-depsc', [file '.eps'])
    print(h, '-dpdf', [file '.pdf'])
end

for t = 2:4
    tissue = tissues{t};
    myfloss = floss(:, t);
    h = figure
    hist(myfloss, 10)
    hold on
    line([pureLoss(t), pureLoss(t)], [0 ,10])
    xlim([min(min(myfloss), pureLoss(t))-2 max(myfloss)+2])
    title(tissue)
    xlabel('#drop-out functions')
    figFolder = '~/resultsAndFigures/TSlinks/figures/'
    file = sprintf('%s%sRandomFunctionLoss_hist', figFolder, tissue);
    print(h, '-depsc', [file '.eps'])
    print(h, '-dpdf', [file '.pdf'])
end

% now which one of them I lose when I pretend I don't have sample
% pure links. 

% 32. Count of functions lost by removing random sets of genes (to
% the count of TS genes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('~/resultsAndFigures/TSgenes.mat')
load('~/resultsAndFigures/expInOneTissueGenes.mat')

% >>>>>>>>> pick one!
% 1. TSgenes:getting the count of EI links from the TSNs
for t = 1:5
    tsgenes = TSGenes.GeneList(TSGenes.tissue == t);
    sib = finalTable(t).wholeNet;
    sib(tsgenes, :) = 0;
    sib(:, tsgenes) = 0;
    eic(t) = sum(sum(finalTable(t).wholeNet)) - sum(sum(sib));
    pc(t) = sum(sum(pureNet{t}));
end

% 2. genes marked as expressed only in one tissue
for t = 1:5
    tsgenes = logical(TSExpMarkedGenes(:, t));
    sib = finalTable(t).wholeNet;
    sib(tsgenes, :) = 0;
    sib(:, tsgenes) = 0;
    eic(t) = sum(sum(finalTable(t).wholeNet)) - sum(sum(sib));
    %    pc(t) = sum(sum(pureNet{t}));
end
TSgeneCounts = sum(TSExpMarkedGenes)

for sit = 1:100
    gc = zeros(5, fCount);
    bgExp = zeros(5, fCount);
    obCount = zeros(5, fCount);
    LRs = zeros(5, fCount);
    hgps = zeros(5, fCount);
    binps = zeros(5, fCount);

    sit
    for t = 1:5
        
        % >>>>> for the TSN
        % myNet = finalTable(t).wholeNet;
        % baseNet = atn{t};
        
        % >>>>> for the pure
        baseNet = finalTable(t).wholeNet;
        
        inGenes = sum(baseNet + baseNet') > 0;
        inGenesInds = find(inGenes);
        
        % sampling count of TS genes
        
        % count of ts genes to be removed (pick one)
        % 1. 
        % tsGeneCount = sum(TSGenes.tissue == t);
        % 2.
        tsGeneCount = TSgeneCounts(t);
        sampleGenes = datasample(inGenesInds, tsGeneCount, 'Replace', ...
                                 false);
        
        % removing ts genes from the network
        myNet = baseNet;
        
        myNet(sampleGenes, :) = 0;
        myNet(:, sampleGenes) = 0;
       
        for f = 1:fCount
            myTerm = fMat(:, f);
            sum(myTerm);

            % how many genes are marked as expressed in this term?
            eGenesCount = sum((myTerm + expMat(:, t)) == 2);
            eGenes = ((myTerm + expMat(:, t)) == 2);
            gc(t, f) = eGenesCount;
            % for ppin >>>>>>>>>>>
            % eGenesCount = sum((myTerm + expArr') == 2);
            % eGenes = ((myTerm + expArr') == 2);
            % <<<<<<<<<<<<<<<<<<<<

            % 81 genes with the term, 50 genes expressed. 
            % if this FSim was fully represented in this network, it would have
            % all the 50*49/2 links 

            % now let's see how many links are there, and how is it bigger than
            % the background expectations. 
            % the background expectation for the brain network is 0.02 

            if(eGenesCount >= 20)

                bgExp(t, f) = (eGenesCount * (eGenesCount-1)/2) * netDensities(t);
                litNet = myNet(eGenes, eGenes);
                obCount(t, f) = sum(sum(litNet));
                LRs(t, f) = (sum(sum(litNet))/((eGenesCount * (eGenesCount-1)/2) ...
                                               * netDensities(t) + 1));
                
                %            tGenes = sum(expMat(:, tissueID));
                N = (eGenesCount * (eGenesCount-1)/2);
                %            M = tGenes * (tGenes - 1) / 2;
                %            k = M * netDensity;
                x = sum(sum(litNet));
                %        hgps(f) = 1 - hygecdf(x, M, k,N);
                binps(t, f) = 1- binocdf(x, N, netDensities(t));
            else 
                bgExp(t, f) = -1;
                obCount(t, f) = -1;
                LRs(t, f) = -1;
                binps(t, f) = -1;
            end
            % so this function is not adding to that expectation, but which one
            % is? 
        end
    end

    % selected Fs for each of the tissues: 
    sigTerms = zeros(5, fCount);
    for t = 1:5
        selected = binps(t, :) > -1;
        sLRs = LRs(t, selected);
        sbinps = binps(t, selected);
        sterms = find(selected);

        % correcting pvalues
        sbinps = sbinps * sum(selected);
        logps = log10(sbinps);
        logps(logps < -10) = -10;
        logps = logps .* -1;
        %    plot(logRs, logps, '.')

        temp1 = (logps >= 2);

        sigTerms(t, sterms(temp1)) = 1;
    end
    sampleRemainSigTerms{sit} = sigTerms;
end

% this is the file for where genes are expressed based on if they
% are marked as expressed or not
save('~/resultsAndFigures/functionsRepresented/sampleMarkedTSGenesRemainSigTerms.mat', ...
     'sampleRemainSigTerms')

% this is the file for where genes are tissue-specific (based on
% FDR and pvalue)
save('~/resultsAndFigures/functionsRepresented/sampleTSGenesRemainSigTerms.mat', ...
     'sampleRemainSigTerms')

load(['~/resultsAndFigures/functionsRepresented/' ...
      'sampleTSGenesRemainSigTerms.mat'])

% I need to load the functions for original TSN and the minus-pure links
load('~/resultsAndFigures/functionsRepresented/tsNonPureRes.mat')
load('~/resultsAndFigures/functionsRepresented/taRes.mat')
load('~/resultsAndFigures/functionsRepresented/tsRes.mat')

% to get the ts exclusive terms
tsExpSigTerms = zeros(5, 3507);
for t = 1:5
    tsst = tsRes.sigTerms;
    tast = taRes.sigTerms;
    tsst(t, :) = [];
    tast(t, :) = [];
    book = tsRes.sigTerms(t, :);
    sib = sum(tsst) + sum(tast);
    final = (book - sib) > 0;
    tsExpSigTerms(t, :) = final;
end

floss = zeros(100, 5);
for i = 1:100
    book = sampleRemainSigTerms{i};
    kado = (tsExpSigTerms - book) >0;
    floss(i, :) = sum((kado>0)');
end

save('~/resultsAndFigures/functionsRepresented/affectedTerms_sampleGeneSelection_TSGeneCounts.mat', ...
     'floss')

load(['~/resultsAndFigures/functionsRepresented/' ...
      'affectedTerms_sampleGeneSelection_TSGeneCounts.mat'])

indFloss = zeros(5, 3507);
for i = 1:100
    book = sampleRemainSigTerms{i};
    kado = (tsExpSigTerms .* book);
    indFloss = indFloss + kado;
end
save('~/resultsAndFigures/functionsRepresented/timesTermsMissed_sampleLinkSelection_TSGenesCount.mat', ...
     'indFloss')

book = tsNonPureRes.sigTerms;
kado = (tsExpSigTerms - book) >0;
pureLoss = sum((kado>0)');

for t = 2:4
    tissue = tissues{t};
    myfloss = floss(:, t);
    [f, x] = ksdensity(floss(:, t));
    h = figure
    plot(x, f)
    hold on
    title(tissue)
    xlabel('#drop-out functions')
    figFolder = '~/resultsAndFigures/TSlinks/figures/'
    file = sprintf('%s%sRandomFunctionLoss_TSGeneCount', figFolder, tissue);
    print(h, '-depsc', [file '.eps'])
    print(h, '-dpdf', [file '.pdf'])
end


for t = 2:4
    tissue = tissues{t};
    myfloss = floss(:, t);
    h = figure
    hist(myfloss, 10)
    hold on
    xlim([(min(myfloss))-2 max(myfloss)+10])
    title(tissue)
    xlabel('#drop-out functions')
    figFolder = '~/resultsAndFigures/TSlinks/figures/'
    file = sprintf('%s%sRandomFunctionLoss_TSGeneCount_hist', figFolder, tissue);
    print(h, '-depsc', [file '.eps'])
    print(h, '-dpdf', [file '.pdf'])
end

% 33. This is another response to the reviewer comments, on the
% enrichment of the disease for genes with pure links
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% I will have different types of genes with pure links:
% 1. common genes which are hubs in TAN : I can let this go, as it
% has no pure link focus
% 2. common genes which are hubs in TSN : I can let this go too, as
% it has no pure link focus
% 3. common genes with *any* pure links (in any of the five tissues)
% 4. common genes which are hubs in pure (in one tissue)
% 5. common genes with multi-tissue pure links
% 6. Does any of the above group of genes, seem to be enriched in
% any kind of disease -tissue specific or not- ?
% what are the conclusions from this? 

% Are the pure genes more likely to have disease associated with
% them than the average commonly expressed genes? 

% how about functions

% how about correlation

% TODO: 
% 1. get the file 
% 2. get the annotations
% 3. get the enrichments

% read the file with 8 columns 
clear
fileName = ['/space/gemmaData/PhenocartaExport/' ...
            'PhenocartaExport_04-23-2018/AnnotationsWithOMIM/' ...
            'AllPhenocartaAnnotations.tsv']

linum = 0;
fid = fopen(fileName)
data = textscan(fid, repmat('%s', 1, 11), 'headerlines', linum, ...
                  'Delimiter', '\t')

% we have 38409 annotations and column 7 is the Disease Ontology
% link/ID, so we go for that.Column 3 is the gene symbols.
% column 5 is the disease name

% get the gene, make the matrix. 
dataFolder = '~/data/general/'
load([dataFolder 'GPL570GemmaMapNEW.mat'])

[a, b] = ismember(data{2}, gpl570.ncbiID);
%[a, b] = ismember(gpl570.uniqueSymbols, data{3});
% I have 35,862 records

% extracting unique DOIDs:
IDs = ones(1, 40615); 
Dnames = cell(1, 40615);
c = 1
idLengths = zeros(1, length(data{7}));
nameLengths = zeros(1, length(data{7}));
sInd = zeros(1, length(data{7}));
cc = 0;
for i = 1:length(data{7})
    i
    mystr = data{7}{i};
    sind = regexp(mystr, '/DOID_(\d{2,12})', 'tokens');
    tempNames = strsplit(data{5}{i}, ';');
    idLengths(i) = length(tempNames);
    nameLengths(i) = length(sind);
    if length(tempNames) == length(sind)
        cc = cc + 1;
        sInd(i) = cc;
        for j = 1:length(sind)
            IDs(c) = str2num(sind{j}{1});
            Dnames{c} = tempNames{j};
            c = c + 1;
        end
    end
end

length(unique(Dnames)) % 2303
length(unique(IDs)) % 2289

DOidNameMap = containers.Map(IDs, Dnames)

dMat = zeros(18494, length(unique(IDs)));
uniqueDOIDs = unique(IDs);
gs = 0;
for i = 1:length(data{7})
    i
    % find the gene
    [ag, bg] = ismember(data{2}{i}, gpl570.uniqueIDs);

    % find the IDs
    mystr = data{7}{i};
    sind = regexp(mystr, '/DOID_(\d{2,12})', 'tokens');
    tempIDs = zeros(1, length(sind));
    for j = 1:length(sind)
        tempIDs(j) = str2num(sind{j}{1});
    end
    [ad, bd] = ismember(tempIDs , uniqueDOIDs);
    
    % fill the mat
    
    if ag
        gs = gs + 1;
        dMat(bg, bd(ad)) = dMat(bg, bd(ad)) + 1;
    end
end

phenoDs.mat = dMat;
phenoDs.DOIDs = uniqueDOIDs;
phenoDs.map = DOidNameMap;
save('~/resultsAndFigures/diseaseRepresented/dMat_pheno.mat', ...
     'phenoDs')

load('~/resultsAndFigures/diseaseRepresented/dMat_pheno.mat')
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% now place the questions: for the OMIM
sib = unique(data{1});

omimAnnInds = ismember(data{1}, 'OMIM');

omimAnn.genes = data{3}(omimAnnInds);
omimAnn.DOID = data{7}(omimAnnInds);
omimAnn.DName = data{5}(omimAnnInds);
omimAnn.uniqueDOID = unique(omimAnn.DOID);

length(unique(omimAnn.DOID)) % 1389

dMat = zeros(18494, 1389);
for i = 1:length(omimAnn.genes)
    i
    [a, b] = ismember(omimAnn.genes(i), gpl570.uniqueSymbols);
    if (a)
        [aa bb] = ismember(omimAnn.DOID(i), omimAnn.uniqueDOID);
        dMat(b, bb) = 1;
    end
end

sib = sum(dMat);
hist(sib(sib > 0))

dCount = sum(dMat');

% genes which have pure links are more likely to be associated with
% a disease 

sumPure = pureNet{1};
for t = 2:5
    sumPure = sumPure + pureNet{t};
end

pnds = sum(sumPure) + sum(sumPure');
% I have 6658 pure genes. 

% total count of expressed genes
commonGenes  = sum(expMat') == 5;

% count of pure genes with a disease 
pureGenes = pnds > 0;

% count of non pure - common genes with a disease 
dCounts = sum(dMat');

sum(dCounts(commonGenes) > 0)/sum(commonGenes)
actualDL = sum(dCounts(pureGenes)> 0)/sum(pureGenes)

s% what is the p-value of that
fCounts = sum(fMat');
sum(fCounts(commonGenes) >0)/sum(commonGenes)
actualFL = sum(fCounts(pureGenes)>0)/sum(pureGenes)

corr(dCounts(pureGenes)', pnds(pureGenes)')

% compare it to random selection
commonGenesInd = find(commonGenes);
randomFL = zeros(1, 1000);
randomDL = zeros(1, 1000);
for i = 1:1000
    sampleGenes = datasample(commonGenesInd, sum(pureGenes), 'Replace', ...
                             false);
    randomDL(i) = sum(dCounts(sampleGenes)> 0)/length(sampleGenes);
    randomFL(i) = sum(fCounts(sampleGenes)>0)/length(sampleGenes);
end

% compare it to the genes without pure links
nonPureGenes = logical(commonGenes - pureGenes);
nonPureFL = sum(fCounts(nonPureGenes) > 0)/sum(nonPureGenes)
nonPureDL = sum(dCounts(nonPureGenes) > 0)/sum(nonPureGenes) 
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
qs = quantile(pnds(pnds > 0), 10)

sib = corr(pnds', dCount');
plot(pnds, dCount, '.')

% pure hub genes enriched with any disease in any of the networks


% are any diseases enriched in the network. 
% how is that enrichment affected by the pure links.

% 34. Load the new omim and do the same thing as network enrichment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
file1 = '/space/databases/omim/genemap.txt'
file2 = '/space/databases/omim/genemap2.txt'

linum = 4;
fid = fopen(file1)
dOne = textscan(fid, repmat('%s', 1, 13), 'headerlines', linum, ...
                  'Delimiter', '\t')

fid = fopen(file1)
temps = fgets(fid) % header is in line 4
headersOne = strsplit(temps, '\t')

% get the genes and disease associations
% geneSym: 6, disease: 12

% getting the count of annotations.
dCounts = 0;
max = 0;
ins = zeros(1, length(dOne{1}));
for i = 1:length(dOne{1})
    myS = dOne{12}{i};
    startInds = regexp(myS, '\d{6,6}');
    if startInds
        dCounts = dCounts + 1;
    end
    sib = length(startInds);
    ins(i) = sib;
    if sib > max
        max = sib;
    end
end

% getting the unique DIDs
dIDs = zeros(1, sum(ins));
c = 1;
for i = 1:length(dOne{1})
    myS = dOne{12}{i};
    startInds = regexp(myS, '\d{6,6}');
    if length(startInds)
        for j = 1:length(startInds)
            dIDs(c) = str2num(myS(startInds(j):startInds(j)+5));
            c = c + 1;
        end
    end
end

uniqueDIDs = unique(dIDs);
dMat = zeros(18494, length(uniqueDIDs));

% I have a list of unique diseases, I have a "list" of gene symbols
% associated with them - this doesn't allow me to do one time
% mapping, I have to parse each gene list separately. As far as I know, each of the gene symbols
% is associated with the disease. 

dataFolder = '~/data/general/'
load([dataFolder 'GPL570GemmaMapNEW.mat'])
for i = 1:length(dOne{1})
    thisDs = dOne{12}{i};
    startInds = regexp(thisDs, '\d{6,6}');
    if length(startInds)
        i
        % get the disease IDs
        myDIDs = zeros(1, length(startInds));
        for j = 1:length(startInds)
            myDIDs(j) = str2num(thisDs(startInds(j):(startInds(j)+5)));
        end
        % get the disease Indexes 
        [ad, bd] = ismember(myDIDs, uniqueDIDs);
        dids = bd(ad);
        
        % get the gene Syms
        myGenes = strsplit(dOne{6}{i}, ', ');
        [ag, bg] = ismember(myGenes, gpl570.uniqueSymbols);
        gids = bg(ag);
        
        for j = 1:length(gids)
            for k = 1:length(dids)
                dMat(gids(j), dids(k)) = 1;
            end
        end
    end
end

tempIDs = find(sum(dMat));
d.Mat = dMat(:, tempIDs);
d.IDs = uniqueDIDs(tempIDs);
% get the disease name and add it
file = '/space/databases/omim/mimTitles.txt'
linum = 3;
fid = fopen(file)
titles = textscan(fid, ['%s%d' repmat('%s', 1, 11)], 'headerlines', linum, ...
                  'Delimiter', '\t')

titleMap = containers.Map(titles{2}, titles{3});

d.titleMap = titleMap;

save('~/resultsAndFigures/diseaseRepresented/dMat.mat', ...
     'd')

% Ok I got the disease Mat! I have 4804 diseases and 3684 genes and
% 5604 annotations.
% now that enrichment thing:

% 35. the same function enrichment for disease
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I got the results from phenocarta

clear

%load('~/resultsAndFigures/diseaseRepresented/dMat.mat')
load('~/resultsAndFigures/diseaseRepresented/dMat_pheno.mat')

FC = '3'
FDR = '0012'
load(sprintf(['~/resultsAndFigures/TSlinks/finalTables/' ...
              'finalTable_CG13_FC%slog_FDR%s.mat'], FC, FDR))

% >>>> load the five ATN and TSN networks
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

tsnLC = sum(tsnds')./2;

expMat = zeros(18494, 5);
for t = 1:5
    tissue = tissues{t};
    load(['~/data/general/tissueExpGenes/' tissue ...
          'ExpGenes0.8.mat'])
    expMat(:, t) = expGenesInd;
end

load('~/resultsAndFigures/TSgenes.mat')
for t = 1:5
    tsNet = finalTable(t).wholeNet;
    mytgenes = TSGenes.GeneList(TSGenes.tissue == t);
    tsNet(:, mytgenes) = 0;
    tsNet(mytgenes, :) = 0;
    tsNetNotTSgenes{t} = tsNet;
end

% 2. remove the ExpInduced links
load('~/resultsAndFigures/expInOneTissueGenes.mat')
for t = 1:5
    tsNet = finalTable(t).wholeNet;
    tsNet(:, logical(TSExpMarkedGenes(:, t))) = 0;
    tsNet(logical(TSExpMarkedGenes(:, t)), :) = 0;
    tsNetNonExpInd{t} = tsNet;
end

% 3. remove the pure link
for t = 1:5
    tsNet = finalTable(t).wholeNet;
    pureNet = (finalTable(t).cgNet + finalTable(t).noTSNet) == 2;
    tsNetNonPure{t} = tsNet - pureNet;
end

dCount = length(phenoDs.DOIDs)
gc = zeros(5, dCount);
bgExp = zeros(5, dCount);
obCount = zeros(5, dCount);
LRs = zeros(5, dCount);
hgps = zeros(5, dCount);
binps = zeros(5, dCount);
netDensities = [.007, .020, .008, .012, .006];
netDensities = [.0004, .004, .0011, .0005,.0004];

% each of this is one round of sampling...
dMat = phenoDs.mat > 0 + 0;
for t = 1:5
    %  myNet = finalTable(t).wholeNet;
    % myNet = atn{t};
     myNet = tsNetNonPure{t};
    %myNet = tsNetNonExpInd{t};
    
    pureNet = (finalTable(t).cgNet + finalTable(t).noTSNet) == 2;
    tsNet = finalTable(t).wholeNet;

    for f = 1:dCount
        myTerm = dMat(:, f);
        sum(myTerm)

        % how many genes are marked as expressed in this term?
        eGenesCount = sum((myTerm + expMat(:, t)) == 2)
        eGenes = ((myTerm + expMat(:, t)) == 2);
        gc(t, f) = eGenesCount;
        % for ppin >>>>>>>>>>>
        % eGenesCount = sum((myTerm + expArr') == 2);
        % eGenes = ((myTerm + expArr') == 2);
        % <<<<<<<<<<<<<<<<<<<<

        if(eGenesCount >= 5)

            bgExp(t, f) = (eGenesCount * (eGenesCount-1)/2) * netDensities(t);
            litNet = myNet(eGenes, eGenes);
            obCount(t, f) = sum(sum(litNet));
            LRs(t, f) = (sum(sum(litNet))/((eGenesCount * (eGenesCount-1)/2) ...
                                           * netDensities(t) + 1));
            
            %            tGenes = sum(expMat(:, tissueID));
            N = (eGenesCount * (eGenesCount-1)/2);
            %            M = tGenes * (tGenes - 1) / 2;
            %            k = M * netDensity;
            x = sum(sum(litNet));
            %        hgps(f) = 1 - hygecdf(x, M, k,N);
            binps(t, f) = 1- binocdf(x, N, netDensities(t));
            
            pNetLC(t, f) = sum(sum(pureNet(eGenes, eGenes)));
            tsnLC(t, f) = sum(sum(tsNet(eGenes, eGenes)));
        else 
            bgExp(t, f) = -1;
            obCount(t, f) = -1;
            LRs(t, f) = -1;
            binps(t, f) = -1;
            pNetLC(t, f) = -1;
            tsnLC(t, f) = -1;
        end
        % so this function is not adding to that expectation, but which one
        % is? 
    end
end

% selected Fs for each of the tissues: 
sigTerms = zeros(5, dCount);
for t = 1:5
    selected = binps(t, :) > -1;
    sum(selected)
    sLRs = LRs(t, selected);
    sbinps = binps(t, selected);
    sterms = find(selected);

    % correcting pvalues
    sbinps = sbinps * sum(selected);
    logps = log10(sbinps);
    logps(logps < -10) = -10;
    logps = logps .* -1;
    %    plot(logps, '.')
    temp1 = (logps >= 2);
    sigTerms(t, sterms(temp1)) = 1;
end

tsnNonExpIndDiseaseRes.bgExp = bgExp;
tsnNonExpIndDiseaseRes.obCount = obCount;
tsnNonExpIndDiseaseRes.LRs = LRs;
tsnNonExpIndDiseaseRes.hgps = hgps;
tsnNonExpIndDiseaseRes.binps = binps;
tsnNonExpIndDiseaseRes.sigTerms = sigTerms;
tsnNonExpIndDiseaseRes.pNetLC = pNetLC;
tsnNonExpIndDiseaseRes.tsnLC = tsnLC;
tsnNonExpIndDiseaseRes.termIDs = phenoDs.DOIDs;
tsnNonExpIndDiseaseRes.gc = gc;

save('~/resultsAndFigures/diseaseRepresented/tsnNonEIDiseaseRes_pheno.mat', ...
     'tsnNonExpIndDiseaseRes')

tsnNonPureDiseaseRes.bgExp = bgExp;
tsnNonPureDiseaseRes.obCount = obCount;
tsnNonPureDiseaseRes.LRs = LRs;
tsnNonPureDiseaseRes.hgps = hgps;
tsnNonPureDiseaseRes.binps = binps;
tsnNonPureDiseaseRes.sigTerms = sigTerms;
tsnNonPureDiseaseRes.pNetLC = pNetLC;
tsnNonPureDiseaseRes.tsnLC = tsnLC;
tsnNonPureDiseaseRes.termIDs = phenoDs.DOIDs;
tsnNonPureDiseaseRes.gc = gc;

save('~/resultsAndFigures/diseaseRepresented/tsnNonPureDiseaseRes_pheno.mat', ...
     'tsnNonPureDiseaseRes')

tsnDiseaseRes.bgExp = bgExp;
tsnDiseaseRes.obCount = obCount;
tsnDiseaseRes.LRs = LRs;
tsnDiseaseRes.hgps = hgps;
tsnDiseaseRes.binps = binps;
tsnDiseaseRes.sigTerms = sigTerms;
tsnDiseaseRes.pNetLC = pNetLC;
tsnDiseaseRes.tsnLC = tsnLC;
tsnDiseaseRes.termIDs = phenoDs.DOIDs;
tsnDiseaseRes.gc = gc;

save('~/resultsAndFigures/diseaseRepresented/tsnDiseaseRes_pheno.mat', ...
     'tsnDiseaseRes')
load(['~/resultsAndFigures/diseaseRepresented/' ...
      'tsnDiseaseRes_pheno.mat'])

tanDiseaseRes.bgExp = bgExp;
tanDiseaseRes.obCount = obCount;
tanDiseaseRes.LRs = LRs;
tanDiseaseRes.hgps = hgps;
tanDiseaseRes.binps = binps;
tanDiseaseRes.pNetLC = pNetLC;
tanDiseaseRes.tsnLC = tsnLC;
tanDiseaseRes.sigTerms = sigTerms;
tanDiseaseRes.termIDs = phenoDs.DOIDs;
tanDiseaseRes.gc = gc;

save('~/resultsAndFigures/diseaseRepresented/tanDiseaseRes_pheno.mat', ...
     'tanDiseaseRes')
load(['~/resultsAndFigures/diseaseRepresented/' ...
      'tanDiseaseRes_pheno.mat'])

save('~/resultsAndFigures/diseaseRepresented/tanDiseaseRes_pheno.mat', ...
     'tanDiseaseRes')


% to get the ts exclusive terms
tsExpSigTerms = zeros(5, dCount);
for t = 1:5
    tsst = tsnDiseaseRes.sigTerms;
    tast = tanDiseaseRes.sigTerms;
    tsst(t, :) = [];
    tast(t, :) = [];
    book = tsnDiseaseRes.sigTerms(t, :);
    sib = sum(tsst) + sum(tast);
    final = (book - sib) > 0;
    tsExpSigTerms(t, :) = final;
end

% only 4 terms are exclusively enriched: 
% I want to print:
% termID
% term title, 
% term presence in TAN, presence in TSN, TAN links, TSN Links, pure links

clear
load('~/resultsAndFigures/diseaseRepresented/tsnDiseaseRes_pheno.mat')
load(['~/resultsAndFigures/diseaseRepresented/' ...
      'tanDiseaseRes_pheno.mat'])

file = sprintf(['~/resultsAndFigures/diseaseRepresented/' ...
                'enrichedDiseaseTerms_pheno_TAN_TSN_table.tsv'])

fid = fopen(file, 'w')
fprintf(fid, ['DOID\t'...
%'OMIM_ID\tOMIM_title\t' ...
              'blood_TANenriched\tblood_TSNenriched\t'...
              'blood_TAN_linkCount\tblood_TSN_linkCount\tblood_pure_linkCount\t' ...
              'brain_TANenriched\tbrain_TSNenriched\t'...
              'brain_TAN_linkCount\tbrain_TSN_linkCount\tbrain_pure_linkCount\t' ...
              'liver_TANenriched\tliver_TSNenriched\t'...
              'liver_TAN_linkCount\tliver_TSN_linkCount\tliver_pure_linkCount\t' ...
              'lung_TANenriched\tlung_TSNenriched\t'...
              'lung_TAN_linkCount\tlung_TSN_linkCount\tlung_pure_linkCount\t' ...
              'SMuscle_TANenriched\tSMuscle_TSNenriched\t'...
              'SMuscle_TAN_linkCount\tSMuscle_TSN_linkCount\tSMuscle_pure_linkCount\t'])


dsInd = tanDiseaseRes.sigTerms + tsnDiseaseRes.sigTerms;
[a, b, c] = find(sum(dsInd));
IDs = phenoDs.DOIDs(b)

IDs = d.IDs(b);
titles = cell(1, length(IDs));

for i = 1:length(titles)
    titles{i} = titleMap(IDs(i));
end

for i = 1:length(IDs)
    %    fprintf(fid, ['%d\t%s\t' repmat('%d\t', 1, 25) '\n'], ...
        fprintf(fid, ['%d\t' repmat('%d\t', 1, 25) '\n'], ...
            IDs(i),...  
        %            titles{i},...
            tanDiseaseRes.sigTerms(1, b(i)),...
            tsnDiseaseRes.sigTerms(1, b(i)),...
            tanDiseaseRes.obCount(1, b(i)),...
            tsnDiseaseRes.obCount(1, b(i)),...
            tsnDiseaseRes.pNetLC(1, b(i)),...
            
            tanDiseaseRes.sigTerms(2, b(i)),...
            tsnDiseaseRes.sigTerms(2, b(i)),...
            tanDiseaseRes.obCount(2, b(i)),...
            tsnDiseaseRes.obCount(2, b(i)),...
            tsnDiseaseRes.pNetLC(2, b(i)),...
            
            tanDiseaseRes.sigTerms(3, b(i)),...
            tsnDiseaseRes.sigTerms(3, b(i)),...
            tanDiseaseRes.obCount(3, b(i)),...
            tsnDiseaseRes.obCount(3, b(i)),...
            tsnDiseaseRes.pNetLC(3, b(i)),...
            
            tanDiseaseRes.sigTerms(4, b(i)),...
            tsnDiseaseRes.sigTerms(4, b(i)),...
            tanDiseaseRes.obCount(4, b(i)),...
            tsnDiseaseRes.obCount(4, b(i)),...
            tsnDiseaseRes.pNetLC(4, b(i)),...
            
            tanDiseaseRes.sigTerms(5, b(i)),...
            tsnDiseaseRes.sigTerms(5, b(i)),...
            tanDiseaseRes.obCount(5, b(i)),...
            tsnDiseaseRes.obCount(5, b(i)),...
            tsnDiseaseRes.pNetLC(5, b(i)))
end

% 36. everything for the disease stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

% >>>>> loading essentials
load(['~/resultsAndFigures/diseaseRepresented/' ...
      'tsnNonPureDiseaseRes_pheno.mat'])
load('~/resultsAndFigures/diseaseRepresented/tsnNonEIDiseaseRes_pheno.mat')
load('~/resultsAndFigures/diseaseRepresented/tanDiseaseRes_pheno.mat')
load('~/resultsAndFigures/diseaseRepresented/tsnDiseaseRes_pheno.mat')

% to get the ts exclusive terms
tsExpSigTerms = zeros(5, 2289);
for t = 1:5
    tsst = tsnDiseaseRes.sigTerms;
    tast = tanDiseaseRes.sigTerms;
    tsst(t, :) = [];
    tast(t, :) = [];
    book = tsnDiseaseRes.sigTerms(t, :);
    sib = sum(tsst) + sum(tast);
    final = (book - sib) > 0;
    tsExpSigTerms(t, :) = final;
end

% >>>> getting the floss and indFloss for three sammple groups:
% >>> 1. sampleEIlink_countOfPure
load(['~/resultsAndFigures/diseaseRepresented/' ...
      'sampleEIlink_countOfPure_RemainSigTerms_pheno.mat'])

floss = zeros(100, 5);
for i = 1:100
    book = sampleRemainSigTerms{i};
    kado = (tsExpSigTerms - book) >0;
    floss(i, :) = sum((kado>0)');
end

save(['~/resultsAndFigures/diseaseRepresented/' ...
      'affectedTerms_sampleEIlink_countOfPure_pheno.mat', 'floss'])

indFloss = zeros(5, 2289);
for i = 1:100
    book = sampleRemainSigTerms{i};
    kado = (tsExpSigTerms .* book);
    indFloss = indFloss + kado;
end

save('~/resultsAndFigures/functionsRepresented/timesTermsMissed_sampleEIlink_countOfPure_pheno.mat', ...
     'indFloss')

eiPureCountIndFloss = indFloss;
clear sampleRemainSigTerms;
clear indFloss

% >>> 2. sampleEIlink
load(['~/resultsAndFigures/diseaseRepresented/' ...
      'sampleEIlinkRemainSigTerms_pheno.mat'])

floss = zeros(100, 5);
for i = 1:100
    book = sampleRemainSigTerms{i};
    kado = (tsExpSigTerms - book) >0;
    floss(i, :) = sum((kado>0)');
end

save(['~/resultsAndFigures/diseaseRepresented/' ...
      'affectedTerms_sampleEIlink_pheno.mat', 'floss'])

indFloss = zeros(5, 2289);
for i = 1:100
    book = sampleRemainSigTerms{i};
    kado = (tsExpSigTerms .* book);
    indFloss = indFloss + kado;
end

save('~/resultsAndFigures/functionsRepresented/timesTermsMissed_sampleEIlink_pheno.mat', ...
     'indFloss')

eiRandomSigTerms = sampleRemainSigTerms;
clear sampleRemainSigTerms;
eiIndFloss = indFloss;
clear indFloss

% >>>>> 3. samplePureLink
load(['~/resultsAndFigures/diseaseRepresented/' ...
      'samplePureRemainSigTerms_Pheno.mat'])

floss = zeros(100, 5);
for i = 1:100
    book = sampleRemainSigTerms{i};
    kado = (tsExpSigTerms - book) >0;
    floss(i, :) = sum((kado>0)');
end

save(['~/resultsAndFigures/diseaseRepresented/' ...
      'affectedTerms_samplePureRemainSigTerms_pheno.mat', 'floss'])

indFloss = zeros(5, 2289);
for i = 1:100
    book = sampleRemainSigTerms{i};
    kado = (tsExpSigTerms .* book);
    indFloss = indFloss + kado;
end

save('~/resultsAndFigures/functionsRepresented/timesTermsMissed_samplePureRemainSigTerms_pheno.mat', ...
     'indFloss')

pureIndFloss = indFloss;
clear indFloss
pureRandomSigTerms = sampleRemainSigTerms;
clear randomSigTerms

load('~/resultsAndFigures/diseaseRepresented/dMat_pheno.mat')

tissues = {'blood', 'brain', 'liver', 'lung', 'skeletalMuscle'};
totalContributingEICount = zeros(1, 5);
totalContributingPureCount = zeros(1, 5);
for t = 2:4
    file = sprintf(['~/resultsAndFigures/diseaseRepresented/' ...
    '%sTSNExclusiveTerms_table_withLog10Pvalues_withCounts_withEIPureCountRemoval_Ds.txt'], tissues{t})
    fid = fopen(file, 'w')
    fprintf(fid, ['DOID\tName\tgeneCount\' ...
                  'tlinkCount(lc)\tlc_minusPurLinks\' ...
                  'tlc_minusExpressionInducedLinks\t1-ratioOfPureLinks' ...
                  '\t1-RatioOfExpressionInducedLinks\' ...
                  'tSigAfterPureRemoval\' ...
                  'tSigAfterExpressionInducedRemoval\t'...
                  'log10Pvalue_plus1e-20\t'...
                  'log10PvalueAfterPureRemoval_plus1e-20\t'...
                  'log10PvalueAfterEIremovalplus1e-20\t'...
                  'significantCountIn100_randomPureLinkRemoval\t'...
    %                  'significantCountIn100_randomEILinkRemoval_pureCount\t'...
                  'significantCountIn100_randomEILinkRemoval\n'])
    termCount = sum(tsExpSigTerms(t, :));
    terms = find(tsExpSigTerms(t, :));
    for i = 1:termCount
        thisTerm = terms(i);
        nonPureObRatio = tsnNonPureDiseaseRes.obCount(t, thisTerm) / ...
            tsnDiseaseRes.obCount(t, thisTerm);
        nonTsObRatio = tsnNonExpIndDiseaseRes.obCount(t, thisTerm) / ...
            tsnDiseaseRes.obCount(t, thisTerm);
        
        totalContributingEICount(t) = totalContributingEICount(t) + ...
            (tsnDiseaseRes.obCount(t, thisTerm) - tsnNonExpIndDiseaseRes.obCount(t, thisTerm));
        totalContributingPureCount(t) = totalContributingPureCount(t) + ...
                    (tsnDiseaseRes.obCount(t, thisTerm) - tsnNonPureDiseaseRes.obCount(t, ...
                                                          thisTerm));
        
        DOID = tsnDiseaseRes.termIDs(thisTerm);
        term = phenoDs.map(DOID)
        % for the first few I can use tsRes or tsNonPureRes or
        % tsNonExpIndRes. they have the same values
        fprintf(fid, ['%d\t%s\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%d\t'...
                      '%d\t%.2f\t%.2f\t%.2f\t%d\t',...
        %                      '%d\t',...
                      '%d\n'], ...
                tsnDiseaseRes.termIDs(thisTerm),...  
                term,...
                tsnDiseaseRes.gc(t, thisTerm),...
                tsnDiseaseRes.obCount(t, thisTerm),...
                tsnNonPureDiseaseRes.obCount(t, thisTerm),...
                tsnNonExpIndDiseaseRes.obCount(t, thisTerm),...
                nonPureObRatio,...
                nonTsObRatio,...
                tsnNonPureDiseaseRes.sigTerms(t, thisTerm),...
                tsnNonExpIndDiseaseRes.sigTerms(t, thisTerm),...
                log10(tsnDiseaseRes.binps(t, thisTerm) + 1e-20),...
                log10(tsnNonPureDiseaseRes.binps(t, thisTerm) + 1e-20),...
                log10(tsnNonExpIndDiseaseRes.binps(t, thisTerm) + 1e-20),...
                pureIndFloss(t, thisTerm),...
        %                eiPureCountIndFloss(t, thisTerm),...
                eiIndFloss(t, thisTerm));
    end
    fclose(fid)
end

% 37. pure genes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nds = zeros(18494, 5);
for t = 1:5
    myNet = pureNet{t};
    ndplus = sum( myNet + myNet') > 0;
    nds(:, t) = ndplus;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dOneCount = length(unique(dOne{12}))
uniqueDOnes = unique(dOne{12});
uniqueDOnes(1) = [];
d1Mat = zeros(18494, dOneCount);
[a, b] = ismember(dOne{12}, uniqueDOnes);

linum = 4;
fid = fopen(file2)
dTwo = textscan(fid, repmat('%s', 1, 14), 'headerlines', linum, ...
                  'Delimiter', '\t')

fid = fopen(file2)
temps = fgets(fid) % header is in line 4
headersTwo = strsplit(temps, '\t')

% get the genes and disease associations
% geneSym: 7, disease: 13
dTwoCount = length(unique(dTwo{13}))
uniqueDTwos = unique(dTwo{13});
uniqueDTwos(1) = [];
d2Mat = zeros(18494, dTwoCount);
[a, b] = ismember(dTwo{13}, uniqueDTwos);

% extract pheynotypes associated with each gene: 
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
phenoLength = 0;

% get the annotations from the first file. do the analysis, report
% it, finalize the file. 

